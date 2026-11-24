use crate::deseq::{ALPHA, GeneSummary, Significance};
use crate::{Args, Sample};
use anyhow::{Context, Result, ensure};
use nalgebra::DMatrix;
use resvg::tiny_skia::{Pixmap, Transform};
use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::BufWriter;
use std::path::Path;
use svg2pdf::{ConversionOptions, PageOptions, usvg};

const INK: &str = "#262626";
const GRID: &str = "#e8e8e8";
const VOLCANO: [(Significance, &str); 5] = [
    (Significance::Down, "#a6cee3"),
    (Significance::NonSignificant, "#999999"),
    (Significance::SignificantDown, "#377eb8"),
    (Significance::SignificantUp, "#e41a1c"),
    (Significance::Up, "#fb9a99"),
];
const GROUP_COLORS: [&str; 8] = [
    "#0072b2", "#d55e00", "#009e73", "#cc79a7", "#e69f00", "#56b4e9", "#7c5ba7", "#545454",
];

pub fn write(
    out: &Path,
    tag: &str,
    args: &Args,
    samples: &[Sample],
    vst: &[f64],
    rows: &[GeneSummary],
) -> Result<()> {
    let sample_count = samples.len();
    let genes = vst.len() / sample_count;
    ensure!(genes >= 2, "PCA needs at least two genes after filtering");
    let mut fonts = usvg::Options::default();
    let fontdb = fonts.fontdb_mut();
    fontdb.load_system_fonts();
    if let Some(prefix) = std::env::var_os("CONDA_PREFIX") {
        fontdb.load_fonts_dir(Path::new(&prefix).join("share/fonts"));
    }
    let font = [
        "Arial",
        "Liberation Sans",
        "DejaVu Sans",
        "DejaVu Sans Mono",
    ]
    .into_iter()
    .find(|name| {
        fontdb
            .faces()
            .any(|face| face.families.iter().any(|f| f.0 == *name))
    })
    .map(str::to_owned)
    .or_else(|| {
        fontdb
            .faces()
            .find_map(|face| face.families.first().map(|f| f.0.clone()))
    })
    .context(
        "No font available for the plots. Install a system font or activate the Conda environment.",
    )?;
    fontdb.set_sans_serif_family(&font);
    fonts.font_family = font;

    let means: Vec<_> = vst
        .chunks_exact(sample_count)
        .map(|row| row.iter().sum::<f64>() / sample_count as f64)
        .collect();
    let centered = DMatrix::from_fn(sample_count, genes, |sample, gene| {
        vst[gene * sample_count + sample] - means[gene]
    });
    let svd = centered.svd(true, false);
    let u = svd.u.context("PCA scores could not be calculated")?;
    let d = svd.singular_values;
    let variance = d.iter().map(|value| value * value).sum::<f64>();
    let explained = |component: usize| {
        if variance > 0.0 {
            d[component].powi(2) * 100.0 / variance
        } else {
            0.0
        }
    };
    let scores: Vec<_> = (0..sample_count)
        .map(|index| (u[(index, 0)] * d[0], u[(index, 1)] * d[1]))
        .collect();
    pca(samples, &scores, (explained(0), explained(1))).save_pdf(&out.join("PCA.pdf"), &fonts)?;

    let sample_means: Vec<_> = (0..sample_count)
        .map(|i| {
            vst.chunks_exact(sample_count)
                .map(|row| row[i])
                .sum::<f64>()
                / genes as f64
        })
        .collect();
    let mut cross = vec![vec![0.0; sample_count]; sample_count];
    for row in vst.chunks_exact(sample_count) {
        for (i, (&value, mean)) in row.iter().zip(&sample_means).enumerate() {
            let delta = value - mean;
            for j in 0..=i {
                cross[i][j] += delta * (row[j] - sample_means[j]);
            }
        }
    }
    let correlation: Vec<Vec<_>> = (0..sample_count)
        .map(|i| {
            (0..sample_count)
                .map(|j| {
                    (cross[i.max(j)][i.min(j)] / (cross[i][i] * cross[j][j]).sqrt())
                        .clamp(-1.0, 1.0)
                })
                .collect()
        })
        .collect();
    let mut csv = csv::Writer::from_path(out.join("sample_correlation.csv"))?;
    csv.write_record(std::iter::once("").chain(samples.iter().map(|sample| sample.id.as_str())))?;
    for (sample, row) in samples.iter().zip(&correlation) {
        csv.serialize((
            &sample.id,
            row.iter()
                .map(|v| {
                    if v.is_finite() {
                        format!("{v:e}")
                    } else {
                        "NA".to_owned()
                    }
                })
                .collect::<Vec<_>>(),
        ))?;
    }
    csv.flush()?;
    heatmap(samples, &correlation).save_pdf(&out.join("sample_correlation.pdf"), &fonts)?;
    ma(rows).save_pdf(&out.join("MA_plot.pdf"), &fonts)?;
    let plot = volcano(args, rows);
    let tree = plot.tree(&fonts)?;
    save_pdf(&tree, &out.join(format!("volcano_DESeq2_{tag}.pdf")))?;
    save_png(&tree, &out.join(format!("volcano_DESeq2_{tag}.png")))?;
    Ok(())
}

fn pca(samples: &[Sample], scores: &[(f64, f64)], explained: (f64, f64)) -> Svg {
    let mut svg = Svg::new(432, 360);
    let chart = Chart {
        left: 62.0,
        top: 48.0,
        width: 258.0,
        height: 246.0,
        x: padded(scores.iter().map(|p| p.0), 0.13),
        y: padded(scores.iter().map(|p| p.1), 0.13),
    };
    chart.draw(
        &mut svg,
        "VST PCA",
        &format!("PC1 ({:.1}%)", explained.0),
        &format!("PC2 ({:.1}%)", explained.1),
        false,
    );
    let mut levels: Vec<_> = samples
        .iter()
        .map(|sample| sample.condition.as_str())
        .collect();
    levels.sort_unstable();
    levels.dedup();
    for (index, condition) in levels.iter().enumerate() {
        let color = GROUP_COLORS[index % GROUP_COLORS.len()];
        let row = 58.0 + index as f64 * 13.0;
        svg.point(334.0, row, 2.6, color, 1.0);
        svg.text(342.0, row + 3.0, condition, 9.0, "start", INK);
    }
    let mut labels = Vec::new();
    for (sample, &(x, y)) in samples.iter().zip(scores) {
        let (px, py) = chart.point(x, y);
        let index = levels
            .binary_search(&sample.condition.as_str())
            .unwrap_or_default();
        let color = GROUP_COLORS[index % GROUP_COLORS.len()];
        svg.point(px, py, 3.2, color, 1.0);
        label(
            &mut svg,
            &chart,
            &mut labels,
            (px, py),
            &sample.id,
            INK,
            8.0,
        );
    }
    svg
}

fn ma(rows: &[GeneSummary]) -> Svg {
    let mut svg = Svg::new(432, 360);
    let chart = Chart {
        left: 70.0,
        top: 48.0,
        width: 331.0,
        height: 246.0,
        x: padded(
            rows.iter()
                .filter(|g| g.base_mean > 0.0)
                .map(|g| g.base_mean.log10()),
            0.04,
        ),
        y: (-8.0, 8.0),
    };
    chart.draw(
        &mut svg,
        "MA plot",
        "Mean of normalized counts",
        "Log2 fold change",
        true,
    );
    let (_, zero) = chart.point(chart.x.0, 0.0);
    svg.line((chart.left, zero), (chart.right(), zero), "#636363", 0.85);
    for row in rows {
        let Some(lfc) = row.lfc.filter(|_| row.base_mean > 0.0) else {
            continue;
        };
        let color = if row.padj.is_some_and(|p| p < ALPHA) {
            "#377eb8"
        } else {
            "#8f8f8f"
        };
        let (x, y) = chart.point(row.base_mean.log10(), lfc.clamp(chart.y.0, chart.y.1));
        if lfc.abs() > chart.y.1 {
            let tip = if lfc > 0.0 { y } else { y - 3.5 };
            let base = if lfc > 0.0 { y + 3.5 } else { y };
            let _ = write!(
                svg.body,
                r#"<path d="M{x} {tip} l-2 {dy} h4 Z" fill="{color}"/>"#,
                dy = base - tip
            );
        } else {
            svg.point(x, y, 1.65, color, 0.82);
        }
    }
    svg
}

fn volcano(args: &Args, rows: &[GeneSummary]) -> Svg {
    let mut svg = Svg::new(936, 648);
    let y_max = rows
        .iter()
        .map(|gene| gene.neg_log_padj)
        .fold(0.0, f64::max);
    let chart = Chart {
        left: 86.0,
        top: 76.0,
        width: 610.0,
        height: 490.0,
        x: padded(rows.iter().filter_map(|gene| gene.lfc), 0.07),
        y: (0.0, (y_max * 1.08).max(1.0)),
    };
    let title = format!(
        "Volcano Plot of Gene Expression Differences ({} vs {})",
        args.condition_a, args.condition_b
    );
    chart.draw(
        &mut svg,
        &title,
        "Log2 fold change",
        "-Log10 adjusted p value",
        false,
    );
    svg.text(
        720.0,
        253.0,
        "Differential Significance",
        14.0,
        "start",
        INK,
    );
    for (index, (group, color)) in VOLCANO.iter().enumerate() {
        let y = 278.0 + index as f64 * 21.0;
        svg.point(727.0, y - 3.0, 3.6, color, 1.0);
        svg.text(739.0, y, group.label(), 12.0, "start", INK);
    }
    for row in rows.iter().rev() {
        if let Some(lfc) = row.lfc {
            let (x, y) = chart.point(lfc, row.neg_log_padj);
            let color = VOLCANO
                .iter()
                .find(|(kind, _)| *kind == row.significance)
                .map_or("#999999", |(_, c)| c);
            svg.point(x, y, 2.6, color, 0.9);
        }
    }
    let mut labels = Vec::new();
    for row in rows.iter().filter(|row| row.show_label) {
        if let Some(lfc) = row.lfc {
            label(
                &mut svg,
                &chart,
                &mut labels,
                chart.point(lfc, row.neg_log_padj),
                &row.label,
                "#b21420",
                9.0,
            );
        }
    }
    svg
}

fn heatmap(samples: &[Sample], values: &[Vec<f64>]) -> Svg {
    let mut svg = Svg::new(432, 360);
    svg.text(216.0, 25.0, "Sample correlation (VST)", 15.0, "middle", INK);
    let (order, merges) = cluster(values);
    let extent = 196.0;
    let cell = extent / samples.len() as f64;
    let size = (0.56 * cell).clamp(3.0, 10.0);
    let left = 83.0;
    let top = 89.0;
    for (i, &row) in order.iter().enumerate() {
        let y = top + i as f64 * cell;
        svg.text(
            left + extent + 5.0,
            y + cell / 2.0 + size * 0.32,
            &samples[row].id,
            size,
            "start",
            INK,
        );
        for (j, &column) in order.iter().enumerate() {
            let x = left + j as f64 * cell;
            svg.rect((x, y), (cell, cell), &heat_color(values[row][column]));
        }
        let x = left + (i as f64 + 0.5) * cell;
        let y = top + extent + 5.0;
        let _ = write!(svg.body, r#"<g transform="rotate(-90 {x} {y})">"#);
        svg.text(x, y + size * 0.34, &samples[row].id, size, "end", INK);
        svg.body.push_str("</g>");
    }
    let mut centers = vec![(0.0, 0.0); samples.len()];
    for (index, &leaf) in order.iter().enumerate() {
        centers[leaf].0 = (index as f64 + 0.5) * cell;
    }
    let height = merges.last().map_or(1.0, |m| m.2).max(f64::EPSILON);
    for (first, second, distance) in merges {
        let (a, a_depth) = centers[first];
        let (b, b_depth) = centers[second];
        let depth = distance.sqrt() / height.sqrt() * 39.0;
        for (start, end) in [
            ((a, a_depth), (a, depth)),
            ((a, depth), (b, depth)),
            ((b, depth), (b, b_depth)),
        ] {
            svg.line(
                (left + start.0, top - start.1),
                (left + end.0, top - end.1),
                INK,
                0.7,
            );
            svg.line(
                (left - start.1, top + start.0),
                (left - end.1, top + end.0),
                INK,
                0.7,
            );
        }
        centers.push(((a + b) / 2.0, depth));
    }
    for i in 0..100 {
        let y = top + extent / 100.0 * i as f64;
        svg.rect(
            (357.0, y),
            (10.0, extent / 100.0 + 0.02),
            &heat_color(1.0 - 2.0 * f64::from(i) / 99.0),
        );
    }
    for (offset, value) in [(0.0, "1"), (extent / 2.0, "0"), (extent, "-1")] {
        svg.text(372.0, top + offset + 3.0, value, 9.0, "start", INK);
    }
    svg
}

type Merge = (usize, usize, f64);

fn cluster(values: &[Vec<f64>]) -> (Vec<usize>, Vec<Merge>) {
    let n = values.len();
    let mut distances = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..i {
            let distance = values[i]
                .iter()
                .zip(&values[j])
                .filter(|(a, b)| a.is_finite() && b.is_finite())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>();
            distances[i][j] = distance;
            distances[j][i] = distance;
        }
    }
    let mut clusters: Vec<_> = (0..n).map(|i| vec![i]).collect();
    let mut active: Vec<_> = (0..n).collect();
    let mut nodes = active.clone();
    let mut merges = Vec::with_capacity(n - 1);
    while active.len() > 1 {
        let (i, j) = active
            .iter()
            .enumerate()
            .flat_map(|(index, &i)| active[index + 1..].iter().map(move |&j| (i, j)))
            .min_by(|&(i, j), &(k, l)| distances[i][j].total_cmp(&distances[k][l]))
            .expect("at least two active clusters");
        let other = std::mem::take(&mut clusters[j]);
        clusters[i].extend(other);
        merges.push((nodes[i], nodes[j], distances[i][j]));
        nodes[i] = n + merges.len() - 1;
        active.retain(|&index| index != j);
        for &k in &active {
            let distance = distances[i][k].max(distances[j][k]);
            distances[i][k] = distance;
            distances[k][i] = distance;
        }
    }
    (std::mem::take(&mut clusters[active[0]]), merges)
}

fn heat_color(value: f64) -> String {
    if !value.is_finite() {
        return "#dedede".to_owned();
    }
    let endpoint = if value >= 0.0 {
        [178, 24, 43]
    } else {
        [33, 102, 172]
    };
    let channel = |i: usize| (255.0 + value.abs() * (f64::from(endpoint[i]) - 255.0)).round() as u8;
    format!("#{:02x}{:02x}{:02x}", channel(0), channel(1), channel(2))
}

struct Chart {
    left: f64,
    top: f64,
    width: f64,
    height: f64,
    x: (f64, f64),
    y: (f64, f64),
}

impl Chart {
    fn right(&self) -> f64 {
        self.left + self.width
    }
    fn bottom(&self) -> f64 {
        self.top + self.height
    }
    fn point(&self, x: f64, y: f64) -> (f64, f64) {
        (
            self.left + self.width * (x - self.x.0) / (self.x.1 - self.x.0),
            self.bottom() - self.height * (y - self.y.0) / (self.y.1 - self.y.0),
        )
    }

    fn draw(&self, svg: &mut Svg, title: &str, x_label: &str, y_label: &str, log_x: bool) {
        let x_center = self.left + self.width / 2.0;
        let y_center = self.top + self.height / 2.0;
        let large = self.width > 500.0;
        let scale = if large { 1.35 } else { 1.0 };
        svg.text(
            svg.width as f64 / 2.0,
            self.top - 21.0,
            title,
            14.0 * scale,
            "middle",
            INK,
        );
        for x in ticks(self.x, 5.0) {
            let px = self.point(x, self.y.0).0;
            svg.line((px, self.top), (px, self.bottom()), GRID, 0.6);
            let label = if log_x {
                tick_label(10f64.powf(x))
            } else {
                tick_label(x)
            };
            svg.text(px, self.bottom() + 15.0, &label, 9.0 * scale, "middle", INK);
        }
        for y in ticks(self.y, 5.0) {
            let py = self.point(self.x.0, y).1;
            svg.line((self.left, py), (self.right(), py), GRID, 0.6);
            svg.text(
                self.left - 7.0,
                py + 3.0,
                &tick_label(y),
                9.0 * scale,
                "end",
                INK,
            );
        }
        svg.line((self.left, self.top), (self.left, self.bottom()), INK, 0.75);
        svg.line(
            (self.left, self.bottom()),
            (self.right(), self.bottom()),
            INK,
            0.75,
        );
        svg.text(
            x_center,
            self.bottom() + 36.0,
            x_label,
            11.0 * scale,
            "middle",
            INK,
        );
        let _ = write!(svg.body, r#"<g transform="rotate(-90 19 {y_center})">"#);
        svg.text(19.0, y_center, y_label, 11.0 * scale, "middle", INK);
        svg.body.push_str("</g>");
    }
}

fn padded(points: impl Iterator<Item = f64>, padding: f64) -> (f64, f64) {
    let (min, max) = points
        .filter(|v| v.is_finite())
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), v| {
            (min.min(v), max.max(v))
        });
    if min < max {
        let margin = (max - min) * padding;
        (min - margin, max + margin)
    } else if min.is_finite() {
        (min - 1.0, max + 1.0)
    } else {
        (-1.0, 1.0)
    }
}

fn ticks((min, max): (f64, f64), count: f64) -> impl Iterator<Item = f64> {
    let rough = (max - min) / count;
    let magnitude = 10f64.powf(rough.log10().floor());
    let step = [1.0, 2.0, 5.0, 10.0]
        .into_iter()
        .map(|v| v * magnitude)
        .find(|v| *v >= rough)
        .unwrap_or(rough);
    (0..=10)
        .map(move |i| ((min / step).ceil() + f64::from(i)) * step)
        .take_while(move |&x| x <= max)
}

fn tick_label(value: f64) -> String {
    if value.abs() >= 1e5 || (value != 0.0 && value.abs() < 0.001) {
        format!("{value:.1e}")
    } else {
        format!("{value:.3}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_owned()
    }
}

fn label(
    svg: &mut Svg,
    chart: &Chart,
    occupied: &mut Vec<[f64; 4]>,
    (x, y): (f64, f64),
    text: &str,
    color: &str,
    size: f64,
) {
    let width = (text.chars().count() as f64 * size * 0.59).min(chart.width);
    for distance in [11.0, 19.0, 29.0, 42.0, 60.0] {
        for (dx, dy) in [
            (0.0, -1.0),
            (1.0, -0.5),
            (-1.0, -0.5),
            (1.0, 0.8),
            (-1.0, 0.8),
            (0.0, 1.0),
        ] {
            let left = (x + dx * (distance + width / 2.0) - width / 2.0)
                .clamp(chart.left, chart.right() - width);
            let baseline = (y + dy * distance).clamp(chart.top + size, chart.bottom());
            let bounds = [
                left - 2.0,
                baseline - size,
                left + width + 2.0,
                baseline + 3.0,
            ];
            if occupied.iter().any(|b| {
                bounds[0] < b[2] && bounds[2] > b[0] && bounds[1] < b[3] && bounds[3] > b[1]
            }) {
                continue;
            }
            let anchor = x.clamp(left, left + width);
            svg.line((x, y), (anchor, baseline - size * 0.4), color, 0.45);
            svg.text(left, baseline, text, size, "start", color);
            occupied.push(bounds);
            return;
        }
    }
}

struct Svg {
    width: u32,
    body: String,
}

impl Svg {
    fn new(width: u32, height: u32) -> Self {
        Self {
            width,
            body: format!(
                r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}"><rect width="{width}" height="{height}" fill="white"/>"#
            ),
        }
    }

    fn line(&mut self, (x1, y1): (f64, f64), (x2, y2): (f64, f64), color: &str, width: f64) {
        let _ = write!(
            self.body,
            r#"<line x1="{x1:.3}" x2="{x2:.3}" y1="{y1:.3}" y2="{y2:.3}" stroke="{color}" stroke-width="{width}"/>"#
        );
    }

    fn point(&mut self, x: f64, y: f64, radius: f64, color: &str, opacity: f64) {
        let _ = write!(
            self.body,
            r#"<circle cx="{x:.3}" cy="{y:.3}" r="{radius}" fill="{color}" fill-opacity="{opacity}"/>"#
        );
    }

    fn rect(&mut self, (x, y): (f64, f64), (width, height): (f64, f64), color: &str) {
        let _ = write!(
            self.body,
            r#"<rect x="{x:.3}" y="{y:.3}" width="{width:.3}" height="{height:.3}" fill="{color}"/>"#
        );
    }

    fn text(&mut self, x: f64, y: f64, text: &str, size: f64, anchor: &str, color: &str) {
        let _ = write!(
            self.body,
            r#"<text x="{x:.3}" y="{y:.3}" font-family="sans-serif" font-size="{size}" text-anchor="{anchor}" fill="{color}">"#
        );
        for c in text.chars() {
            match c {
                '<' => self.body.push_str("&lt;"),
                '>' => self.body.push_str("&gt;"),
                '&' => self.body.push_str("&amp;"),
                _ if !c.is_control() => self.body.push(c),
                _ => {}
            }
        }
        self.body.push_str("</text>");
    }

    fn tree(mut self, fonts: &usvg::Options<'_>) -> Result<usvg::Tree> {
        self.body.push_str("</svg>");
        usvg::Tree::from_str(&self.body, fonts).context("Building plot")
    }

    fn save_pdf(self, path: &Path, fonts: &usvg::Options<'_>) -> Result<()> {
        save_pdf(&self.tree(fonts)?, path)
    }
}

fn save_pdf(tree: &usvg::Tree, path: &Path) -> Result<()> {
    let pdf = svg2pdf::to_pdf(
        tree,
        ConversionOptions::default(),
        PageOptions { dpi: 72.0 },
    )
    .map_err(|error| anyhow::anyhow!("Rendering {}: {error:?}", path.display()))?;
    fs::write(path, pdf)?;
    Ok(())
}

fn save_png(tree: &usvg::Tree, path: &Path) -> Result<()> {
    let mut pixmap = Pixmap::new(3900, 2700).context("Allocating volcano image")?;
    let scale = 300.0 / 72.0;
    resvg::render(
        tree,
        Transform::from_scale(scale, scale),
        &mut pixmap.as_mut(),
    );
    let mut encoder = png::Encoder::new(BufWriter::new(File::create(path)?), 3900, 2700);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    encoder.set_pixel_dims(Some(png::PixelDimensions {
        xppu: 11811,
        yppu: 11811,
        unit: png::Unit::Meter,
    }));
    let mut writer = encoder.write_header()?;
    writer.write_image_data(pixmap.data())?;
    writer.finish()?;
    Ok(())
}
