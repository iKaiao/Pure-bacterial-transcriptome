#!/usr/bin/env bash
set -euo pipefail

export CARGO_HOME="$SRC_DIR/.cargo-build"
export CARGO_TARGET_DIR="$SRC_DIR/.cargo-target"
target=${CARGO_BUILD_TARGET:-$(rustc -vV | sed -n 's/^host: //p')}
linker_var="CARGO_TARGET_$(printf '%s' "$target" | tr '[:lower:]-' '[:upper:]_')_LINKER"
export "$linker_var=$CC"

cargo build --locked --release --target "$target"
mkdir -p "$PREFIX/bin"
install -m755 "$CARGO_TARGET_DIR/$target/release/pure-bacterial-transcriptome" \
  "$PREFIX/bin/pure-bacterial-transcriptome"
