#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREFIX="${PREFIX:-$ROOT_DIR/.local}"
BIN_DIR="$PREFIX/bin"
SRC_DIR="$PREFIX/src"
BUILD_DIR="$PREFIX/build"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)
      PREFIX="$2"
      BIN_DIR="$PREFIX/bin"
      SRC_DIR="$PREFIX/src"
      BUILD_DIR="$PREFIX/build"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

mkdir -p "$BIN_DIR" "$SRC_DIR" "$BUILD_DIR"

for cmd in curl tar make; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "[ERROR] Required command not found: $cmd" >&2
    exit 1
  fi
done

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "[WARN] No active conda environment detected."
  echo "[WARN] Please activate your pipeline environment first (e.g., conda activate cbdating)."
fi

echo "[INFO] Installing Python dependencies (wLogDate + MD-Cat requirements)..."
python3 -m pip install --upgrade pip
python3 -m pip install "wlogdate==1.0.2" "dendropy>=5.0.0"

if [[ -n "${CONDA_PREFIX:-}" ]] && command -v conda >/dev/null 2>&1; then
  echo "[INFO] Installing TreePL build dependencies into active conda env..."
  conda install -y -c conda-forge adol-c nlopt
fi

fetch_and_extract() {
  local url="$1"
  local out_dir="$2"
  local tarball="$BUILD_DIR/$(basename "$out_dir").tar.gz"

  rm -rf "$out_dir"
  mkdir -p "$out_dir"
  curl -fsSL "$url" -o "$tarball"
  tar -xzf "$tarball" -C "$out_dir" --strip-components=1
}

install_mdcat_wrapper() {
  local mdcat_src="$SRC_DIR/MD-Cat"
  echo "[INFO] Downloading MD-Cat source..."
  fetch_and_extract "https://github.com/uym2/MD-Cat/archive/refs/heads/master.tar.gz" "$mdcat_src"

  cat > "$BIN_DIR/md_cat.py" <<EOF
#!/usr/bin/env bash
exec python3 "$mdcat_src/md_cat.py" "\$@"
EOF
  chmod +x "$BIN_DIR/md_cat.py"
}

install_wlogdate_wrapper() {
  local launch_path
  launch_path="$(command -v launch_wLogDate.py || true)"

  if [[ -z "$launch_path" ]]; then
    echo "[ERROR] launch_wLogDate.py was not installed by pip." >&2
    exit 1
  fi

  ln -sf "$launch_path" "$BIN_DIR/launch_wLogDate.py"
}

install_lsd2() {
  local lsd2_src="$SRC_DIR/lsd2"
  echo "[INFO] Downloading and building LSD2..."
  fetch_and_extract "https://github.com/tothuhien/lsd2/archive/refs/heads/master.tar.gz" "$lsd2_src"
  make -C "$lsd2_src/src"
  install -m 0755 "$lsd2_src/src/lsd2" "$BIN_DIR/lsd2"
}

install_treepl() {
  local treepl_src="$SRC_DIR/treePL"
  echo "[INFO] Downloading and building TreePL..."
  fetch_and_extract "https://github.com/blackrim/treePL/archive/refs/heads/master.tar.gz" "$treepl_src"

  pushd "$treepl_src/src" >/dev/null
  chmod +x configure

  local cppflags=""
  local ldflags=""
  if [[ -n "${CONDA_PREFIX:-}" ]]; then
    cppflags="-I$CONDA_PREFIX/include"
    ldflags="-L$CONDA_PREFIX/lib"
  fi

  ./configure CPPFLAGS="$cppflags" LDFLAGS="$ldflags"
  make
  popd >/dev/null

  install -m 0755 "$treepl_src/src/treePL" "$BIN_DIR/treePL"
}

install_mdcat_wrapper
install_wlogdate_wrapper
install_lsd2
install_treepl

echo

echo "[INFO] Bootstrap completed."
echo "[INFO] Add these tools to your PATH:"
echo "       export PATH=\"$BIN_DIR:\$PATH\""
