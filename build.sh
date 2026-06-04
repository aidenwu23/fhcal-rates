#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$ROOT_DIR/build"
MODULE="background-rates"

usage() {
  cat <<EOF2
Usage: $0 [clean] [fast] [cmake-opts...]

Options:
  clean   Remove the build directory before configuring
  fast    Skip configure, just build/install

Env:
  EIC_SHELL_PREFIX        Hint for dependencies (EIC env)
  BUILD_NPROC             Parallel build jobs (default: nproc or 1)
  LFHCAL_BUILD_TYPE       CMAKE_BUILD_TYPE (default: RelWithDebInfo)
  LFHCAL_INSTALL_PREFIX   Install prefix (default: build/install)
EOF2
}

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

clean=0
fast=0
extra_opts=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    clean) clean=1 ;;
    fast) fast=1 ;;
    *) extra_opts+=("$1") ;;
  esac
  shift
done

install_prefix="${LFHCAL_INSTALL_PREFIX:-${BUILD_DIR}}"
nproc="${BUILD_NPROC:-}"
if [[ -z "$nproc" ]]; then
  if command -v nproc >/dev/null 2>&1; then
    nproc="$(nproc)"
  else
    nproc=1
  fi
fi

existing_prefix_path="${CMAKE_PREFIX_PATH:-}"
existing_prefix_path="${existing_prefix_path//:/;}"
cmake_prefix_entries=()
[[ -n "${EIC_SHELL_PREFIX:-}" ]] && cmake_prefix_entries+=("$EIC_SHELL_PREFIX")
cmake_prefix_entries+=("$install_prefix")
if [[ -z "${EIC_SHELL_PREFIX:-}" ]]; then
  cmake_prefix_entries+=("/opt/local")
fi
[[ -n "$existing_prefix_path" ]] && cmake_prefix_entries+=("$existing_prefix_path")
export CMAKE_PREFIX_PATH="$(IFS=';'; echo "${cmake_prefix_entries[*]}")"

build_type="${LFHCAL_BUILD_TYPE:-RelWithDebInfo}"

cmake_gen=(
  cmake
  -S "$ROOT_DIR"
  -B "$BUILD_DIR"
  -DCMAKE_BUILD_TYPE="$build_type"
  -DCMAKE_INSTALL_PREFIX="$install_prefix"
  -DCMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH"
  -DCMAKE_FIND_DEBUG_MODE=OFF
)
cmake_build=(cmake --build "$BUILD_DIR" -j"$nproc")
cmake_install=(cmake --install "$BUILD_DIR")

cat <<EOF2

BUILDING:
module    = $MODULE
clean     = $clean
fast      = $fast
prefix    = $install_prefix
nproc     = $nproc
extraOpts = ${extra_opts[*]:-(none)}
EOF2

if [[ $clean -eq 1 ]]; then
  rm -rf "$BUILD_DIR"
fi

[[ $fast -eq 0 ]] && "${cmake_gen[@]}" "${extra_opts[@]}"
"${cmake_build[@]}"
"${cmake_install[@]}"
