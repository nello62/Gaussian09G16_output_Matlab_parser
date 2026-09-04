#!/usr/bin/env bash
# Compiles G_gaussian_viewer.m into a standalone macOS app via MATLAB Compiler (mcc).
#
# Usage:
#   ./build_standalone.sh
#   MATLAB_ROOT=/Applications/MATLAB_R2024b.app ./build_standalone.sh
#
# Output goes to ./build (gitignored) as G_gaussian_viewer.app plus the
# mcc-generated run_G_gaussian_viewer.sh launcher and readme.txt.
#
# Requires a licensed MATLAB Compiler (license('test','Compiler') == 1).
# The build only runs on the platform it's built on (mcc does not
# cross-compile) -- run this script again on Windows/Linux to get a build
# for those platforms.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MATLAB_ROOT="${MATLAB_ROOT:-/Applications/MATLAB_R2023a.app}"
MCC="$MATLAB_ROOT/bin/mcc"
ENTRY="$SCRIPT_DIR/G_gaussian_viewer.m"
BUILD_DIR="$SCRIPT_DIR/build"

if [[ ! -x "$MCC" ]]; then
    echo "error: mcc not found at $MCC" >&2
    echo "Set MATLAB_ROOT to your MATLAB installation, e.g.:" >&2
    echo "  MATLAB_ROOT=/Applications/MATLAB_R2024b.app ./build_standalone.sh" >&2
    exit 1
fi

if [[ ! -f "$ENTRY" ]]; then
    echo "error: entry point not found: $ENTRY" >&2
    exit 1
fi

rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"

echo "Building with: $MCC"
echo "Entry point:   $ENTRY"
echo "Output dir:    $BUILD_DIR"
echo

"$MCC" -m "$ENTRY" \
    -I "$SCRIPT_DIR/../../G09" \
    -I "$SCRIPT_DIR/../../G16" \
    -I "$SCRIPT_DIR/../MOSurface" \
    -d "$BUILD_DIR" \
    -o G_gaussian_viewer

echo
# Both files always contain a fixed header row even with zero real entries
# (1 line for unresolvedSymbols.txt, 2 for mccExcludedFiles.log), so a
# plain non-empty check always fires -- compare against that baseline.
if [[ "$(wc -l < "$BUILD_DIR/unresolvedSymbols.txt")" -gt 1 ]]; then
    echo "warning: unresolved symbols detected -- see $BUILD_DIR/unresolvedSymbols.txt" >&2
fi
if [[ "$(wc -l < "$BUILD_DIR/mccExcludedFiles.log")" -gt 2 ]]; then
    echo "warning: files excluded from build -- see $BUILD_DIR/mccExcludedFiles.log" >&2
fi

echo "Build complete: $BUILD_DIR/G_gaussian_viewer.app"
echo "Run it with:    $BUILD_DIR/run_G_gaussian_viewer.sh $MATLAB_ROOT"
echo "(or any MATLAB Runtime R2023a installation in place of $MATLAB_ROOT)"
