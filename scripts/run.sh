#!/usr/bin/env bash
# ==============================================================================
# scripts/run.sh
# Full pipeline: equilibrate -> single impact -> analyze -> plot
#
# Usage:
#   bash scripts/run.sh [v0] [angle] [T]
#
# Defaults:
#   v0    = 2.0  (incident speed, LJ units)
#   angle = 0.0  (degrees from normal)
#   T     = 0.1  (surface temperature, eps/kB)
# ==============================================================================

set -euo pipefail

V0="${1:-2.0}"
ANGLE="${2:-0.0}"
T="${3:-0.1}"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

echo "========================================"
echo "  Inelastic Rebound MD – Full Pipeline"
echo "========================================"
echo "  v0    = $V0 sigma/tau"
echo "  angle = $ANGLE degrees"
echo "  T*    = $T eps/kB"
echo "========================================"

# --- Step 1: Surface equilibration -------------------------------------------
echo ""
echo "[1/3] Equilibrating surface..."
python -m python.runner equilibrate \
    --T "$T" \
    --Nx 8 --Ny 8 --Nz 6

# --- Step 2: Impact simulation ------------------------------------------------
echo ""
echo "[2/3] Running impact simulation..."
python -m python.runner impact \
    --v0 "$V0" \
    --angle "$ANGLE" \
    --T "$T"

# --- Step 3: Analysis and plotting -------------------------------------------
DATA_DIR="data/outputs/v0_${V0}_ang_${ANGLE}"

echo ""
echo "[3/3] Generating figures..."
python -m python.runner plot trajectory \
    --data "${DATA_DIR}/rebound_data.dat" \
    --output-dir "data/outputs/figures"

echo ""
echo "========================================"
echo "  Pipeline complete."
echo "  Results: ${DATA_DIR}/"
echo "  Figures: data/outputs/figures/"
echo "========================================"
