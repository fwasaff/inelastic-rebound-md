#!/usr/bin/env bash
# ==============================================================================
# scripts/analyze.sh
# Analyse an existing rebound_data.dat and generate all figures.
#
# Usage:
#   bash scripts/analyze.sh <data_dir> <v0> [angle] [T]
#
# Example:
#   bash scripts/analyze.sh data/outputs/v0_2.00_ang_0.0 2.0 0.0 0.1
# ==============================================================================

set -euo pipefail

DATA_DIR="${1:?Usage: analyze.sh <data_dir> <v0> [angle] [T]}"
V0="${2:?Provide v0}"
ANGLE="${3:-0.0}"
T="${4:-0.1}"

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

DATA_FILE="${DATA_DIR}/rebound_data.dat"

if [[ ! -f "$DATA_FILE" ]]; then
    echo "ERROR: Data file not found: $DATA_FILE" >&2
    exit 1
fi

echo "Analysing: $DATA_FILE"

python -m python.runner analyze \
    --data "$DATA_FILE" \
    --v0 "$V0" \
    --angle "$ANGLE" \
    --T "$T"

python -m python.runner plot trajectory \
    --data "$DATA_FILE" \
    --output-dir "${DATA_DIR}/figures"

echo "Figures saved to: ${DATA_DIR}/figures/"
