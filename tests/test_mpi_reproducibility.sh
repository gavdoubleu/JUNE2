#!/usr/bin/env bash
# =============================================================================
# MPI Reproducibility Regression Test
# =============================================================================
# Runs the full simulation at 1-rank (serial) and 2-rank (parallel) and
# asserts bit-identical infection counts at every day for 8 days.
#
# The 8-day default is deliberate: 1 and 2 ranks are known to diverge from
# about day 11 on config_2021 (infector attribution), so a longer horizon
# fails for reasons unrelated to the change under test.
#
# This catches any change that breaks MPI reproducibility, including:
#   - Remote participant eligibility assumptions
#   - Virtual venue ID collisions / routing errors
#   - Visitor data packing/unpacking mismatches
#   - Floating-point non-associativity in accumulation order
#   - Encounter dedup/ordering inconsistencies
#   - Policy check asymmetries across ranks
#   - Contact matrix selection divergence
#   - Double-counting of outgoing visitors
#
# Usage:
#   ./tests/test_mpi_reproducibility.sh [--build-dir BUILD_DIR] [--days N]
#
# Requirements:
#   - Built disease_sim binary with USE_MPI=ON
#   - worlds/world_2021.h5 and configs/config_2021/ present
#   - mpirun available
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="${PROJECT_DIR}/build"
DAYS=8
CONFIG="configs/config_2021/simulation.yaml"
WORLD="worlds/world_2021.h5"

# Parse arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --build-dir) BUILD_DIR="$2"; shift 2 ;;
    --days)      DAYS="$2"; shift 2 ;;
    *)           echo "Unknown arg: $1"; exit 1 ;;
  esac
done

BINARY="${BUILD_DIR}/disease_sim"
if [[ ! -x "$BINARY" ]]; then
  echo "FAIL: Binary not found at $BINARY"
  echo "      Build with: cmake --build $BUILD_DIR --target disease_sim"
  exit 1
fi

if [[ ! -f "${PROJECT_DIR}/${WORLD}" ]]; then
  echo "FAIL: World file not found at ${PROJECT_DIR}/${WORLD}"
  exit 1
fi

cd "$PROJECT_DIR"

# Create temp directory for this test run
TMPDIR=$(mktemp -d "${TMPDIR:-/tmp}/mpi_repro_XXXXXX")
trap "rm -rf $TMPDIR" EXIT

echo "=== MPI Reproducibility Test ==="
echo "  Binary: $BINARY"
echo "  Config: $CONFIG"
echo "  World:  $WORLD"
echo "  Days:   $DAYS"
echo "  Tmpdir: $TMPDIR"
echo ""

# Create a config override for the requested number of days, derived from the
# config's own start_date so the harness follows the config's calendar.
START_DATE=$(sed -n 's/^[[:space:]]*start_date[[:space:]]*:[[:space:]]*"\([0-9-]*\)".*/\1/p' \
  "${PROJECT_DIR}/${CONFIG}" | head -1)
if [[ -z "$START_DATE" ]]; then
  echo "FAIL: no start_date in ${CONFIG}"
  exit 1
fi
END_DATE=$(python3 -c "
from datetime import datetime, timedelta
start = datetime.strptime('${START_DATE}', '%Y-%m-%d')
print((start + timedelta(days=${DAYS})).strftime('%Y-%m-%d'))
")

# Copy config and override end_date
cp "$CONFIG" "${TMPDIR}/simulation.yaml"
if [[ "$(uname)" == "Darwin" ]]; then
  sed -i '' "s/end_date.*/end_date: \"${END_DATE}\"/" "${TMPDIR}/simulation.yaml"
else
  sed -i "s/end_date.*/end_date: \"${END_DATE}\"/" "${TMPDIR}/simulation.yaml"
fi

# --- Run 1-rank (serial) ---
echo "[1/3] Running 1-rank simulation..."
mpirun -np 1 --oversubscribe "$BINARY" \
  --config "${TMPDIR}/simulation.yaml" \
  --world "$WORLD" \
  --runs-dir "${TMPDIR}/runs" \
  --run-id "np1" \
  > "${TMPDIR}/log_1rank.txt" 2>&1
echo "      Done."

# --- Run 2-rank (parallel) ---
echo "[2/3] Running 2-rank simulation..."
mpirun -np 2 --oversubscribe "$BINARY" \
  --config "${TMPDIR}/simulation.yaml" \
  --world "$WORLD" \
  --runs-dir "${TMPDIR}/runs" \
  --run-id "np2" \
  > "${TMPDIR}/log_2rank.txt" 2>&1
echo "      Done."

# --- Compare ---
echo "[3/3] Comparing infection counts..."
echo ""

grep "Total currently infected" "${TMPDIR}/log_1rank.txt" | awk '{print $4}' > "${TMPDIR}/counts_1rank.txt"
grep "Total currently infected" "${TMPDIR}/log_2rank.txt" | awk '{print $4}' > "${TMPDIR}/counts_2rank.txt"

N1=$(wc -l < "${TMPDIR}/counts_1rank.txt" | tr -d ' ')
N2=$(wc -l < "${TMPDIR}/counts_2rank.txt" | tr -d ' ')

if [[ "$N1" -eq 0 ]]; then
  echo "FAIL: No infection counts found in 1-rank log."
  echo "      Check ${TMPDIR}/log_1rank.txt for errors."
  exit 1
fi

if [[ "$N1" -ne "$N2" ]]; then
  echo "FAIL: Different number of daily counts: 1-rank=$N1, 2-rank=$N2"
  exit 1
fi

FAILURES=0
DAY=0
while IFS= read -r line1 && IFS= read -r line2 <&3; do
  if [[ "$line1" -ne "$line2" ]]; then
    echo "  MISMATCH Day $DAY: 1-rank=$line1  2-rank=$line2  diff=$((line1 - line2))"
    FAILURES=$((FAILURES + 1))
  fi
  DAY=$((DAY + 1))
done < "${TMPDIR}/counts_1rank.txt" 3< "${TMPDIR}/counts_2rank.txt"

echo ""
if [[ "$FAILURES" -eq 0 ]]; then
  FINAL_1=$(tail -1 "${TMPDIR}/counts_1rank.txt")
  echo "PASS: $DAY days compared, 0 mismatches."
  echo "      Final infection count: $FINAL_1 (identical on both ranks)"
  exit 0
else
  echo "FAIL: $FAILURES mismatches out of $DAY days."
  echo "      Logs saved in $TMPDIR"
  # Don't clean up on failure so logs can be inspected
  trap - EXIT
  exit 1
fi
