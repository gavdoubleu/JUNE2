#!/usr/bin/env bash
# =============================================================================
# MPI determinism: full disease_sim end-to-end invariants
# =============================================================================
# Diff the HDF5 datasets that exercise coordinated encounters, follows, infections,
# and symptom transitions. The config_2021 fixture does not enable OOE
# relationship generation, so /events/relationships is not a required input
# to this gate.
#
# It also compares seeded counts and the existing daily infection-count logs.
#
# Strategy: run disease_sim at np=1/2/3 into three HDF5 outputs, then use
# a small python+h5py helper to dump each dataset to a sorted canonical
# form (row order varies across runs due to per-rank event buffering;
# determinism is about the *set* of records, not file-offset order).
# A failure prints the first diverging records with full field context.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-${PROJECT_DIR}/build}"
DAYS="${DAYS:-10}"
SEED="${SEED:-12345}"
CONFIG="${CONFIG:-configs/config_2021/simulation.yaml}"
WORLD="${WORLD:-worlds/world_2021.h5}"
# Optional infection-seeds file, overriding the one the config names. Set it to
# exercise a seeding variant (e.g. structured seeds above the partition level)
# without a second config directory.
SEEDS="${SEEDS:-}"
# Space-separated list of rank counts to compare. Must include 1 as the
# reference.
NPS="${NPS:-1 2 3}"
# Further datasets to compare, written as DATASETS entries and separated by
# spaces, for a config that switches on more output (vaccination, say). Each
# one is required to be present and non-empty, so leave it empty for a config
# that does not write them.
EXTRA_DATASETS="${EXTRA_DATASETS:-}"

BINARY="${BUILD_DIR}/disease_sim"
if [[ ! -x "$BINARY" ]]; then
  echo "FAIL: binary not found at $BINARY"
  exit 1
fi
# The world file is not committed, so CI has no way to run this. Skip (ctest
# SKIP_RETURN_CODE 77) rather than fail, matching test_checkpoint_determinism.
if [[ ! -f "${PROJECT_DIR}/${WORLD}" ]]; then
  echo "SKIP: world not found at ${PROJECT_DIR}/${WORLD}"
  exit 77
fi

cd "$PROJECT_DIR"

TMP=$(mktemp -d "${TMPDIR:-/tmp}/mpi_full_XXXXXX")
CLEAN_ON_EXIT=1
trap '[[ $CLEAN_ON_EXIT -eq 1 ]] && rm -rf "$TMP"' EXIT

echo "=== MPI full-reproducibility test ==="
echo "  binary: $BINARY"
echo "  days:   $DAYS"
echo "  nps:    $NPS"
echo "  seed:   $SEED"
[[ -n "$SEEDS" ]] && echo "  seeds:  $SEEDS"
echo "  tmp:    $TMP"
echo ""

# --- config override for day count -------------------------------------------
# Derived from the config's own start_date, so the harness follows the config's
# calendar rather than pinning one of its own.
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
cp "$CONFIG" "$TMP/simulation.yaml"
if [[ "$(uname)" == "Darwin" ]]; then
  sed -i '' "s/end_date.*/end_date: \"${END_DATE}\"/" "$TMP/simulation.yaml"
else
  sed -i "s/end_date.*/end_date: \"${END_DATE}\"/" "$TMP/simulation.yaml"
fi

# --- run disease_sim at each np ----------------------------------------------
# Each run writes into $TMP/runs/np${NP}/. We then move the merged HDF5
# to $TMP/sim_np${NP}.h5 for the canonical-diff step below.
for NP in $NPS; do
  echo "[run] np=$NP"
  RUN_ID="np${NP}"
  RUN_DIR="$TMP/runs/${RUN_ID}"
  OUT="$TMP/sim_np${NP}.h5"
  SEEDS_ARG=()
  [[ -n "$SEEDS" ]] && SEEDS_ARG=(--infection_seeds "$SEEDS")
  mpirun -np "$NP" --oversubscribe "$BINARY" \
    --config "$TMP/simulation.yaml" \
    --world "$WORLD" \
    --seed "$SEED" \
    ${SEEDS_ARG[@]+"${SEEDS_ARG[@]}"} \
    --runs-dir "$TMP/runs" \
    --run-id "$RUN_ID" \
    > "$TMP/log_np${NP}.txt" 2>&1 || {
      echo "FAIL: disease_sim exited non-zero at np=$NP"
      echo "      last 20 log lines:"
      tail -20 "$TMP/log_np${NP}.txt"
      CLEAN_ON_EXIT=0
      exit 1
    }
  if [[ ! -f "$RUN_DIR/simulation_events.h5" ]]; then
    echo "FAIL: no HDF5 output produced at np=$NP (expected $RUN_DIR/simulation_events.h5)"
    CLEAN_ON_EXIT=0
    exit 1
  fi
  mv "$RUN_DIR/simulation_events.h5" "$OUT"
done

# --- canonicalize every dataset of interest at every np ----------------------
# The Python helper is shared with checkpoint_determinism_check.py.
HDF5_HELPER="$SCRIPT_DIR/hdf5_determinism.py"
# Each entry: "kind:path[|fields=a,b,c]". Every field is compared unless a
# fields= list narrows it.
DATASETS=(
  "dataset:/events/coordinated_encounters"
  "dataset:/events/follows"
  "dataset:/events/infections"
  "dataset:/events/symptom_changes"
)
# shellcheck disable=SC2206  # split on spaces on purpose
DATASETS+=($EXTRA_DATASETS)
# config_2021 does not configure the optional OOE relationship producer, so
# /events/relationships is absent by design. It is reported as a fixture
# limitation rather than treated as a required, non-exercised dataset.
# NOTE: /lookups/profile_assignments is NOT diffed here — the multi-rank
# event merger (main.cpp:476) discards per-rank lookup tables, so the
# merged file only contains rank 0's slice at np>1. Profile-assignment
# determinism is therefore not directly asserted by this HDF5 comparison.
# Profile drift that affects attendance, encounters, follows, infections, or
# symptoms is still observable in the canonical event streams above.

FAIL=0
for SPEC in "${DATASETS[@]}"; do
  MAIN="${SPEC%%|*}"
  EXTRA=""
  if [[ "$SPEC" == *"|"* ]]; then EXTRA="${SPEC#*|}"; fi
  KIND="${MAIN%%:*}"
  PATH_="${MAIN#*:}"
  SAFE=$(echo "$PATH_" | tr '/' '_')
  for NP in $NPS; do
    CANON="$TMP/canon${SAFE}_np${NP}.txt"
    python3 "$HDF5_HELPER" dump "$TMP/sim_np${NP}.h5" "$KIND" "$PATH_" \
      "$CANON" $EXTRA
    if grep -q '^# MISSING ' "$CANON"; then
      echo "FAIL: expected dataset ${PATH_} is missing at np=${NP}"
      CLEAN_ON_EXIT=0
      FAIL=1
    elif grep -q '^# EMPTY ' "$CANON"; then
      # Every dataset in DATASETS is required to exercise this fixture's
      # transmission, encounter, seeding, or cross-domain paths.
      echo "FAIL: expected non-empty dataset ${PATH_} is empty at np=${NP}"
      CLEAN_ON_EXIT=0
      FAIL=1
    fi
  done
done

# --- diff against np=1 reference ---------------------------------------------
REF_NP=$(echo "$NPS" | awk '{print $1}')
if [[ "$REF_NP" != "1" ]]; then
  echo "WARN: first np in NPS is not 1 — using np=$REF_NP as reference"
fi

for SPEC in "${DATASETS[@]}"; do
  MAIN="${SPEC%%|*}"
  EXTRA=""
  if [[ "$SPEC" == *"|"* ]]; then EXTRA="${SPEC#*|}"; fi
  PATH_="${MAIN#*:}"
  SAFE=$(echo "$PATH_" | tr '/' '_')
  REF="$TMP/canon${SAFE}_np${REF_NP}.txt"
  REF_LINES=$(wc -l < "$REF" | tr -d ' ')
  SPEC_FAIL=0
  for NP in $NPS; do
    [[ "$NP" == "$REF_NP" ]] && continue
    CUR="$TMP/canon${SAFE}_np${NP}.txt"
    DIFF="$TMP/diff${SAFE}_np${NP}.txt"
    if ! python3 "$HDF5_HELPER" compare "$TMP/sim_np${REF_NP}.h5" \
         "$TMP/sim_np${NP}.h5" "$PATH_" $EXTRA > "$DIFF"; then
      echo ""
      echo "FAIL: ${PATH_} diverges between np=${REF_NP} and np=${NP}"
      echo "      ref ($REF_LINES lines): $REF"
      echo "      cur: $CUR"
      sed -n '1,30p' "$DIFF"
      CLEAN_ON_EXIT=0
      FAIL=1
      SPEC_FAIL=1
    fi
  done
  if [[ $SPEC_FAIL -eq 0 ]]; then
    echo "  PASS ${PATH_}: ${REF_LINES} records identical across $NPS"
  fi
done

# --- also: infection count line-by-line (cheaper sanity check) ---------------
echo ""
echo "=== infection-count log diff ==="
for NP in $NPS; do
  grep "Total currently infected" "$TMP/log_np${NP}.txt" \
    | awk '{print $4}' > "$TMP/infect_np${NP}.txt"
  # A structured seed places an absolute count, so the per-step seeded totals
  # are an invariant in their own right: they diverge before the day counts do
  # when a seed resolves per rank rather than globally.
  grep "\[INFECTION SEED\] Seeded" "$TMP/log_np${NP}.txt" \
    | awk '{print $4}' > "$TMP/seeded_np${NP}.txt"
done
# Two empty files compare equal, so a run that seeds nobody would pass the
# check below without comparing anything. Every config this harness runs
# seeds someone; if the reference run did not, that is the failure.
if [[ ! -s "$TMP/seeded_np${REF_NP}.txt" ]]; then
  echo "FAIL: np=${REF_NP} logged no '[INFECTION SEED] Seeded' lines, so the"
  echo "      seeded-count comparison would pass without comparing anything"
  CLEAN_ON_EXIT=0
  FAIL=1
fi
for NP in $NPS; do
  [[ "$NP" == "$REF_NP" ]] && continue
  if ! diff -q "$TMP/seeded_np${REF_NP}.txt" "$TMP/seeded_np${NP}.txt" > /dev/null; then
    echo "FAIL: seeded counts diverge between np=${REF_NP} and np=${NP}"
    diff -u "$TMP/seeded_np${REF_NP}.txt" "$TMP/seeded_np${NP}.txt" | head -20
    CLEAN_ON_EXIT=0
    FAIL=1
  else
    echo "  PASS seeded counts: identical np=${REF_NP} vs np=${NP}"
  fi
  if ! diff -q "$TMP/infect_np${REF_NP}.txt" "$TMP/infect_np${NP}.txt" > /dev/null; then
    echo "FAIL: infection counts diverge between np=${REF_NP} and np=${NP}"
    diff -u "$TMP/infect_np${REF_NP}.txt" "$TMP/infect_np${NP}.txt" | head -20
    CLEAN_ON_EXIT=0
    FAIL=1
  else
    DAYS_CMP=$(wc -l < "$TMP/infect_np${REF_NP}.txt" | tr -d ' ')
    echo "  PASS infection counts: ${DAYS_CMP} days identical np=${REF_NP} vs np=${NP}"
  fi
done

echo ""
if [[ $FAIL -eq 0 ]]; then
  echo "ALL PASS — coordinated encounters, follows, infections, symptom changes, and infection counts are identical across $NPS"
  exit 0
else
  echo "FAIL — reproducibility invariant violated. Temp files kept at $TMP"
  exit 1
fi
