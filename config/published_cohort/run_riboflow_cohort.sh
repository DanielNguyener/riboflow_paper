#!/bin/bash

if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
else
    for c in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3"; do
        [ -f "$c/etc/profile.d/conda.sh" ] && source "$c/etc/profile.d/conda.sh" && break
    done
fi
conda activate ribo_genome

set -uo pipefail

REPO="${RIBOFLOW_REPO:?set RIBOFLOW_REPO to the RiboFlow_v2 checkout}"
YAML_DIR="$REPO/RiboFlow_YAMLs_cohort" 
PROFILE=${NXF_PROFILE-lonestar6}
FORCE=${FORCE:-0}

DONE_GLOB='ribo/*.ribo'

cd "$REPO"

sample_done() {
    local sample=$1
    local base="$REPO/output/$LABEL/$sample"
    [ -d "$base" ] && find "$base"/$DONE_GLOB -size +0c 2>/dev/null | grep -q .
}

fastqs_ready() {
    local yaml=$1 fq fqs
    fqs=$(grep -oE '\$\{FASTQ_DIR\}/[^[:space:]]+\.fastq\.gz' "$yaml" \
          | sed "s#\\\${FASTQ_DIR}#${FASTQ_DIR:?set FASTQ_DIR to the input FASTQ root}#" \
          | sort -u)
    for fq in $fqs; do
        [ -s "$fq" ] || return 1
        [ "$(od -An -tx1 -N2 "$fq" 2>/dev/null | tr -d ' \n')" = "1f8b" ] || return 1
    done
    return 0
}

YAMLS=("$YAML_DIR"/*.yaml)
TOTAL=${#YAMLS[@]}
COUNT=0
DONE=(); SKIPPED=(); FAILED=()

PROFILE_ARG=(); [ -n "$PROFILE" ] && PROFILE_ARG=(-profile "$PROFILE")

for yaml in "${YAMLS[@]}"; do
    COUNT=$((COUNT + 1))
    SAMPLE=$(basename "$yaml" .yaml)
    echo ""
    echo "  [$COUNT/$TOTAL]  $SAMPLE   $(date)"

    if [ "$FORCE" != "1" ] && sample_done "$SAMPLE"; then
        DONE+=("$SAMPLE"); continue
    fi

    if ! fastqs_ready "$yaml"; then
        SKIPPED+=("$SAMPLE"); continue
    fi

    nextflow run main.nf "${PROFILE_ARG[@]}" -params-file "$yaml" -resume \
        || { echo "ERROR: $SAMPLE failed (exit $?), queued for retry"; FAILED+=("$yaml"); }
done

if [[ ${#FAILED[@]} -gt 0 ]]; then
    echo ""
    echo "  Retrying ${#FAILED[@]} failed sample(s) — $(date)"
    for yaml in "${FAILED[@]}"; do
        SAMPLE=$(basename "$yaml" .yaml)
        echo "  [retry] $SAMPLE  $(date)"
        nextflow run main.nf "${PROFILE_ARG[@]}" -params-file "$yaml" -resume \
            || echo "ERROR: $SAMPLE failed on retry (exit $?)"
    done
fi

echo ""
echo "Pass complete — $(date)"
echo "  finished/skipped-done : ${#DONE[@]}  ${DONE[*]:-}"
echo "  skipped (no fastq)    : ${#SKIPPED[@]}  ${SKIPPED[*]:-}"
echo "  failed (after retry)  : ${#FAILED[@]}  $(printf '%s ' "${FAILED[@]##*/}")"
