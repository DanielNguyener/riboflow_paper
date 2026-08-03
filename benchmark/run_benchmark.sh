#!/usr/bin/env bash
set -uo pipefail

BENCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "$BENCH_DIR/.." && pwd)"   # --pipeline overrides

PROFILE="conda,default"
REPLICATES=""          # empty => use the per-scenario value from scenarios.tsv
SCENARIO_FILTER=""
DRY_RUN=0
OUTDIR=""
CPUS_OVERRIDE=""       # --cpus N       : pin executor to N cores (default: nproc)
MEM_OVERRIDE=""        # --memory 64.GB : pin executor memory (default: 3/4 of host RAM)
TASK_CPUS_OVERRIDE=""  # --task-cpus N   : default per-process cpu request
TASK_MEM_OVERRIDE=""   # --task-memory N : default per-process memory request, GB
STAR_CPUS_OVERRIDE=""  # --star-cpus N   : cpus for STAR_ALIGN/_RNASEQ/_INDEX
STAR_MEM_OVERRIDE=""   # --star-memory N : memory for STAR, GB
EXTRA_CONFIG=""        # --extra-config <file> : layered after benchmark.config
REFERENCES="$BENCH_DIR/references/human_server.yaml"

BIND_MODE="auto"       # auto | on | off

while [[ $# -gt 0 ]]; do
    case "$1" in
        --profile)     PROFILE="$2"; shift 2 ;;
        --pipeline)    PIPELINE_ROOT="$(cd "$2" && pwd)"; shift 2 ;;
        --cpus)        CPUS_OVERRIDE="$2"; shift 2 ;;
        --memory)      MEM_OVERRIDE="$2"; shift 2 ;;
        --task-cpus)   TASK_CPUS_OVERRIDE="$2"; shift 2 ;;
        --task-memory) TASK_MEM_OVERRIDE="$2"; shift 2 ;;
        --star-cpus)   STAR_CPUS_OVERRIDE="$2"; shift 2 ;;
        --star-memory) STAR_MEM_OVERRIDE="$2"; shift 2 ;;
        --extra-config) EXTRA_CONFIG="$2"; shift 2 ;;
        --references)  REFERENCES="$2"; shift 2 ;;
        --bind)        BIND_MODE="on"; shift ;;
        --no-bind)     BIND_MODE="off"; shift ;;
        --replicates)  REPLICATES="$2"; shift 2 ;;
        --scenarios)   SCENARIO_FILTER="$2"; shift 2 ;;
        --outdir)      OUTDIR="$2"; shift 2 ;;
        --dry-run)     DRY_RUN=1; shift ;;
        -h|--help)     sed -n '2,17p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

cd "$PIPELINE_ROOT" || exit 1

exec < /dev/null

die() { echo "ERROR: $*" >&2; exit 1; }
log() { printf '\n\033[1m==> %s\033[0m\n' "$*"; }

write_inputs_json() {
    python3 - "$1" "$2" <<'PY'
import json, os, sys, yaml

params = yaml.safe_load(open(sys.argv[1])) or {}
entries, missing = [], []

for kind, block in (("ribo", params.get("input") or {}),
                    ("rnaseq", params.get("rnaseq") or {})):
    base = block.get("fastq_base") or ""
    for sample, lanes in (block.get("fastq") or {}).items():
        for lane in (lanes if isinstance(lanes, list) else [lanes]):
            for path in (lane if isinstance(lane, list) else [lane]):
                if not isinstance(path, str):
                    continue
                full = os.path.join(base, path)
                if not os.path.exists(full):
                    missing.append(full)
                    continue
                entries.append({"kind": kind, "sample": sample,
                                "path": full, "bytes": os.path.getsize(full)})

if missing:
    print("declared FASTQ files do not exist:", file=sys.stderr)
    for m in missing:
        print(f"    {m}", file=sys.stderr)
    sys.exit(1)

json.dump({"n_files": len(entries),
           "n_ribo_files": sum(e["kind"] == "ribo" for e in entries),
           "n_rnaseq_files": sum(e["kind"] == "rnaseq" for e in entries),
           "n_samples": len({e["sample"] for e in entries}),
           "input_bytes": sum(e["bytes"] for e in entries),
           "files": entries},
          open(sys.argv[2], "w"), indent=2)
PY
}

nextflow_works() { nextflow -version >/dev/null 2>&1; }

if ! nextflow_works; then
    set +u
    for hook in "$HOME/miniconda3" /opt/miniconda3 "$HOME/mambaforge"; do
        [[ -f "$hook/etc/profile.d/conda.sh" ]] || continue
        source "$hook/etc/profile.d/conda.sh" 2>/dev/null || continue
        conda activate ribo_genome 2>/dev/null || continue
        nextflow_works && break
    done
    set -u
fi
export NXF_OPTS="${NXF_OPTS:--Xms1g -Xmx4g}"

log "Preflight"
[[ -f "$PIPELINE_ROOT/main.nf" ]] || die "no main.nf under $PIPELINE_ROOT (pass --pipeline <dir>)"
command -v nextflow >/dev/null || die "nextflow not on PATH (activate the ribo_genome env)"
nextflow -version >/dev/null 2>&1 || die "nextflow present but cannot start - check JAVA_HOME (needs JDK 17+)"
command -v python3  >/dev/null || die "python3 not on PATH"
python3 -c 'import yaml' 2>/dev/null || die "PyYAML missing (pip install pyyaml)"

[[ -r /proc/meminfo ]] || die "no /proc: peak_rss and %cpu would come back empty.
       This harness requires Linux."

if [[ "$PROFILE" != *conda* && "$PROFILE" != *apptainer* && "$PROFILE" != *docker* ]]; then
    for tool in STAR bowtie2 samtools cutadapt bedtools umi_tools umicollapse \
                bamCoverage ribopy rfc; do
        command -v "$tool" >/dev/null || die "'$tool' not on PATH. Either activate the ribo_genome env
       (conda env create -f environment.yaml && conda activate ribo_genome)
       or use --profile conda,default to let Nextflow build the env itself."
    done
fi

AVAIL_GB=$(df -Pk . | awk 'NR==2{printf "%d", $4/1048576}')
[[ "${AVAIL_GB:-0}" -ge 20 ]] || die "only ${AVAIL_GB}GB free; need >= 20GB"
echo "  pipeline  $PIPELINE_ROOT"
echo "  nextflow  $(nextflow -version 2>&1 | grep -oE 'version [0-9.]+' | head -1)"
echo "  free disk ${AVAIL_GB}GB"

HOST_CPUS=$(nproc)
HOST_MEM_B=$(( $(awk '/MemTotal/{print $2}' /proc/meminfo) * 1024 ))
export BENCH_EXECUTOR_CPUS="${CPUS_OVERRIDE:-$HOST_CPUS}"
if [[ -n "$MEM_OVERRIDE" ]]; then
    export BENCH_EXECUTOR_MEMORY="$MEM_OVERRIDE"
else
    BENCH_MEM_GB=$(( HOST_MEM_B / 1073741824 ))
    [[ "$BENCH_MEM_GB" -lt 4 ]] && BENCH_MEM_GB=8
    export BENCH_EXECUTOR_MEMORY="$(( BENCH_MEM_GB * 3 / 4 )).GB"
fi
export BENCH_TASK_CPUS="${TASK_CPUS_OVERRIDE:-}"
export BENCH_TASK_MEMORY="${TASK_MEM_OVERRIDE:-}"

BIND_CMD=""
BIND_STATUS="off"
if [[ "$BIND_MODE" != "off" ]]; then
    _last_core=$(( BENCH_EXECUTOR_CPUS - 1 ))
    if command -v taskset >/dev/null 2>&1 && [[ "$_last_core" -ge 0 ]]; then
        BIND_CMD="taskset -c 0-${_last_core}"
        BIND_STATUS="taskset:0-${_last_core}"
    else
        BIND_STATUS="unavailable:no-taskset"
    fi
    if [[ -z "$BIND_CMD" && "$BIND_MODE" == "on" ]]; then
        die "--bind requested but unavailable ($BIND_STATUS).
       Re-run with --no-bind to proceed unbound, and do not describe the result as
       an exclusive ${BENCH_EXECUTOR_CPUS}-core allocation."
    fi
fi
if [[ -n "$BIND_CMD" ]]; then
    echo "  cpu binding ON  ($BIND_STATUS) -- run is confined to $BENCH_EXECUTOR_CPUS cores"
else
    echo "  cpu binding OFF ($BIND_STATUS) -- executor.cpus is an ADMISSION limit only;"
    echo "                  tasks may exceed their request. Recorded in host.json."
fi
export BENCH_CPU_BINDING="$BIND_STATUS"

[[ -f "$REFERENCES" ]] || die "reference set not found: $REFERENCES"
if [[ -n "$EXTRA_CONFIG" ]]; then
    [[ -f "$EXTRA_CONFIG" ]] || die "extra config not found: $EXTRA_CONFIG"
    EXTRA_CONFIG="$(cd "$(dirname "$EXTRA_CONFIG")" && pwd)/$(basename "$EXTRA_CONFIG")"
    echo "  extra config $(basename "$EXTRA_CONFIG") (layered after benchmark.config)"
fi

_bad_refs=$(python3 - "$REFERENCES" <<'PY'
import glob, os, sys, yaml
refs = (yaml.safe_load(open(sys.argv[1])) or {}).get('input', {}).get('reference', {})
for key, path in refs.items():
    if not isinstance(path, str):
        continue
    if any(ch in path for ch in '*?['):
        if not glob.glob(path):
            print(f"    {key}: no files match {path}")
    elif not os.path.exists(path):
        print(f"    {key}: missing {path}")
PY
)
if [[ -n "$_bad_refs" ]]; then
    die "reference paths in $(basename "$REFERENCES") do not resolve:
$_bad_refs
       Edit that file, or pass a different set with --references."
fi
echo "  reference paths resolve"

_ref_star_cpus=$(python3 -c "import yaml,sys;print((yaml.safe_load(open(sys.argv[1])) or {}).get('_bench',{}).get('star_cpus',''))" "$REFERENCES" 2>/dev/null)
_ref_star_mem=$(python3 -c "import yaml,sys;print((yaml.safe_load(open(sys.argv[1])) or {}).get('_bench',{}).get('star_memory',''))" "$REFERENCES" 2>/dev/null)
_ref_label=$(python3 -c "import yaml,sys;print((yaml.safe_load(open(sys.argv[1])) or {}).get('_bench',{}).get('label','(unlabelled)'))" "$REFERENCES" 2>/dev/null)
export BENCH_STAR_CPUS="${STAR_CPUS_OVERRIDE:-$_ref_star_cpus}"
export BENCH_STAR_MEMORY="${STAR_MEM_OVERRIDE:-$_ref_star_mem}"
echo "  references  $_ref_label"
echo "              $(basename "$REFERENCES")"
[[ -n "$BENCH_STAR_CPUS$BENCH_STAR_MEMORY" ]] && \
    echo "  STAR sized  ${BENCH_STAR_CPUS:-8} cpus / ${BENCH_STAR_MEMORY:-48}.GB"
echo "  executor pinned to ${BENCH_EXECUTOR_CPUS} cpus / ${BENCH_EXECUTOR_MEMORY}"
[[ -n "$CPUS_OVERRIDE" ]] && echo "  (cpus set explicitly; host reports ${HOST_CPUS})"

_task_cpus="${BENCH_TASK_CPUS:-1}"
_task_mem="${BENCH_TASK_MEMORY:-2}"
_mem_gb="${BENCH_EXECUTOR_MEMORY%%.*}"
_min_mem=16                                   # RIBOPY_CREATE in fixed_resources.config
[[ -n "$BENCH_TASK_MEMORY" && "$BENCH_TASK_MEMORY" -gt "$_min_mem" ]] && _min_mem="$BENCH_TASK_MEMORY"
[[ -n "$BENCH_STAR_MEMORY" && "$BENCH_STAR_MEMORY" -gt "$_min_mem" ]] && _min_mem="$BENCH_STAR_MEMORY"
if [[ "$_mem_gb" -lt "$_min_mem" ]]; then
    die "executor memory ${BENCH_EXECUTOR_MEMORY} is below the largest single task
       (${_min_mem}.GB). Nextflow aborts rather than queueing. Raise --memory."
fi
_by_cpu=$(( BENCH_EXECUTOR_CPUS / _task_cpus )); [[ "$_by_cpu" -lt 1 ]] && _by_cpu=1
_by_mem=$(( _mem_gb / _task_mem ));             [[ "$_by_mem" -lt 1 ]] && _by_mem=1
_conc=$_by_cpu; [[ "$_by_mem" -lt "$_conc" ]] && _conc=$_by_mem
echo "  default-task concurrency ~${_conc}  (cpu-bound ${_by_cpu}, mem-bound ${_by_mem};" \
     "default task = ${_task_cpus} cpus / ${_task_mem}.GB)"

TIME_CMD=""
if [[ -x /usr/bin/time ]]; then
    if /usr/bin/time -v true >/dev/null 2>&1; then TIME_CMD="/usr/bin/time -v"      # GNU
    elif /usr/bin/time -l true >/dev/null 2>&1; then TIME_CMD="/usr/bin/time -l"    # BSD
    fi
fi
[[ -z "$TIME_CMD" ]] && echo "  note: no usable /usr/bin/time; driver-overhead stats skipped"

STAMP="$(date +%Y%m%d-%H%M%S)"
RESULTS="${OUTDIR:-$BENCH_DIR/results/$STAMP}"

declare -a IDS BASES OVERLAYS REPS DESCS
while IFS=$'\t' read -r id base overlay reps desc; do
    [[ -z "$id" || "$id" == \#* ]] && continue
    if [[ -n "$SCENARIO_FILTER" && ",$SCENARIO_FILTER," != *",$id,"* ]]; then continue; fi
    IDS+=("$id"); BASES+=("$base")
    OVERLAYS+=("$overlay"); REPS+=("${REPLICATES:-$reps}"); DESCS+=("$desc")
done < "$BENCH_DIR/scenarios.tsv"

[[ ${#IDS[@]} -gt 0 ]] || die "no scenarios selected"

MAXREP=0
for r in "${REPS[@]}"; do (( r > MAXREP )) && MAXREP=$r; done

declare -a ORDER_ID ORDER_REP
for (( rep=1; rep<=MAXREP; rep++ )); do
    for i in "${!IDS[@]}"; do
        (( rep <= REPS[i] )) && { ORDER_ID+=("$i"); ORDER_REP+=("$rep"); }
    done
done
TOTAL=${#ORDER_ID[@]}

log "Plan: ${#IDS[@]} scenarios, $TOTAL runs, profile '$PROFILE'"
for i in "${!IDS[@]}"; do
    printf '  %-22s x%s  %s\n' "${IDS[i]}" "${REPS[i]}" "${DESCS[i]}"
done
echo "  results -> $RESULTS"

if [[ $DRY_RUN -eq 1 ]]; then
    log "Dry run - commands that would execute"
    for k in $(seq 0 $((TOTAL-1))); do
        i=${ORDER_ID[$k]}; rep=${ORDER_REP[$k]}
        echo "  [$((k+1))/$TOTAL] nextflow run main.nf -profile $PROFILE \\"
        echo "        -params-file $RESULTS/runs/${IDS[i]}.rep${rep}/params.yaml \\"
        echo "        -c $BENCH_DIR/benchmark.config \\"
        [[ -n "$EXTRA_CONFIG" ]] && echo "        -c $EXTRA_CONFIG \\"
        echo "        -with-trace ... -with-report ... -with-timeline ..."
    done

    log "Dry run - params build + input check (one per scenario, replicates identical)"
    _dry=$(mktemp -d); _dry_bad=0
    for i in "${!IDS[@]}"; do
        overlay_arg=()
        [[ "${OVERLAYS[i]}" != "-" && -n "${OVERLAYS[i]}" ]] && \
            overlay_arg=(--overlay "$BENCH_DIR/${OVERLAYS[i]}")
        base_path="${BASES[i]}"
        [[ "$base_path" == scenarios/* ]] && base_path="$BENCH_DIR/$base_path"
        _d="$_dry/${IDS[i]}"; mkdir -p "$_d"
        if ! python3 "$BENCH_DIR/make_params.py" \
                --base "$base_path" ${overlay_arg[@]+"${overlay_arg[@]}"} \
                --references "$REFERENCES" \
                --run-dir "$_d" --out "$_d/params.yaml"; then
            echo "  ${IDS[i]}: params build FAILED"; _dry_bad=$((_dry_bad+1)); continue
        fi
        if ! write_inputs_json "$_d/params.yaml" "$_d/inputs.json"; then
            echo "  ${IDS[i]}: input check FAILED"; _dry_bad=$((_dry_bad+1)); continue
        fi
        python3 - "$_d/inputs.json" "${IDS[i]}" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
print(f"  {sys.argv[2]:<22} {d['n_files']} files "
      f"({d['n_ribo_files']} ribo + {d['n_rnaseq_files']} rnaseq), "
      f"{d['n_samples']} samples, {d['input_bytes']/2**30:.2f} GB")
PY
        echo "    params -> $_d/params.yaml"
    done
    [[ $_dry_bad -eq 0 ]] || die "$_dry_bad scenario(s) would not run; see above"
    echo
    echo "  merged params kept for inspection under $_dry"
    exit 0
fi

mkdir -p "$RESULTS/runs"

_j() { printf '%s' "$*" | tr -d '\n\r\t' | sed 's/\\/\\\\/g; s/"/\\"/g'; }

{
  echo "{"
  echo "  \"timestamp\": \"$(_j "$(date -u +%Y-%m-%dT%H:%M:%SZ)")\","
  echo "  \"cpu_model\": \"$(_j "$(awk -F': ' '/model name/{print $2; exit}' /proc/cpuinfo 2>/dev/null || echo unknown)")\","
  echo "  \"host_cpus\": $HOST_CPUS,"
  echo "  \"host_memory_bytes\": $HOST_MEM_B,"
  echo "  \"arch\": \"$(_j "$(uname -m)")\","
  echo "  \"os\": \"$(_j "$(uname -sr)")\","
  echo "  \"profile\": \"$PROFILE\","
  echo "  \"executor_cpus\": $BENCH_EXECUTOR_CPUS,"
  echo "  \"executor_memory\": \"$BENCH_EXECUTOR_MEMORY\","
  echo "  \"cpu_binding\": \"$(_j "${BENCH_CPU_BINDING:-off}")\","
  echo "  \"memory_enforced\": false,"
  echo "  \"task_cpus\": \"${BENCH_TASK_CPUS:-1}\","
  echo "  \"task_memory_gb\": \"${BENCH_TASK_MEMORY:-2}\","
  echo "  \"extra_config\": \"$(_j "$( [[ -n "$EXTRA_CONFIG" ]] && basename "$EXTRA_CONFIG" || echo none )")\","
  echo "  \"nextflow\": \"$(_j "$(nextflow -version 2>&1 | grep -oE 'version [0-9.]+' | head -1 | cut -d' ' -f2)")\","
  echo "  \"nxf_opts\": \"$NXF_OPTS\","
  echo "  \"references\": \"$(_j "$(basename "$REFERENCES")")\","
  echo "  \"references_label\": \"$_ref_label\","
  echo "  \"star_cpus\": \"${BENCH_STAR_CPUS:-8}\","
  echo "  \"star_memory_gb\": \"${BENCH_STAR_MEMORY:-48}\","
  echo "  \"git_pipeline\": \"$(_j "$(git rev-parse --short HEAD 2>/dev/null || echo unknown)")\","
  echo "  \"git_pipeline_dirty\": $( [[ -n "$(git status --porcelain 2>/dev/null)" ]] && echo true || echo false )"
  echo "}"
} > "$RESULTS/host.json"

SWEEP_START=$(date +%s)
FAILED=0

for k in $(seq 0 $((TOTAL-1))); do
    i=${ORDER_ID[$k]}; rep=${ORDER_REP[$k]}
    id="${IDS[i]}"; run_id="${id}.rep${rep}"
    run_dir="$RESULTS/runs/$run_id"
    mkdir -p "$run_dir"

    log "[$((k+1))/$TOTAL] $run_id"

    overlay_arg=()
    [[ "${OVERLAYS[i]}" != "-" && -n "${OVERLAYS[i]}" ]] && \
        overlay_arg=(--overlay "$BENCH_DIR/${OVERLAYS[i]}")
    base_path="${BASES[i]}"
    [[ "$base_path" == scenarios/* ]] && base_path="$BENCH_DIR/$base_path"
    python3 "$BENCH_DIR/make_params.py" \
        --base "$base_path" ${overlay_arg[@]+"${overlay_arg[@]}"} \
        --references "$REFERENCES" \
        --run-dir "$run_dir" --out "$run_dir/params.yaml" \
        || { echo "params build failed" >&2; FAILED=$((FAILED+1)); continue; }

    write_inputs_json "$run_dir/params.yaml" "$run_dir/inputs.json" \
        || { echo "input check failed" >&2; FAILED=$((FAILED+1)); continue; }

    rm -rf work .nextflow .nextflow.log* 2>/dev/null

    echo "{\"run_id\":\"$run_id\",\"scenario\":\"$id\",\"replicate\":$rep,\"run_order_index\":$k,\"description\":\"${DESCS[i]}\"}" \
        > "$run_dir/manifest.json"

    extra_c=""
    [[ -n "$EXTRA_CONFIG" ]] && extra_c="-c '$EXTRA_CONFIG'"

    run_start=$(date +%s)
    ${TIME_CMD:-} bash -c "
        ${BIND_CMD} nextflow run main.nf \
            -profile '$PROFILE' \
            -params-file '$run_dir/params.yaml' \
            -c '$BENCH_DIR/benchmark.config' \
            $extra_c \
            -with-trace    '$run_dir/trace.txt' \
            -with-report   '$run_dir/report.html' \
            -with-timeline '$run_dir/timeline.html' \
            < /dev/null > '$run_dir/nextflow.log' 2>&1
    " < /dev/null 2> "$run_dir/wrapper_time.txt"
    exit_code=$?
    run_end=$(date +%s)

    echo "$exit_code" > "$run_dir/exit_code"
    echo "$((run_end - run_start))" > "$run_dir/wall_seconds"

    if [[ $exit_code -ne 0 ]]; then
        echo "  FAILED (exit $exit_code) - see $run_dir/nextflow.log"
        FAILED=$((FAILED+1))
    else
        ntasks=$(( $(wc -l < "$run_dir/trace.txt" 2>/dev/null || echo 1) - 1 ))
        printf '  ok  %ss, %s tasks\n' "$((run_end - run_start))" "$ntasks"
    fi
done

log "Sweep finished in $(( ($(date +%s) - SWEEP_START) / 60 )) min ($FAILED failed of $TOTAL)"

log "Aggregating"
python3 "$BENCH_DIR/aggregate.py" "$RESULTS" || die "aggregation failed"
echo
echo "Report:  $RESULTS/REPORT.md"
echo "Tables:  $RESULTS/{table1,averages,summary,stages,reads,dedup_steps,efficiency}.csv"
[[ $FAILED -gt 0 ]] && exit 1
exit 0
