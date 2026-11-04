#!/usr/bin/env bash
# Run steps 02_01a (total TAS) and 02_01b (direct TAS) in parallel: one task per
# step x temperature variable (TS, TA) x site, each running the step's script for
# one site (TAS_SITE), as a Slurm array; then combine the tasks' results into the
# files the serial run of the scripts writes.
#
# From the repository root:
#   TAS_DIR_RAWDATA=<step-01 outputs> bash workflows/02_01_tas_parallel.sh submit
#   bash workflows/02_01_tas_parallel.sh combine     # once every task has finished
# Without Slurm:
#   TAS_DIR_RAWDATA=<...> bash workflows/02_01_tas_parallel.sh tasks   # write the task list
#   TAS_DIR_RAWDATA=<...> bash workflows/02_01_tas_parallel.sh run 3   # run task (line) 3
#
# Settings (environment variables):
#   TAS_DIR_RAWDATA     directory with RespirationData/ and the ERA5 soil-water file (required
#                       for tasks, submit and run)
#   TAS_OUT_DIR         results of the tasks, their logs and the task list (default _tas_parallel)
#   TAS_SITES           space-separated site IDs (default: every site in data/site_info.csv)
#   TAS_STEPS           "a b" (default): 02_01a and 02_01b
#   TAS_TEMP_VARS       "TS TA" (default)
#   TAS_PIXI_MANIFEST   pixi.toml of an environment with R, brms and cmdstan to activate in
#                       each task (e.g. refactor-v2's), with pixi from PATH, ~/.pixi/bin or PIXI;
#                       default: use the Rscript on PATH
#   TAS_PARTITION (day), TAS_TIME (23:30:00), TAS_MEM (8G), TAS_CPUS (4: the 4 chains of each
#   brms fit), TAS_MAX_CONCURRENT (200): the Slurm array's resources
#   TAS_COMBINE_DIR     where combine writes the combined tables (default data)
set -euo pipefail

repo_root="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$repo_root"
out_dir="${TAS_OUT_DIR:-_tas_parallel}"
mkdir -p "$out_dir"
out_dir="$(cd "$out_dir" && pwd)"
tasks_file="$out_dir/tasks.tsv"

script_of() {
  case "$1" in
    a) echo "workflows/02_01a_estimate_total_TAS_use_moving_window.R" ;;
    b) echo "workflows/02_01b_estimate_direct_TAS_use_moving_window.R" ;;
    *) echo "unknown step: $1" >&2; exit 1 ;;
  esac
}

activate_env() {
  if [[ -n "${TAS_PIXI_MANIFEST:-}" ]]; then
    # --as-is: use the installed environment, so that concurrent tasks never resolve or install it;
    # its activation scripts use unset variables
    set +u
    eval "$("${PIXI:-$(command -v pixi || echo "$HOME/.pixi/bin/pixi")}" shell-hook --as-is --manifest-path "$TAS_PIXI_MANIFEST")"
    set -u
  fi
}

require_rawdata() {
  if [[ -z "${TAS_DIR_RAWDATA:-}" || ! -d "$TAS_DIR_RAWDATA/RespirationData" ]]; then
    echo "TAS_DIR_RAWDATA must name a directory with RespirationData/ (got '${TAS_DIR_RAWDATA:-}')" >&2
    exit 1
  fi
}

write_tasks() {
  local sites
  if [[ -n "${TAS_SITES:-}" ]]; then
    sites="$TAS_SITES"
  else
    # first column of the site table, without its header (and BOM, quotes, CR)
    sites="$(tail -n +2 data/site_info.csv | cut -d, -f1 | tr -d '"\r')"
  fi
  : >| "$tasks_file"
  for step in ${TAS_STEPS:-a b}; do
    script_of "$step" > /dev/null
    for temp_var in ${TAS_TEMP_VARS:-TS TA}; do
      for site in $sites; do
        printf '%s\t%s\t%s\n' "$step" "$temp_var" "$site" >> "$tasks_file"
      done
    done
  done
  echo "$(wc -l < "$tasks_file" | tr -d ' ') tasks in $tasks_file"
}

run_task() {
  local n="${1:-${SLURM_ARRAY_TASK_ID:-}}"
  [[ -n "$n" ]] || { echo "run: give a task number (a line of $tasks_file)" >&2; exit 1; }
  local step temp_var site
  IFS=$'\t' read -r step temp_var site < <(sed -n "${n}p" "$tasks_file")
  [[ -n "${site:-}" ]] || { echo "no task $n in $tasks_file" >&2; exit 1; }
  local script task_dir
  script="$(script_of "$step")"
  task_dir="$out_dir/${temp_var}_${step}/$site"
  mkdir -p "$task_dir"
  # a failed earlier attempt must not leave results behind
  rm -f "$task_dir"/*.csv
  activate_env
  echo "task $n: step 02_01$step, $temp_var, $site -> $task_dir ($(hostname), $(date))"
  # each task its own plot file, instead of every task writing Rplots.pdf in the repository
  TAS_SITE="$site" TAS_TEMP_VAR="$temp_var" TAS_OUT_DIR="$task_dir" \
    Rscript -e "pdf(file.path('$task_dir', 'plots.pdf')); source('$script')"
  echo "task $n done ($(date))"
}

case "${1:-}" in
  tasks)
    require_rawdata
    write_tasks
    ;;
  submit)
    require_rawdata
    command -v sbatch > /dev/null || { echo "sbatch not found: submit on the cluster" >&2; exit 1; }
    write_tasks
    n_tasks="$(wc -l < "$tasks_file" | tr -d ' ')"
    mkdir -p "$out_dir/logs"
    sbatch --job-name=tas-02_01 \
      --array="1-${n_tasks}%${TAS_MAX_CONCURRENT:-200}" \
      --partition="${TAS_PARTITION:-day}" --time="${TAS_TIME:-23:30:00}" \
      --cpus-per-task="${TAS_CPUS:-4}" --mem="${TAS_MEM:-8G}" \
      --output="$out_dir/logs/%A_%a.log" \
      --export="ALL,TAS_OUT_DIR=$out_dir,TAS_DIR_RAWDATA=$TAS_DIR_RAWDATA,OPENBLAS_NUM_THREADS=1,OMP_NUM_THREADS=1" \
      "$repo_root/workflows/02_01_tas_parallel.sh" run
    ;;
  run)
    require_rawdata
    run_task "${2:-}"
    ;;
  combine)
    activate_env
    Rscript workflows/02_01_tas_parallel_combine.R "$tasks_file" "$out_dir" "${TAS_COMBINE_DIR:-data}"
    ;;
  *)
    sed -n '2,/^set -euo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'
    exit 1
    ;;
esac
