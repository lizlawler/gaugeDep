#!/bin/zsh
#
# Local (zsh) script: computes both-angular joint-model BMA weights in parallel across all 200
# datasets per scenario using GNU parallel, with a per-job memory throttle (mem_throttle.sh) and
# a pre-flight single-job runtime benchmark. Set MAX_JOBS for your machine; DRYRUN=1 to preview.
#

MAX_JOBS=6	
DRYRUN=0  # Set to 1 for dry-run

dep_types=("gauss" "logistic")
dep_levels=("low" "mid" "high")
data_num=($(seq 1 200))

#####################################################
# PRE-FLIGHT RUNTIME ESTIMATE (single-job benchmark)
#####################################################

echo "Benchmarking one test job to estimate runtime..."
start=$(date +%s)
# Use a representative dataset (data_num=1)
Rscript --vanilla extraction_scripts/extract_weights_joint_both_ang.R gauss low 1
end=$(date +%s)
single_runtime=$(( end - start ))
echo "Single job took $single_runtime seconds (~$(( single_runtime/60 )) minutes)."

#####################################################
# GNU PARALLEL FLAGS (Zsh-safe)
#####################################################

RUNNER_FLAGS=(-j $MAX_JOBS --bar)
if (( DRYRUN == 1 )); then
    RUNNER_FLAGS+=("--dry-run")
    echo "\n>>> DRY RUN ENABLED — COMMANDS WILL NOT EXECUTE <<<\n"
fi

#####################################################
# PARALLEL RUNNER FUNCTION
#####################################################

run_parallel() {
    local cmd_template=$1
    local num_jobs=$2
    shift 2
    local estimated_time=$(( single_runtime * num_jobs / MAX_JOBS ))
    echo "Estimated runtime for this model ≈ $(( estimated_time / 60 )) minutes."

    # Each job first runs the mem_throttle.sh script
    parallel "${RUNNER_FLAGS[@]}" \
        './shell_scripts/local_machine/mem_throttle.sh && '"$cmd_template" \
        "$@"
}

#####################################################
# MODEL EXECUTION
#####################################################

TOTAL_JOBS=$(( ${#dep_types[@]} * ${#dep_levels[@]} * ${#data_num[@]} ))
echo "Running ($TOTAL_JOBS jobs)..."
run_parallel 'Rscript --vanilla extraction_scripts/extract_weights_joint_both_ang.R {1} {2} {3} >/dev/null 2>&1' \
    $TOTAL_JOBS \
    ::: "${dep_types[@]}" ::: "${dep_levels[@]}" ::: "${data_num[@]}"

echo "All $TOTAL_JOBS simulations completed."
