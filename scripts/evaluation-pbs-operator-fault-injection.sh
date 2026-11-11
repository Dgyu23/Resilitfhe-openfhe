#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   scripts/evaluation-pbs-operator-fault-injection.sh <benchmark> <result-directory> \
#       [trials] [iterations] [faults] [parallel]
#
# Positional arguments:
#   benchmark: executable evaluation-pbs-operator-fault-injection binary.
#   result-directory: root directory for per-operator, per-configuration
#       trial logs, progress.log, and summary.log.
#   trials: independent process/context repetitions for every
#       (operator, security parameter, plaintext modulus) configuration
#       (default: 10). Each trial constructs a new crypto context, secret key,
#       bootstrapping key, and switching key inside the benchmark.
#   iterations: fault-injected PBS invocations for each fault site in one
#       trial (default: 10000). For example, VCMM has four sites, so one VCMM
#       trial executes 4 * iterations injected PBS calls.
#   faults: number of independently selected transient bit flips per PBS;
#       supported values are 1 (default) and 2. When faults=2, the benchmark
#       selects two distinct dynamic locations within the selected site.
#   parallel: maximum number of independent trial processes running at once
#       (default: 4). It is process-level parallelism; OpenMP remains disabled
#       by the build configuration.
#
# The script creates all 45 configuration/trial tasks in this fixed order:
# VCMM, NTT, INTT, Decompose, KeySwitch; then STD128/STD192/STD256; then
# p=2/p=4/p=8; then trial number. A dynamic task pool keeps at most "parallel"
# jobs active and launches the next queued task whenever any job finishes.
# For p=8, the benchmark selects the corresponding STDxxx_4 parameter set.
#
# summary.log is written only after every task succeeds, so it never races an
# unfinished trial log. Its CSV rows follow the same operator/security/p order.
# For one row, the script first averages each fault site's DER over "trials",
# then takes the arithmetic mean across that operator's sites. Individual site
# DERs and live progress are retained in trial-*.log and progress.log.
if [[ $# -lt 2 || $# -gt 6 ]]; then
    echo "usage: $0 <benchmark> <result-directory> [trials=10] [iterations=10] [faults=1] [parallel=4]" >&2
    exit 2
fi
benchmark=$1
result_dir=$2
trials=${3:-10}
iterations=${4:-10000}
faults=${5:-1}
parallel=${6:-4}
# Change this value here when a different deterministic random sequence is wanted.
base_seed=7

positive_integer() {
    [[ $2 =~ ^[0-9]+$ ]] && ((10#$2 > 0)) || { echo "$1 must be a positive integer: $2" >&2; exit 2; }
}
positive_integer trials "$trials"
positive_integer iterations "$iterations"
positive_integer parallel "$parallel"
[[ $faults == 1 || $faults == 2 ]] || { echo "faults must be 1 or 2" >&2; exit 2; }
[[ -x $benchmark ]] || { echo "benchmark is not executable: $benchmark" >&2; exit 2; }

# Extracts one configuration's operator DER from all completed trial logs.
# Each benchmark trial prints one "site=... error_rate=..." line per site. The
# result gives every site equal weight, rather than weighting a site by its
# number of dynamic injection candidates.
average_der() {
    awk -F'[ =]' -v expected_trials="$trials" '
        /error_rate=/ {
            site = ""
            rate = ""
            for (i = 1; i <= NF; ++i) {
                if ($i == "site") site = $(i + 1)
                if ($i == "error_rate") rate = $(i + 1)
            }
            if (site != "" && rate != "") {
                sum[site] += rate * 100.0
                count[site] += 1
            }
        }
        END {
            for (site in count) {
                if (count[site] != expected_trials) exit 1
                operator_sum += sum[site] / count[site]
                site_count += 1
            }
            if (site_count == 0) exit 1
            printf "%.6f", operator_sum / site_count
        }
    ' "$1"/trial-*.log
}

mkdir -p "$result_dir"
summary_log="$result_dir/summary.log"
progress_log="$result_dir/progress.log"
printf 'start=%s trials=%s iterations_per_site=%s faults=%s parallel=%s base_seed=%s\n' \
    "$(date --iso-8601=seconds)" "$trials" "$iterations" "$faults" "$parallel" "$base_seed" > "$progress_log"

# Required CSV order. Site remains a dimension in individual logs but not in
# summary.log, whose DER is the equal-weight mean described above.
operators=(
    'VCMM vcmm'
    'NTT ntt'
    'INTT intt'
    'Decompose decompose'
    'KeySwitch keyswitch'
)

task_operator_labels=()
task_operators=()
task_securities=()
task_parameter_sets=()
task_plaintext_moduli=()
task_trials=()
task_seeds=()
task_logs=()

# Reap one completed task. A failed child is recorded but does not prevent
# already-running tasks from finishing; summary.log is withheld afterwards.
wait_for_one() {
    local status=0
    if wait -n; then
        :
    else
        status=$?
        failed_jobs=$((failed_jobs + 1))
        printf '%s TASK_FAILED exit_status=%s failed_jobs=%s\n' \
            "$(date --iso-8601=seconds)" "$status" "$failed_jobs" >> "$progress_log"
    fi
    running_jobs=$((running_jobs - 1))
}

# Build the complete task table before starting any benchmark. This keeps
# parallelism useful even when trials=1: different operators/parameters can
# occupy the same batch, and the construction order defines reproducible seeds.
task_index=0
for entry in "${operators[@]}"; do
    read -r operator_label operator <<< "$entry"
    for security in STD128 STD192 STD256; do
        for plaintext_modulus in 2 4 8; do
            parameter_set=$security
            [[ $plaintext_modulus == 8 ]] && parameter_set="${security}_4"
            config_dir="$result_dir/${operator_label}/${security}_p${plaintext_modulus}"
            mkdir -p "$config_dir"
            for trial in $(seq -w 1 "$trials"); do
                task_operator_labels+=("$operator_label")
                task_operators+=("$operator")
                task_securities+=("$security")
                task_parameter_sets+=("$parameter_set")
                task_plaintext_moduli+=("$plaintext_modulus")
                task_trials+=("$trial")
                task_seeds+=("$((base_seed + task_index))")
                task_logs+=("$config_dir/trial-${trial}.log")
                task_index=$((task_index + 1))
            done
        done
    done
done

printf 'tasks=%s\n' "${#task_operators[@]}" >> "$progress_log"
running_jobs=0
failed_jobs=0
for task_index in "${!task_operators[@]}"; do
    (
        printf '%s START task=%s/%s operator=%s parameter_set=%s p=%s trial=%s/%s\n' \
            "$(date --iso-8601=seconds)" "$((task_index + 1))" "${#task_operators[@]}" \
            "${task_operator_labels[task_index]}" "${task_parameter_sets[task_index]}" \
            "${task_plaintext_moduli[task_index]}" "${task_trials[task_index]}" "$trials" >> "$progress_log"
        if "$benchmark" "${task_operators[task_index]}" all "${task_securities[task_index]}" \
                "${task_plaintext_moduli[task_index]}" "$iterations" "$faults" "${task_seeds[task_index]}" \
                > "${task_logs[task_index]}" 2>&1; then
            :
        else
            status=$?
            printf '%s FAIL task=%s/%s operator=%s parameter_set=%s p=%s trial=%s/%s exit_status=%s log=%s\n' \
                "$(date --iso-8601=seconds)" "$((task_index + 1))" "${#task_operators[@]}" \
                "${task_operator_labels[task_index]}" "${task_parameter_sets[task_index]}" \
                "${task_plaintext_moduli[task_index]}" "${task_trials[task_index]}" "$trials" "$status" \
                "${task_logs[task_index]}" >> "$progress_log"
            exit "$status"
        fi
        printf '%s DONE task=%s/%s operator=%s parameter_set=%s p=%s trial=%s/%s log=%s\n' \
            "$(date --iso-8601=seconds)" "$((task_index + 1))" "${#task_operators[@]}" \
            "${task_operator_labels[task_index]}" "${task_parameter_sets[task_index]}" \
            "${task_plaintext_moduli[task_index]}" "${task_trials[task_index]}" "$trials" \
            "${task_logs[task_index]}" >> "$progress_log"
    ) &
    running_jobs=$((running_jobs + 1))
    if ((running_jobs >= parallel)); then
        wait_for_one
    fi
done
while ((running_jobs > 0)); do
    wait_for_one
done

if ((failed_jobs > 0)); then
    printf 'aborted=%s failed_jobs=%s summary_not_updated=%s\n' \
        "$(date --iso-8601=seconds)" "$failed_jobs" "$summary_log" >> "$progress_log"
    exit 1
fi

# All background trials have completed, so every configuration can now be
# scanned without racing an unfinished log file. The loop also deliberately
# preserves the CSV order documented at the top of this script.
printf 'Operator,ParameterSet,PlaintextModulus,Decryption Error Rate (%%)\n' > "$summary_log"
for entry in "${operators[@]}"; do
    read -r operator_label _ <<< "$entry"
    for security in STD128 STD192 STD256; do
        for plaintext_modulus in 2 4 8; do
            parameter_set=$security
            [[ $plaintext_modulus == 8 ]] && parameter_set="${security}_4"
            config_dir="$result_dir/${operator_label}/${security}_p${plaintext_modulus}"
            der=$(average_der "$config_dir")
            printf '%s,%s,%s,%s\n' "$operator_label" "$parameter_set" "$plaintext_modulus" "$der" >> "$summary_log"
        done
    done
done
printf 'finished=%s summary=%s\n' "$(date --iso-8601=seconds)" "$summary_log" >> "$progress_log"
