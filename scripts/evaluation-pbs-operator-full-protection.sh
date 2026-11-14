#!/usr/bin/env bash
set -euo pipefail

# Run the TFHE/GINX PBS full-protection experiment for all security levels.
#
# Usage:
#   scripts/evaluation-pbs-operator-full-protection.sh <benchmark> <result-directory> \
#       [trials] [pbs-samples-per-trial] [inter-trial-seconds] [parallelism]
#
# Positional arguments:
#   benchmark: evaluation-pbs-operator-full-protection executable built with
#       -DWITH_PBS_OPERATOR_BREAKDOWN=ON.
#   result-directory: destination for per-security logs and CSV summaries.
#   trials: independent process/context repetitions (default: 5).
#   pbs-samples-per-trial: measured PBS calls in one trial (default: 1000). The
#       executable also performs one unmeasured warm-up per implementation.
#   inter-trial-seconds: pause between trial batches, outside measurement
#       (default: 2).
#   parallelism: maximum concurrent trials within one security level (default:
#       2). All trials in one level finish before the next level starts.
#
# The script runs STD128, STD192, then STD256. Raw output is retained in
# <result-directory>/<security>/trial-*.log. Each security directory receives
# its own summary.csv with arithmetic means over exactly that level's logs.
# Full-Protection-PBS overhead is calculated from averaged total PBS times,
# rather than by averaging trial-level percentages.

if [[ $# -lt 2 || $# -gt 6 ]]; then
    echo "usage: $0 <benchmark> <result-directory> [trials] [pbs-samples-per-trial] [inter-trial-seconds] [parallelism]" >&2
    exit 2
fi

benchmark=$1
result_dir=$2
trials=${3:-5}
pbs_samples=${4:-1000}
inter_trial_seconds=${5:-2}
parallelism=${6:-2}

ValidateNonnegativeInteger() {
    local name=$1
    local value=$2
    if [[ ! $value =~ ^[0-9]+$ ]]; then
        echo "$name must be a nonnegative integer: $value" >&2
        exit 2
    fi
}

ValidatePositiveInteger() {
    local name=$1
    local value=$2
    ValidateNonnegativeInteger "$name" "$value"
    if ((10#$value == 0)); then
        echo "$name must be greater than zero" >&2
        exit 2
    fi
}

ValidatePositiveInteger "trials" "$trials"
ValidatePositiveInteger "pbs-samples-per-trial" "$pbs_samples"
ValidateNonnegativeInteger "inter-trial-seconds" "$inter_trial_seconds"
ValidatePositiveInteger "parallelism" "$parallelism"

if [[ ! -x "$benchmark" ]]; then
    echo "benchmark is not executable: $benchmark" >&2
    exit 2
fi

# Average a key=value field from exactly the logs created by this invocation.
AverageField() {
    local match_kind=$1
    local implementation=$2
    local operator=$3
    local field=$4
    shift 4

    awk -v kind="$match_kind" -v implementation="$implementation" -v operator="$operator" -v field="$field" \
        -v expected_trials="$trials" '
        {
            matched = 0
            if (kind == "parameter")
                matched = ($1 ~ /^parameter=/)
            else if (kind == "implementation")
                matched = ($1 == "implementation=" implementation && $2 == "operator=" operator)
            if (!matched)
                next

            for (i = 1; i <= NF; ++i) {
                split($i, pair, "=")
                if (pair[1] == field) {
                    if (pair[2] !~ /^-?[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$/) {
                        printf "invalid %s in %s: %s\\n", field, FILENAME, pair[2] > "/dev/stderr"
                        exit 2
                    }
                    sum += pair[2]
                    count += 1
                    break
                }
            }
        }
        END {
            if (count != expected_trials) {
                printf "expected %d %s values for %s/%s, found %d\\n", expected_trials, field, implementation, operator, count > "/dev/stderr"
                exit 2
            }
            printf "%.6f", sum / count
        }
    ' "$@"
}

mkdir -p "$result_dir"
progress_log="$result_dir/progress.log"

printf 'start=%s security_levels=STD128,STD192,STD256 trials=%d pbs_samples_per_trial=%d inter_trial_seconds=%d parallelism=%d\n' \
    "$(date --iso-8601=seconds)" "$trials" "$pbs_samples" "$inter_trial_seconds" "$parallelism" \
    | tee "$progress_log"

for security in STD128 STD192 STD256; do
    configuration_dir="$result_dir/$security"
    mkdir -p "$configuration_dir"
    summary_csv="$configuration_dir/summary.csv"
    printf 'Record,ParameterSet,Implementation,Operator,AverageTimeNs,RatioPct,OverheadPct\n' > "$summary_csv"
    trial_logs=()
    trial_pids=()
    trial_ids=()

    for trial in $(seq -w 1 "$trials"); do
        trial_log="$configuration_dir/trial-${trial}.log"
        trial_logs+=("$trial_log")
        printf '%s START security=%s trial=%s/%s\n' "$(date --iso-8601=seconds)" "$security" "$trial" "$trials" \
            | tee -a "$progress_log"
        "$benchmark" "$security" "$pbs_samples" >"$trial_log" 2>&1 &
        trial_pids+=("$!")
        trial_ids+=("$trial")

        if ((${#trial_pids[@]} < parallelism)) && [[ "$trial" != "$trials" ]]; then
            continue
        fi

        batch_failed=0
        batch_status=0
        for index in "${!trial_pids[@]}"; do
            batch_trial=${trial_ids[index]}
            batch_log="$configuration_dir/trial-${batch_trial}.log"
            if wait "${trial_pids[index]}"; then
                printf '%s DONE  security=%s trial=%s/%s\n' "$(date --iso-8601=seconds)" "$security" "$batch_trial" "$trials" \
                    | tee -a "$progress_log"
            else
                status=$?
                printf '%s FAIL  security=%s trial=%s/%s status=%d log=%s\n' "$(date --iso-8601=seconds)" "$security" \
                    "$batch_trial" "$trials" "$status" "$batch_log" | tee -a "$progress_log"
                batch_failed=1
                batch_status=$status
            fi
        done
        if ((batch_failed)); then
            exit "$batch_status"
        fi
        trial_pids=()
        trial_ids=()

        if [[ "$security" != "STD256" || "$trial" != "$trials" ]]; then
            printf '%s PAUSE %s seconds before next trial batch\n' "$(date --iso-8601=seconds)" "$inter_trial_seconds" \
                | tee -a "$progress_log"
            sleep "$inter_trial_seconds"
        fi
    done

    no_protection_time=$(AverageField parameter '' '' no_protection_average_time_ns "${trial_logs[@]}")
    full_protection_time=$(AverageField parameter '' '' full_protection_average_time_ns "${trial_logs[@]}")
    total_overhead=$(awk -v no_protection="$no_protection_time" -v full_protection="$full_protection_time" \
        'BEGIN { printf "%.6f", 100.0 * (full_protection - no_protection) / no_protection }')
    printf 'PBS total,%s,No-Protection-PBS,PBSTotal,%s,,\n' "$security" "$no_protection_time" >> "$summary_csv"
    printf 'PBS total,%s,Full-Protection-PBS,PBSTotal,%s,,\n' "$security" "$full_protection_time" >> "$summary_csv"
    printf 'PBS total,%s,Full-Protection-PBS-vs-No-Protection-PBS,PBSTotal,,,%s\n' "$security" "$total_overhead" \
        >> "$summary_csv"

    for implementation in No-Protection-PBS Full-Protection-PBS; do
        for operator in NTT INTT VCMM Decompose KeySwitch Other FiveOperators PBSTotal; do
            average_time=$(AverageField implementation "$implementation" "$operator" average_time_ns "${trial_logs[@]}")
            ratio=$(AverageField implementation "$implementation" "$operator" ratio_pct "${trial_logs[@]}")
            printf 'Breakdown,%s,%s,%s,%s,%s,\n' "$security" "$implementation" "$operator" "$average_time" "$ratio" \
                >> "$summary_csv"
        done
    done

    printf '%s SUMMARY security=%s path=%s\n' "$(date --iso-8601=seconds)" "$security" "$summary_csv" \
        | tee -a "$progress_log"
done

printf 'finish=%s status=success summaries=%s/<STD128|STD192|STD256>/summary.csv\n' "$(date --iso-8601=seconds)" "$result_dir" \
    | tee -a "$progress_log"
