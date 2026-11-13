#!/usr/bin/env bash
set -euo pipefail

# Run every reliability configuration for one protected BinFHE operator.
#
# Usage:
#   scripts/evaluation-operator-reliability.sh <benchmark> <result-directory> \
#       <operator> [samples] [parallel]
#
# The script expands the selected operator into all of its current fault sites,
# STD128/STD192/STD256, and both 1F and 2F.  The executable's default seed (7)
# is used for every configuration.  Complete per-configuration output is kept
# under <result-directory>/<OPERATOR>, beside a CSV summary.
#
# Positional arguments:
#   benchmark: evaluation-operator-reliability executable.
#   result-directory: root directory for this experiment's results.
#   operator: ntt, intt, vcmm, decompose, or keyswitch.
#   samples: fault-injection samples per configuration (default: 1000000).
#   parallel: maximum concurrent benchmark processes (default: 1).

if [[ $# -lt 3 || $# -gt 5 ]]; then
    echo "usage: $0 <benchmark> <result-directory> <ntt|intt|vcmm|decompose|keyswitch> [samples] [parallel]" >&2
    exit 2
fi

benchmark=$1
result_root=$2
operation=$3
samples=${4:-1000000}
parallel=${5:-1}
if [[ ! $samples =~ ^[0-9]+$ ]] || ((10#$samples == 0)); then
    echo "samples must be a positive integer: $samples" >&2
    exit 2
fi

case $operation in
    ntt)
        sites=(NTT-butterfly-load NTT-butterfly-add NTT-butterfly-sub NTT-butterfly-mul NTT-butterfly-store)
        ;;
    intt)
        sites=(INTT-butterfly-load INTT-butterfly-add INTT-butterfly-sub INTT-butterfly-mul INTT-butterfly-store)
        ;;
    vcmm)
        sites=(input multiplication accumulation output)
        ;;
    decompose)
        sites=(digit-extraction intermediate-digit)
        ;;
    keyswitch)
        sites=(key-lookup accumulation update)
        ;;
    *)
        echo "operator must be ntt, intt, vcmm, decompose, or keyswitch: $operation" >&2
        exit 2
        ;;
esac

if [[ ! $parallel =~ ^[0-9]+$ ]] || ((10#$parallel == 0)); then
    echo "parallel must be a positive integer: $parallel" >&2
    exit 2
fi
if [[ ! -x $benchmark ]]; then
    echo "benchmark is not executable: $benchmark" >&2
    exit 2
fi

result_dir="$result_root/${operation^^}"
log_dir="$result_dir/logs"
summary_log="$result_dir/summary.csv"
progress_log="$result_dir/progress.log"
mkdir -p "$log_dir"

printf 'start=%s operator=%s samples=%s parallel=%s seed=default\n' \
    "$(date --iso-8601=seconds)" "$operation" "$samples" "$parallel" > "$progress_log"

task_sites=()
task_securities=()
task_faults=()
task_logs=()
for security in STD128 STD192 STD256; do
    for site in "${sites[@]}"; do
        for faults in 1 2; do
            task_sites+=("$site")
            task_securities+=("$security")
            task_faults+=("$faults")
            task_logs+=("$log_dir/${operation}-${security}-${site}-${faults}F.log")
        done
    done
done

running_jobs=0
failed_jobs=0
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

for task_index in "${!task_sites[@]}"; do
    (
        site=${task_sites[task_index]}
        security=${task_securities[task_index]}
        faults=${task_faults[task_index]}
        log=${task_logs[task_index]}
        printf '%s START task=%s/%s site=%s security=%s faults=%s\n' \
            "$(date --iso-8601=seconds)" "$((task_index + 1))" "${#task_sites[@]}" \
            "$site" "$security" "$faults" >> "$progress_log"
        if "$benchmark" "$operation" "$site" "$security" "$samples" "$faults" > "$log" 2>&1; then
            printf '%s DONE task=%s/%s site=%s security=%s faults=%s log=%s\n' \
                "$(date --iso-8601=seconds)" "$((task_index + 1))" "${#task_sites[@]}" \
                "$site" "$security" "$faults" "$log" >> "$progress_log"
        else
            status=$?
            printf '%s FAIL task=%s/%s site=%s security=%s faults=%s exit_status=%s log=%s\n' \
                "$(date --iso-8601=seconds)" "$((task_index + 1))" "${#task_sites[@]}" \
                "$site" "$security" "$faults" "$status" "$log" >> "$progress_log"
            exit "$status"
        fi
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
    printf 'finish=%s status=failed failed_jobs=%s summary_not_updated=%s\n' \
        "$(date --iso-8601=seconds)" "$failed_jobs" "$summary_log" >> "$progress_log"
    exit 1
fi

printf 'Operator,ParameterSet,Site,FaultsPerSample,Samples,DetectedSamples,UndetectedSamples,MaskedSamples,DetectionRate(%%),MaskingRate(%%)\n' \
    > "$summary_log"
for task_index in "${!task_sites[@]}"; do
    log=${task_logs[task_index]}
    awk '
        BEGIN { count = 0 }
        /^operator=/ {
            delete value
            for (i = 1; i <= NF; ++i) {
                split($i, field, "=")
                value[field[1]] = field[2]
            }
            required = (value["operator"] != "" && value["parameter"] != "" && value["site"] != "" &&
                        value["faults_per_sample"] != "" && value["samples"] != "" &&
                        value["detected_samples"] != "" && value["undetected_samples"] != "" &&
                        value["masked_samples"] != "" && value["detection_rate_pct"] != "" &&
                        value["masking_rate_pct"] != "")
            if (!required) exit 2
            printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", value["operator"], value["parameter"],
                   value["site"], value["faults_per_sample"], value["samples"], value["detected_samples"],
                   value["undetected_samples"], value["masked_samples"], value["detection_rate_pct"],
                   value["masking_rate_pct"]
            count += 1
        }
        END { if (count != 1) exit 2 }
    ' "$log" >> "$summary_log"
done

printf 'finish=%s status=success summary=%s\n' \
    "$(date --iso-8601=seconds)" "$summary_log" >> "$progress_log"
echo "completed=${#task_sites[@]} summary=$summary_log"
