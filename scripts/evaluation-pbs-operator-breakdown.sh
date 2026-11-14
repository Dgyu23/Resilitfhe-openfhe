#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   scripts/evaluation-pbs-operator-breakdown.sh <benchmark> <result-directory> \
#       [trials] [warmups] [iterations] [inter-trial-seconds]
#
# Positional arguments:
#   benchmark: executable evaluation-pbs-operator-breakdown binary.
#   result-directory: directory where per-configuration trial logs are stored.
#   trials: independent process/context repetitions for each parameter setting
#       (default: 50). Each trial regenerates the secret key and bootstrapping
#       key, so this captures variation across independent setups.
#   warmups: PBS invocations before measurement in every trial (default: 50).
#       Warmups are not included in the reported timing and reduce cold-start
#       effects such as cache and runtime initialization.
#   iterations: measured PBS invocations in every trial (default: 1000). The
#       benchmark reports each operator's aggregate timing divided by this value.
#   inter-trial-seconds: pause after each trial, except the final one (default:
#       5). This delay is outside the measurement and can reduce thermal or
#       frequency carry-over between trials.
#
# The script runs all nine configurations: STD128/STD192/STD256 crossed with
# plaintext moduli p=2, p=4, and p=8. For p=8, the benchmark selects STDxxx_4.
# After all trials of each configuration finish, summary.log receives one CSV
# row per operator. Its Ratio(%) value is the arithmetic mean of that
# operator's ratio across the configuration's successful trial logs. A final
# Five operators row records the mean combined ratio of the five categories.

if [[ $# -lt 2 || $# -gt 6 ]]; then
    echo "usage: $0 <benchmark> <result-directory> [trials] [warmups] [iterations] [inter-trial-seconds]" >&2
    exit 2
fi

benchmark=$1
result_dir=$2
trials=${3:-50}
warmups=${4:-50}
iterations=${5:-1000}
inter_trial_seconds=${6:-5}

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
ValidateNonnegativeInteger "warmups" "$warmups"
ValidatePositiveInteger "iterations" "$iterations"
ValidateNonnegativeInteger "inter-trial-seconds" "$inter_trial_seconds"

# Extracts the mean Ratio(%) for one operator from all trial logs in a
# configuration directory. The benchmark prints the ratio as the final field
# of each operator line (for example, "42.31%").
AverageOperatorRatio() {
    local operator_pattern=$1
    local configuration_dir=$2

    awk -v pattern="$operator_pattern" -v expected_trials="$trials" '
        $0 ~ pattern {
            ratio = $NF
            sub(/%$/, "", ratio)
            if (ratio !~ /^[0-9]+([.][0-9]+)?$/) {
                printf "invalid operator ratio in %s: %s\n", FILENAME, $NF > "/dev/stderr"
                exit 2
            }
            sum += ratio
            count += 1
        }
        END {
            if (count != expected_trials) {
                printf "expected %d operator ratios matching %s, found %d\n", expected_trials, pattern, count > "/dev/stderr"
                exit 2
            }
            printf "%.6f", sum / count
        }
    ' "$configuration_dir"/trial-*.log
}

if [[ ! -x "$benchmark" ]]; then
    echo "benchmark is not executable: $benchmark" >&2
    exit 2
fi

mkdir -p "$result_dir"
progress_log="$result_dir/progress.log"
summary_log="$result_dir/summary.log"
printf 'Operator,PlaintextModulus,ParameterSet,Ratio(%%)\n' > "$summary_log"
printf 'start=%s trials=%d warmups=%d measured_pbs=%d inter_trial_seconds=%d\n' \
    "$(date --iso-8601=seconds)" "$trials" "$warmups" "$iterations" "$inter_trial_seconds" | tee "$progress_log"

for security in STD128 STD192 STD256; do
    for plaintext_modulus in 2 4 8; do
        configuration="${security}_p${plaintext_modulus}"
        configuration_dir="$result_dir/$configuration"
        mkdir -p "$configuration_dir"

        for trial in $(seq -w 1 "$trials"); do
            trial_log="$configuration_dir/trial-${trial}.log"
            printf '%s START %s trial=%s/%s\n' "$(date --iso-8601=seconds)" "$configuration" "$trial" "$trials" \
                | tee -a "$progress_log"

            if "$benchmark" "$security" "$plaintext_modulus" "$warmups" "$iterations" >"$trial_log" 2>&1; then
                printf '%s DONE  %s trial=%s/%s\n' "$(date --iso-8601=seconds)" "$configuration" "$trial" "$trials" \
                    | tee -a "$progress_log"
            else
                status=$?
                printf '%s FAIL  %s trial=%s/%s status=%d log=%s\n' "$(date --iso-8601=seconds)" "$configuration" \
                    "$trial" "$trials" "$status" "$trial_log" | tee -a "$progress_log"
                exit "$status"
            fi

            if [[ "$security" != "STD256" || "$plaintext_modulus" != "8" || "$trial" != "$trials" ]]; then
                printf '%s PAUSE %s seconds before next trial\n' "$(date --iso-8601=seconds)" "$inter_trial_seconds" \
                    | tee -a "$progress_log"
                sleep "$inter_trial_seconds"
            fi
        done

        parameter_set=$security
        if [[ "$plaintext_modulus" == "8" ]]; then
            parameter_set="${security}_4"
        fi

        # Keep this order aligned with the five PBS categories reported by the
        # benchmark and append configurations in the nested-loop order above.
        while IFS='|' read -r operator operator_pattern; do
            average_ratio=$(AverageOperatorRatio "$operator_pattern" "$configuration_dir")
            printf '%s,%s,%s,%s\n' "$operator" "$plaintext_modulus" "$parameter_set" "$average_ratio" \
                >> "$summary_log"
        done <<'OPERATORS'
NTT|^NTT[[:space:]]
INTT|^INTT[[:space:]]
VCMM|^VCMM[[:space:]]
Decompose|^Decompose[[:space:]]
KeySwitch|^KeySwitch[[:space:]]
Five operators|^Five[[:space:]]+operators[[:space:]]
OPERATORS
    done
done

printf 'finish=%s status=success\n' "$(date --iso-8601=seconds)" | tee -a "$progress_log"
printf 'summary=%s\n' "$summary_log" | tee -a "$progress_log"
