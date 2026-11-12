#!/usr/bin/env bash
set -euo pipefail

# Run one protected BinFHE operator for STD128, STD192, and STD256.  Each
# security level is executed once, using the same benchmark process for all
# methods of that level.  summary.log records the average nanoseconds per
# invocation and the relative overhead of every protected method.
#
# Usage:
#   scripts/evaluation-operator-protection.sh [operation] [warmup] [measure]
#
# Positional arguments:
#   operation: operator to benchmark: ntt, intt, decompose, keyswitch, or
#       vcmm (default: ntt).
#   warmup: invocations of every method before timing. Warmups are excluded
#       from the reported time and reduce cache and runtime initialization
#       effects (default: 10000; keyswitch default: 100).
#   measure: timed invocations of every method. The benchmark divides each
#       method's aggregate timing by this value (default: 100000; keyswitch
#       default: 5000).
#
# For each security level, the script writes the benchmark output to
# results/operator-protection/<OPERATOR>/<operation>-<security>.log.  The
# corresponding summary.log has one CSV row per method, with No-Protection as
# the overhead baseline. progress.log records start, completion, and output
# locations. A later invocation for the same operator replaces these files.
#
# Environment overrides:
#   EVALUATION_OPERATOR_PROTECTION_BIN: path to the benchmark executable
#       (default: build-release/bin/examples/binfhe/evaluation-operator-protection).
#   EVALUATION_OPERATOR_PROTECTION_DIR: root directory for operator results
#       (default: results/operator-protection).

validate_nonnegative_integer() {
    local name=$1
    local value=$2
    if [[ ! $value =~ ^[0-9]+$ ]]; then
        echo "$name must be a nonnegative integer: $value" >&2
        exit 2
    fi
}

validate_positive_integer() {
    validate_nonnegative_integer "$1" "$2"
    if ((10#$2 == 0)); then
        echo "$1 must be greater than zero" >&2
        exit 2
    fi
}

operation=${1:-ntt}
if [[ $operation != ntt && $operation != intt && $operation != decompose && $operation != keyswitch && $operation != vcmm ]]; then
    echo "operation must be ntt, intt, decompose, keyswitch, or vcmm: $operation" >&2
    exit 2
fi
if [[ $operation == keyswitch ]]; then
    warmup=${2:-100}
    measure=${3:-5000}
else
    warmup=${2:-10000}
    measure=${3:-100000}
fi
validate_nonnegative_integer "warmup" "$warmup"
validate_positive_integer "measure" "$measure"

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$script_dir/.." && pwd)
benchmark=${EVALUATION_OPERATOR_PROTECTION_BIN:-"$repo_root/build-release/bin/examples/binfhe/evaluation-operator-protection"}
result_root=${EVALUATION_OPERATOR_PROTECTION_DIR:-"$repo_root/results/operator-protection"}
operator_label=${operation^^}
summary_operator=$operator_label
if [[ $operation == decompose || $operation == keyswitch || $operation == vcmm ]]; then
    # The operator names used by the benchmark and paper tables.
    summary_operator=$operation
fi
result_dir="$result_root/$operator_label"

if [[ ! -x $benchmark ]]; then
    echo "benchmark is not executable: $benchmark" >&2
    exit 2
fi

mkdir -p "$result_dir"
summary_log="$result_dir/summary.log"
progress_log="$result_dir/progress.log"
printf 'Operator,Method,Params,Time,ProtectionOverhead(%%)\n' > "$summary_log"
printf 'start=%s warmup=%s measure=%s benchmark=%s\n' \
    "$(date --iso-8601=seconds)" "$warmup" "$measure" "$benchmark" > "$progress_log"

# Extract exactly one method's benchmark-reported average from a per-security
# log. The executable emits the value in nanoseconds.
extract_average_time_ns() {
    local log=$1
    local method=$2
    awk '
        index($0, "Operator=" operation ", Method=" requested ",") == 1 && /Average time: / {
            split($0, fields, "Average time: ")
            split(fields[2], value, " ")
            if (value[1] !~ /^[0-9]+([.][0-9]+)?$/) {
                printf "invalid average time in %s: %s\n", FILENAME, value[1] > "/dev/stderr"
                exit 2
            }
            time = value[1]
            count += 1
        }
        END {
            if (count != 1) {
                printf "expected one average time in %s, found %d\n", FILENAME, count > "/dev/stderr"
                exit 2
            }
            print time
        }
    ' operation="$operation" requested="$method" "$log"
}

# Calculate the percentage change from the unprotected baseline.
protection_overhead() {
    awk -v baseline="$1" -v protected="$2" 'BEGIN {
        if (baseline <= 0) {
            print "baseline time must be positive" > "/dev/stderr"
            exit 2
        }
        printf "%.6f", 100.0 * (protected - baseline) / baseline
    }'
}

# Run all methods for one security level and retain its complete stdout/stderr
# in a separate log, so summary extraction never mixes parameter sets.
run_security() {
    local security=$1
    local log="$result_dir/${operation}-${security}.log"

    printf '%s START security=%s methods=%s\n' \
        "$(date --iso-8601=seconds)" "$security" "${methods[*]}" \
        | tee -a "$progress_log" >&2
    "$benchmark" "$operation" "$security" "$warmup" "$measure" > "$log" 2>&1
    printf '%s DONE  security=%s log=%s\n' "$(date --iso-8601=seconds)" "$security" "$log" \
        | tee -a "$progress_log" >&2
}

# Keep this order consistent with the methods emitted by the executable and
# the order of rows in summary.log.
if [[ $operation == decompose || $operation == keyswitch ]]; then
    methods=("No-Protection" "DMR" "Ours")
elif [[ $operation == vcmm ]]; then
    methods=("No-Protection" "DMR" "Naive-ABFT" "Ours")
else
    methods=("No-Protection" "DMR" "Naive-ABFT" "Ours")
fi

for security in STD128 STD192 STD256; do
    run_security "$security"
    log="$result_dir/${operation}-${security}.log"
    baseline_time=$(extract_average_time_ns "$log" "No-Protection")
    printf '%s,No-Protection,%s,%s,-\n' "$summary_operator" "$security" "$baseline_time" >> "$summary_log"

    # Protection overhead is always relative to the corresponding
    # No-Protection result for this security level.
    for method in "${methods[@]:1}"; do
        method_time=$(extract_average_time_ns "$log" "$method")
        method_overhead=$(protection_overhead "$baseline_time" "$method_time")
        printf '%s,%s,%s,%s,%s\n' "$summary_operator" "$method" "$security" "$method_time" "$method_overhead" >> "$summary_log"
    done
done

printf 'finish=%s status=success summary=%s\n' "$(date --iso-8601=seconds)" "$summary_log" | tee -a "$progress_log"
