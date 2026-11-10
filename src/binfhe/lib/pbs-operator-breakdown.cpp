#include "pbs-operator-breakdown.h"

namespace lbcrypto {
namespace {
thread_local bool pbsOperatorBreakdownEnabled = false;
thread_local PBSOperatorBreakdownData pbsOperatorBreakdownData;
}  // namespace

void PBSOperatorBreakdown::Enable(bool enabled) {
    pbsOperatorBreakdownEnabled = enabled;
}

bool PBSOperatorBreakdown::IsEnabled() {
    return pbsOperatorBreakdownEnabled;
}

void PBSOperatorBreakdown::Reset() {
    pbsOperatorBreakdownData = {};
}

void PBSOperatorBreakdown::Add(PBSOperatorBreakdownStage stage, uint64_t nanoseconds) {
    pbsOperatorBreakdownData.nanoseconds[static_cast<std::size_t>(stage)] += nanoseconds;
}

PBSOperatorBreakdownData PBSOperatorBreakdown::Snapshot() {
    return pbsOperatorBreakdownData;
}

}  // namespace lbcrypto
