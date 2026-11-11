#include "utils/pbs-sbf-fault-injection.h"

#include <stdexcept>

namespace lbcrypto {
namespace detail {
// Each benchmark trial is a separate process and the SBF experiment is built
// with WITH_OPENMP=OFF.  These are deliberately ordinary process-local
// variables: Match() is on the PBS hot path and atomic RMW operations at every
// candidate event dominate the measured execution time.
bool g_configured{false};
PBSFaultSite g_site{PBSFaultSite::BUTTERFLY_LOAD};
PBSFaultTransform g_transform{PBSFaultTransform::ANY};
PBSFaultTransform g_currentTransform{PBSFaultTransform::ANY};
std::vector<PBSFaultTarget> g_targets;
uint64_t g_seen{0};
uint64_t g_injected{0};
uint32_t g_observedBitWidth{0};
bool g_active{false};
bool g_blockConfigured{false};
PBSFaultBlockKind g_selectedBlockKind{PBSFaultBlockKind::NONE};
uint32_t g_selectedEvalAcc{0};
uint32_t g_selectedSubBlock{0};
PBSFaultBlockKind g_currentBlockKind{PBSFaultBlockKind::NONE};
uint32_t g_currentEvalAcc{0};
uint32_t g_currentSubBlock{0};
bool g_preExploring{false};
std::vector<PBSFaultBlock> g_preExploreBlocks;
std::array<uint64_t, 42> g_preExploreOffsets{};
}  // namespace detail

void PBSFaultInjection::Configure(PBSFaultSite site, std::vector<PBSFaultTarget> targets, PBSFaultTransform transform) {
    if (targets.empty() || targets.size() > 2)
        throw std::invalid_argument("PBS SBF plan must contain one or two targets");
    std::sort(targets.begin(), targets.end(), [](const auto& a, const auto& b) { return a.occurrence < b.occurrence; });
    for (size_t i = 1; i < targets.size(); ++i)
        if (targets[i - 1].occurrence == targets[i].occurrence)
            throw std::invalid_argument("PBS SBF targets must be distinct");
    detail::g_configured = false;
    detail::g_blockConfigured = false;
    detail::g_site = site;
    detail::g_transform = transform;
    detail::g_targets = std::move(targets);
    detail::g_seen = 0;
    detail::g_injected = 0;
    detail::g_observedBitWidth = 0;
    detail::g_configured = true;
}

void PBSFaultInjection::ConfigureBlock(PBSFaultSite site, std::vector<PBSFaultTarget> targets,
                                       PBSFaultTransform transform, PBSFaultBlockKind blockKind,
                                       uint32_t evalAcc, uint32_t subBlock) {
    Configure(site, std::move(targets), transform);
    detail::g_selectedBlockKind = blockKind;
    detail::g_selectedEvalAcc = evalAcc;
    detail::g_selectedSubBlock = subBlock;
    detail::g_blockConfigured = true;
}

void PBSFaultInjection::Disable() {
    detail::g_configured = false;
    detail::g_blockConfigured = false;
    detail::g_targets.clear();
    detail::g_seen = 0;
    detail::g_injected = 0;
    detail::g_observedBitWidth = 0;
}

bool PBSFaultInjection::IsConfigured() { return detail::g_configured; }
uint64_t PBSFaultInjection::GetInjectedCount() { return detail::g_injected; }
uint64_t PBSFaultInjection::GetSeenCount() { return detail::g_seen; }
uint32_t PBSFaultInjection::GetObservedBitWidth() { return detail::g_observedBitWidth; }
void PBSFaultInjection::SetActive(bool active) { detail::g_active = active; }
bool PBSFaultInjection::IsActive() { return detail::g_active; }
void PBSFaultInjection::SetTransform(PBSFaultTransform transform) { detail::g_currentTransform = transform; }
PBSFaultTransform PBSFaultInjection::GetTransform() { return detail::g_currentTransform; }
bool PBSFaultInjection::IsSiteSelected(PBSFaultSite site) {
    if (!detail::g_configured || detail::g_site != site)
        return false;
    return !detail::g_blockConfigured ||
           (detail::g_selectedBlockKind == detail::g_currentBlockKind &&
            detail::g_selectedEvalAcc == detail::g_currentEvalAcc &&
            detail::g_selectedSubBlock == detail::g_currentSubBlock);
}
bool PBSFaultInjection::IsAnySiteSelected(std::initializer_list<PBSFaultSite> sites) {
    for (const auto site : sites)
        if (IsSiteSelected(site))
            return true;
    return false;
}
bool PBSFaultInjection::IsButterflyTransformSelected(PBSFaultTransform transform) {
    return detail::g_configured && detail::g_site <= PBSFaultSite::BUTTERFLY_STORE &&
           (detail::g_transform == PBSFaultTransform::ANY || detail::g_transform == transform) &&
           (!detail::g_blockConfigured ||
            (detail::g_selectedBlockKind == detail::g_currentBlockKind &&
             detail::g_selectedEvalAcc == detail::g_currentEvalAcc &&
             detail::g_selectedSubBlock == detail::g_currentSubBlock));
}
void PBSFaultInjection::SetCurrentBlock(PBSFaultBlockKind blockKind, uint32_t evalAcc, uint32_t subBlock) {
    detail::g_currentBlockKind = blockKind;
    detail::g_currentEvalAcc = evalAcc;
    detail::g_currentSubBlock = subBlock;
}
void PBSFaultInjection::ClearCurrentBlock() { detail::g_currentBlockKind = PBSFaultBlockKind::NONE; }
uint32_t PBSFaultInjection::GetCurrentEvalAcc() { return detail::g_currentEvalAcc; }
void PBSFaultInjection::BeginPreExplore() {
    detail::g_preExploreBlocks.clear();
    detail::g_preExploreOffsets.fill(0);
    detail::g_preExploring = true;
}
void PBSFaultInjection::EndPreExplore() {
    detail::g_preExploring = false;
}
bool PBSFaultInjection::IsPreExploring() {
    return detail::g_preExploring;
}
void PBSFaultInjection::RegisterPreExploreBlock(PBSFaultSite site, uint64_t candidateCount, uint32_t bitWidth,
                                                 std::array<uint32_t, 4> dimensions) {
    if (!detail::g_preExploring || candidateCount == 0)
        return;
    const auto transform = detail::g_currentTransform;
    const auto slot = static_cast<size_t>(site) * 3 + static_cast<size_t>(transform);
    detail::g_preExploreBlocks.push_back({site, transform, detail::g_preExploreOffsets[slot], candidateCount,
                                          bitWidth, dimensions});
    detail::g_preExploreOffsets[slot] += candidateCount;
}
const std::vector<PBSFaultBlock>& PBSFaultInjection::GetPreExploreBlocks() {
    return detail::g_preExploreBlocks;
}
}  // namespace lbcrypto
