#ifndef OPENFHE_CORE_PBS_SBF_FAULT_INJECTION_H
#define OPENFHE_CORE_PBS_SBF_FAULT_INJECTION_H

#include <algorithm>
#include <array>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <vector>

namespace lbcrypto {

// These names are shared by the end-to-end PBS SBF experiment and its ABFT
// counterpart.  A plan addresses dynamic executions, not static source lines.
enum class PBSFaultSite : uint8_t {
    BUTTERFLY_LOAD, BUTTERFLY_ADD, BUTTERFLY_SUB, BUTTERFLY_MUL, BUTTERFLY_STORE,
    VCMM_INPUT, VCMM_MULTIPLICATION, VCMM_ACCUMULATION, VCMM_OUTPUT,
    DECOMPOSE_DIGIT_EXTRACTION, DECOMPOSE_INTERMEDIATE_DIGIT,
    KEYSWITCH_KEY_LOOKUP, KEYSWITCH_ACCUMULATION, KEYSWITCH_UPDATE,
};

enum class PBSFaultTransform : uint8_t { ANY, NTT, INTT };
enum class PBSFaultBlockKind : uint8_t { NONE, NTT, INTT, DECOMPOSE, VCMM };

struct PBSFaultTarget {
    uint64_t occurrence;
    uint32_t bit;
};

// A compact description of one contiguous dynamic candidate block.  The
// dimensions are operator-specific loop extents and are sufficient to turn a
// sampled rank into a source-level coordinate in the later direct-injection
// implementation.
struct PBSFaultBlock {
    PBSFaultSite site;
    PBSFaultTransform transform;
    uint64_t siteOffset;
    uint64_t candidateCount;
    uint32_t bitWidth;
    std::array<uint32_t, 4> dimensions;
};

namespace detail {
// Defined once in pbs-sbf-fault-injection.cpp.  Match() is intentionally
// header-inline because every coefficient-level injection site calls it.
extern bool g_configured;
extern PBSFaultSite g_site;
extern PBSFaultTransform g_transform;
extern PBSFaultTransform g_currentTransform;
extern std::vector<PBSFaultTarget> g_targets;
extern uint64_t g_seen;
extern uint64_t g_injected;
extern uint32_t g_observedBitWidth;
extern bool g_active;
extern bool g_blockConfigured;
extern PBSFaultBlockKind g_selectedBlockKind;
extern uint32_t g_selectedEvalAcc;
extern uint32_t g_selectedSubBlock;
extern PBSFaultBlockKind g_currentBlockKind;
extern uint32_t g_currentEvalAcc;
extern uint32_t g_currentSubBlock;
extern bool g_preExploring;
extern std::vector<PBSFaultBlock> g_preExploreBlocks;
extern std::array<uint64_t, 42> g_preExploreOffsets;
}  // namespace detail

class PBSFaultInjection {
public:
    // All targets must be distinct occurrences of the same site.  Supplying
    // one or two targets implements the one-/two-fault experiment exactly.
    static void Configure(PBSFaultSite site, std::vector<PBSFaultTarget> targets,
                          PBSFaultTransform transform = PBSFaultTransform::ANY);
    static void ConfigureBlock(PBSFaultSite site, std::vector<PBSFaultTarget> targets,
                               PBSFaultTransform transform, PBSFaultBlockKind blockKind,
                               uint32_t evalAcc, uint32_t subBlock = 0);
    static void Disable();
    static bool IsConfigured();
    static uint64_t GetInjectedCount();
    static uint64_t GetSeenCount();
    static uint32_t GetObservedBitWidth();

    // Returns the chosen bit when this dynamic event is selected.  This must
    // remain inline: calling through the shared library once per polynomial
    // coefficient dominates PBS time for large VCMM candidate spaces.
    static inline bool Match(PBSFaultSite site, uint32_t bitWidth, uint32_t* bit) {
        if (!detail::g_active)
            return false;
        if (!detail::g_configured || detail::g_site != site)
            return false;
        if (detail::g_transform != PBSFaultTransform::ANY &&
            detail::g_transform != detail::g_currentTransform)
            return false;
        const uint64_t occurrence = detail::g_seen++;
        detail::g_observedBitWidth = bitWidth;
        for (const auto& target : detail::g_targets) {
            if (target.occurrence == occurrence) {
                if (target.bit >= bitWidth)
                    return false;
                *bit = target.bit;
                ++detail::g_injected;
                return true;
            }
        }
        return false;
    }
    static void SetActive(bool active);
    static bool IsActive();
    static void SetTransform(PBSFaultTransform transform);
    static PBSFaultTransform GetTransform();
    static bool IsSiteSelected(PBSFaultSite site);
    static bool IsAnySiteSelected(std::initializer_list<PBSFaultSite> sites);
    static bool IsButterflyTransformSelected(PBSFaultTransform transform);
    static void SetCurrentBlock(PBSFaultBlockKind blockKind, uint32_t evalAcc, uint32_t subBlock = 0);
    static void ClearCurrentBlock();
    static uint32_t GetCurrentEvalAcc();
    // Pre-exploration registers loop blocks, never individual coefficients.
    // This API is compiled only with the SBF experiment feature.
    static void BeginPreExplore();
    static void EndPreExplore();
    static bool IsPreExploring();
    static void RegisterPreExploreBlock(PBSFaultSite site, uint64_t candidateCount, uint32_t bitWidth,
                                        std::array<uint32_t, 4> dimensions = {});
    static const std::vector<PBSFaultBlock>& GetPreExploreBlocks();

    template <typename Integer>
    static void Inject(Integer& value, PBSFaultSite site,
                       uint32_t bitWidth = std::numeric_limits<typename Integer::Integer>::digits) {
        uint32_t bit = 0;
        if (!Match(site, bitWidth, &bit))
            return;
        using Raw = typename Integer::Integer;
        if (bit >= std::numeric_limits<Raw>::digits)
            return;
        value = Integer(value.template ConvertToInt<Raw>() ^ (static_cast<Raw>(1) << bit));
    }

    template <typename Raw>
    static void InjectRaw(Raw& value, PBSFaultSite site,
                          uint32_t bitWidth = std::numeric_limits<Raw>::digits) {
        uint32_t bit = 0;
        if (!Match(site, bitWidth, &bit) || bit >= std::numeric_limits<Raw>::digits)
            return;
        value ^= (static_cast<Raw>(1) << bit);
    }

};

class ScopedPBSFaultTransform {
public:
    explicit ScopedPBSFaultTransform(PBSFaultTransform transform) : m_previous(PBSFaultInjection::GetTransform()) {
        PBSFaultInjection::SetTransform(transform);
    }
    ~ScopedPBSFaultTransform() { PBSFaultInjection::SetTransform(m_previous); }
private:
    PBSFaultTransform m_previous;
};

class ScopedPBSFaultInjection {
public:
    explicit ScopedPBSFaultInjection(PBSFaultSite site) : m_previous(PBSFaultInjection::IsActive()) {
        PBSFaultInjection::SetActive(m_previous || PBSFaultInjection::IsSiteSelected(site));
    }
    ScopedPBSFaultInjection(std::initializer_list<PBSFaultSite> sites) : m_previous(PBSFaultInjection::IsActive()) {
        PBSFaultInjection::SetActive(m_previous || PBSFaultInjection::IsAnySiteSelected(sites));
    }
    explicit ScopedPBSFaultInjection(PBSFaultTransform transform) : m_previous(PBSFaultInjection::IsActive()) {
        PBSFaultInjection::SetActive(m_previous || PBSFaultInjection::IsButterflyTransformSelected(transform));
    }
    ~ScopedPBSFaultInjection() {
        PBSFaultInjection::SetActive(m_previous);
    }
    ScopedPBSFaultInjection(const ScopedPBSFaultInjection&) = delete;
    ScopedPBSFaultInjection& operator=(const ScopedPBSFaultInjection&) = delete;

private:
    bool m_previous;
};

}  // namespace lbcrypto

#endif
