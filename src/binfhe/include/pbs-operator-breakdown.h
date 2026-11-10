#ifndef OPENFHE_BINFHE_PBS_OPERATOR_BREAKDOWN_H
#define OPENFHE_BINFHE_PBS_OPERATOR_BREAKDOWN_H

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>

namespace lbcrypto {

enum class PBSOperatorBreakdownStage : std::size_t {
    TOTAL = 0,
    NTT,
    INTT,
    DECOMPOSE,
    VCMM,
    KEYSWITCH,
    COUNT,
};

struct PBSOperatorBreakdownData {
    std::array<uint64_t, static_cast<std::size_t>(PBSOperatorBreakdownStage::COUNT)> nanoseconds{};
};

class PBSOperatorBreakdown {
public:
    static void Enable(bool enabled);
    static bool IsEnabled();
    static void Reset();
    static void Add(PBSOperatorBreakdownStage stage, uint64_t nanoseconds);
    static PBSOperatorBreakdownData Snapshot();
};

class ScopedPBSOperatorBreakdownTimer {
public:
    explicit ScopedPBSOperatorBreakdownTimer(PBSOperatorBreakdownStage stage)
        : m_stage(stage), m_enabled(PBSOperatorBreakdown::IsEnabled()) {
        if (m_enabled)
            m_start = std::chrono::steady_clock::now();
    }

    ~ScopedPBSOperatorBreakdownTimer() {
        if (m_enabled) {
            const auto elapsed = std::chrono::steady_clock::now() - m_start;
            PBSOperatorBreakdown::Add(
                m_stage, static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count()));
        }
    }

private:
    PBSOperatorBreakdownStage m_stage;
    bool m_enabled;
    std::chrono::steady_clock::time_point m_start{};
};

}  // namespace lbcrypto

#endif
