//==================================================================================
// End-to-end protected TFHE/GINX PBS experiment.
//
// The five dominant operators (NTT, INTT, signed decomposition, VCMM, and key
// switching) use the optimized dedicated checks from
// evaluation-operator-protection.cpp. The remaining computation is executed
// twice and compared (DMR); stored data is assumed to be covered by ECC and is
// not duplicated here. This file deliberately keeps the copied PBS
// orchestration outside the production library.
//==================================================================================

// Reuse the exact optimized detector implementations used to produce
// results/operator-protection. Rename that file's command-line entry point so
// this translation unit can provide its own one.
#define main evaluation_operator_protection_main
#include "evaluation-operator-protection.cpp"
#undef main

#include <chrono>
#include <sstream>
#if defined(WITH_PBS_OPERATOR_BREAKDOWN)
#include "pbs-operator-breakdown.h"
#endif

namespace {

using uint128_t = unsigned __int128;

struct CheckedPBSResult {
    LWECiphertext output;
    bool detected = false;
    std::string detector;
};

#if defined(WITH_PBS_OPERATOR_BREAKDOWN)
using BreakdownStage = PBSOperatorBreakdownStage;
using BreakdownData  = PBSOperatorBreakdownData;

void AddBreakdown(BreakdownData& destination, const BreakdownData& source) {
    for (size_t i = 0; i < static_cast<size_t>(BreakdownStage::COUNT); ++i)
        destination.nanoseconds[i] += source.nanoseconds[i];
}

uint64_t BreakdownValue(const BreakdownData& data, BreakdownStage stage) {
    return data.nanoseconds[static_cast<size_t>(stage)];
}

void PrintBreakdownStage(const char* implementation, const char* operation, uint64_t nanoseconds, uint64_t total,
                         uint32_t samples) {
    const double average = static_cast<double>(nanoseconds) / samples;
    const double ratio   = total == 0 ? 0.0 : 100.0 * static_cast<double>(nanoseconds) / total;
    std::cout << "implementation=" << implementation << " operator=" << operation
              << " average_time_ns=" << std::fixed << std::setprecision(3) << average << " ratio_pct=" << ratio
              << '\n';
}

void PrintBreakdown(const char* implementation, const BreakdownData& data, uint32_t samples) {
    const uint64_t total     = BreakdownValue(data, BreakdownStage::TOTAL);
    const uint64_t ntt       = BreakdownValue(data, BreakdownStage::NTT);
    const uint64_t intt      = BreakdownValue(data, BreakdownStage::INTT);
    const uint64_t vcmm      = BreakdownValue(data, BreakdownStage::VCMM);
    const uint64_t decompose = BreakdownValue(data, BreakdownStage::DECOMPOSE);
    const uint64_t keySwitch = BreakdownValue(data, BreakdownStage::KEYSWITCH);
    const uint64_t protectedOperators = ntt + intt + vcmm + decompose + keySwitch;
    const uint64_t other = total > protectedOperators ? total - protectedOperators : 0;
    PrintBreakdownStage(implementation, "NTT", ntt, total, samples);
    PrintBreakdownStage(implementation, "INTT", intt, total, samples);
    PrintBreakdownStage(implementation, "VCMM", vcmm, total, samples);
    PrintBreakdownStage(implementation, "Decompose", decompose, total, samples);
    PrintBreakdownStage(implementation, "KeySwitch", keySwitch, total, samples);
    PrintBreakdownStage(implementation, "Other", other, total, samples);
    PrintBreakdownStage(implementation, "FiveOperators", protectedOperators, total, samples);
    PrintBreakdownStage(implementation, "PBSTotal", total, total, samples);
}

void PrintPBSComparison(double noProtectionTime, double fullProtectionTime) {
    const double overhead = noProtectionTime == 0.0 ? 0.0 : 100.0 * (fullProtectionTime - noProtectionTime) / noProtectionTime;
    std::cout << "comparison=Full-Protection-PBS-vs-No-Protection-PBS operator=PBSTotal overhead_pct="
              << std::fixed << std::setprecision(3) << overhead << '\n';
}
#endif

bool PolysEqual(const NativePoly& lhs, const NativePoly& rhs) {
    return lhs == rhs;
}

bool PolyVectorsEqual(const std::vector<NativePoly>& lhs, const std::vector<NativePoly>& rhs) {
    if (lhs.size() != rhs.size())
        return false;
    for (size_t i = 0; i < lhs.size(); ++i) {
        if (!PolysEqual(lhs[i], rhs[i]))
            return false;
    }
    return true;
}

bool FullCiphertextsEqual(const LWECiphertext& lhs, const LWECiphertext& rhs) {
    if (lhs == nullptr || rhs == nullptr)
        return lhs == rhs;
    return lhs->GetModulus() == rhs->GetModulus() && lhs->GetptModulus() == rhs->GetptModulus() &&
           keyswitch::CiphertextsEqual(lhs, rhs);
}

template <typename T, typename Compute, typename Equal>
T RunDMR(const char* name, Compute&& compute, Equal&& equal, bool& detected, std::string& detector) {
    T primary = compute();
    T replica = compute();
    if (!equal(primary, replica)) {
        detected = true;
        if (detector.empty())
            detector = std::string("DMR:") + name;
    }
    return primary;
}

namespace vcmm {

struct AccumulatingCheckResult {
    bool vcmmDetected = false;
    bool accumulatorDetected = false;
};

// End-to-end PBS variant of the optimized branch-major VCMM check. The four
// branch outputs are consumed while they are hot in the checksum pass, and the
// accumulator modular additions are duplicated in that same pass. Stored
// polynomials remain single-copy because storage faults are covered by ECC.
template <typename CheckAccumulator, bool ReduceCheckBeforeMonomial>
AccumulatingCheckResult VCMMABFTAndAccumulate(const std::vector<NativePoly>& dct,
                                              const ConstRingGSWEvalKey& ek1,
                                              const ConstRingGSWEvalKey& ek2, const NativePoly& monomial,
                                              const NativePoly& monomialNeg, const ChecksumPrecomp& pre,
                                              RLWECiphertext& accumulator) {
    const uint32_t size = dct[0].GetLength();
    const uint64_t modulus = dct[0].GetModulus().ConvertToInt<uint64_t>();
    const NativeInteger nativeModulus(modulus);
    auto output = VCMM(dct, ek1, ek2, monomial, monomialNeg);
    auto& acc = accumulator->GetElements();
    AccumulatingCheckResult result;

    for (uint32_t j = 0; j < size; ++j) {
        CheckAccumulator cin1 = 0;
        CheckAccumulator cin2 = 0;
        for (uint32_t i = 0; i < dct.size(); ++i) {
            const auto digit = static_cast<CheckAccumulator>(dct[i][j].ConvertToInt<uint64_t>());
            cin1 += digit * pre.ev1sum[i][j].ConvertToInt<uint64_t>();
            cin2 += digit * pre.ev2sum[i][j].ConvertToInt<uint64_t>();
        }
        if constexpr (ReduceCheckBeforeMonomial) {
            cin1 %= modulus;
            cin2 %= modulus;
        }
        const auto expected = (cin1 * static_cast<CheckAccumulator>(monomial[j].ConvertToInt<uint64_t>()) +
                               cin2 * static_cast<CheckAccumulator>(monomialNeg[j].ConvertToInt<uint64_t>())) %
                              modulus;
        const uint64_t actual = (output[0][j].ConvertToInt<uint64_t>() + output[1][j].ConvertToInt<uint64_t>() +
                                 output[2][j].ConvertToInt<uint64_t>() + output[3][j].ConvertToInt<uint64_t>()) %
                                modulus;
        result.vcmmDetected |= static_cast<uint64_t>(expected) != actual;

        auto primary0 = acc[0][j];
        auto replica0 = primary0;
        primary0.ModAddFastEq(output[0][j], nativeModulus);
        primary0.ModAddFastEq(output[2][j], nativeModulus);
        replica0.ModAddFastEq(output[0][j], nativeModulus);
        replica0.ModAddFastEq(output[2][j], nativeModulus);

        auto primary1 = acc[1][j];
        auto replica1 = primary1;
        primary1.ModAddFastEq(output[1][j], nativeModulus);
        primary1.ModAddFastEq(output[3][j], nativeModulus);
        replica1.ModAddFastEq(output[1][j], nativeModulus);
        replica1.ModAddFastEq(output[3][j], nativeModulus);

        result.accumulatorDetected |= primary0 != replica0 || primary1 != replica1;
        acc[0][j] = primary0;
        acc[1][j] = primary1;
    }
    return result;
}

}  // namespace vcmm

class ProtectedPBS {
public:
    ProtectedPBS(BinFHEContext& context, BINFHE_PARAMSET level)
        : m_context(context),
          m_level(level),
          m_params(context.GetParams()),
          m_lweParams(m_params->GetLWEParams()),
          m_ringParams(m_params->GetRingGSWParams()),
          m_refreshKey(context.GetRefreshKey()),
          m_switchKey(context.GetSwitchKey()),
          m_keySwitchPre(keyswitch::MakeABFTPrecomp(m_lweParams, m_switchKey)) {
        const auto poly = m_ringParams->GetPolyParams();
        const auto cyclotomicOrder = poly->GetCyclotomicOrder();
        const auto modulus = poly->GetModulus();
        intnat::ChineseRemainderTransformFTTNat<NativeVector>().PreCompute(poly->GetRootOfUnity(), cyclotomicOrder,
                                                                          modulus);
        using CRT = intnat::ChineseRemainderTransformFTTNat<NativeVector>;
        m_roots = CRT::m_rootOfUnityReverseTableByModulus[modulus];
        m_preconRoots = CRT::m_rootOfUnityPreconReverseTableByModulus[modulus];
        m_inverseRoots = CRT::m_rootOfUnityInverseReverseTableByModulus[modulus];
        m_preconInverseRoots = CRT::m_rootOfUnityInversePreconReverseTableByModulus[modulus];
        const uint32_t msb = GetMSB(poly->GetRingDimension() - 1);
        m_nInv = CRT::m_cycloOrderInverseTableByModulus[modulus][msb];
        m_preconNInv = CRT::m_cycloOrderInversePreconTableByModulus[modulus][msb];
        m_transformWeights = ntt::MakeABFTBitReverseWeights(poly->GetRootOfUnity(), cyclotomicOrder, modulus);

        const uint32_t lweDimension = m_lweParams->Getn();
        m_vcmmPre.reserve(lweDimension);
        for (uint32_t i = 0; i < lweDimension; ++i) {
            m_vcmmPre.emplace_back(vcmm::MakeChecksumPrecomp((*m_refreshKey)[0][0][i],
                                                             (*m_refreshKey)[0][1][i], DigitsG2()));
        }
    }

    CheckedPBSResult Bootstrap(ConstLWECiphertext& input) const {
#if defined(WITH_PBS_OPERATOR_BREAKDOWN)
        ScopedPBSOperatorBreakdownTimer totalTimer(BreakdownStage::TOTAL);
#endif
        CheckedPBSResult result;
        if (input == nullptr) {
            result.detected = true;
            result.detector = "invalid-input";
            return result;
        }

        bool detected = false;
        std::string detector;
        const NativeInteger Q = m_lweParams->GetQ();

        // This experiment intentionally supports the normal SMALL_DIM input
        // used by BinFHEContext::Bootstrap.  Its modulus is q, so PBS starts
        // from a ciphertext copy rather than the LARGE_DIM Q-to-q path.
        // Object copies and stored ciphertext data are covered by ECC in the
        // paper's fault model. DMR only duplicates the arithmetic itself.
        auto cct = std::make_shared<LWECiphertextImpl>(*input);
        auto bPrimary = cct->GetB();
        auto bReplica = bPrimary;
        const auto inputConstant = input->GetModulus() >> 2;
        bPrimary.ModAddFastEq(inputConstant, cct->GetModulus());
        bReplica.ModAddFastEq(inputConstant, cct->GetModulus());
        if (bPrimary != bReplica) {
            detected = true;
            detector = "DMR:input-constant-add";
        }
        cct->GetB() = bPrimary;
        if (detected)
            return {nullptr, true, detector};

        auto accumulator = MakeGateAccumulatorDMR(cct, detected, detector);
        if (detected)
            return {nullptr, true, detector};

        if (ForwardNTTChecked(accumulator->GetElements()[1]))
            return {nullptr, true, "NTT"};

        EvalAccChecked(accumulator, cct->GetA(), detected, detector);
        if (detected)
            return {nullptr, true, detector};

        auto& accVec = accumulator->GetElements();
        auto transposeReplica = accVec[0];
        accVec[0] = accVec[0].Transpose();
        transposeReplica = transposeReplica.Transpose();
        if (accVec[0] != transposeReplica)
            return {nullptr, true, "DMR:accumulator-transpose"};

        if (InverseNTTChecked(accVec[0]) || InverseNTTChecked(accVec[1]))
            return {nullptr, true, "INTT"};

        NativeInteger outputBPrimary{Q / (input->GetptModulus() * 2) + 1};
        NativeInteger outputBReplica(outputBPrimary);
        outputBPrimary.ModAddFastEq(accVec[1][0], Q);
        outputBReplica.ModAddFastEq(accVec[1][0], Q);
        if (outputBPrimary != outputBReplica)
            return {nullptr, true, "DMR:output-reconstruction"};
        auto ctExt = std::make_shared<LWECiphertextImpl>(accVec[0].GetValues(), outputBPrimary);
        ctExt->SetptModulus(input->GetptModulus());

        auto output = SwitchCTtoqnChecked(ctExt, detected, detector);
        if (output != nullptr)
            output->SetptModulus(input->GetptModulus());
        return {std::move(output), detected, detector};
    }

private:
    uint32_t DigitsG2() const {
        return (m_ringParams->GetDigitsG() - 1) << 1;
    }

    bool ForwardNTTChecked(NativePoly& poly) const {
#if defined(WITH_PBS_OPERATOR_BREAKDOWN)
        ScopedPBSOperatorBreakdownTimer timer(BreakdownStage::NTT);
#endif
        // NativePoly exposes its coefficients as a const view even for a
        // mutable polynomial. The underlying object is non-const here, so the
        // mutable view is safe and lets the optimized NativeVector kernel run
        // in place without allocation or per-coefficient Poly indirection.
        auto& values = const_cast<NativeVector&>(poly.GetValues());
        const bool detected = ntt::ForwardTransformABFT(m_roots, m_preconRoots, &values, m_transformWeights, m_level);
        poly.OverrideFormat(Format::EVALUATION);
        return detected;
    }

    bool InverseNTTChecked(NativePoly& poly) const {
#if defined(WITH_PBS_OPERATOR_BREAKDOWN)
        ScopedPBSOperatorBreakdownTimer timer(BreakdownStage::INTT);
#endif
        auto& values = const_cast<NativeVector&>(poly.GetValues());
        const bool detected = intt::InverseTransformABFT(m_inverseRoots, m_preconInverseRoots, m_nInv,
                                                         m_preconNInv, &values, m_transformWeights, m_level);
        poly.OverrideFormat(Format::COEFFICIENT);
        return detected;
    }

    LWECiphertext ModSwitchDMR(const NativeInteger& modulus, ConstLWECiphertext& input, bool& detected,
                               std::string& detector, const char* name) const {
        return RunDMR<LWECiphertext>(name, [&]() { return m_context.GetLWEScheme()->ModSwitch(modulus, input); },
                                     FullCiphertextsEqual, detected, detector);
    }

    LWECiphertext SwitchCTtoqnChecked(ConstLWECiphertext& input, bool& detected, std::string& detector) const {
        auto ctMS = ModSwitchDMR(m_lweParams->GetqKS(), input, detected, detector, "pre-keyswitch-modswitch");
        if (detected)
            return nullptr;
        keyswitch::Result checked;
        {
#if defined(WITH_PBS_OPERATOR_BREAKDOWN)
            ScopedPBSOperatorBreakdownTimer timer(BreakdownStage::KEYSWITCH);
#endif
            checked = keyswitch::KeySwitchABFT(m_lweParams, m_switchKey, ctMS, m_keySwitchPre);
        }
        if (checked.detected) {
            detected = true;
            detector = "KeySwitch";
            return nullptr;
        }
        return ModSwitchDMR(m_lweParams->Getq(), checked.output, detected, detector,
                            "post-keyswitch-modswitch");
    }

    RLWECiphertext MakeGateAccumulatorDMR(ConstLWECiphertext& input, bool& detected,
                                          std::string& detector) const {
        return RunDMR<RLWECiphertext>(
            "gate-accumulator-initialization",
            [&]() {
                const NativeInteger q = input->GetModulus();
                const uint32_t qHalf = q.ConvertToInt<uint32_t>() >> 1;
                const auto q1 = m_ringParams->GetGateConst()[static_cast<size_t>(AND)];
                const auto q2 = q1.ModAddFast(NativeInteger(qHalf), q);
                const bool swap = q1 >= q2;
                const auto lb = swap ? q2 : q1;
                const auto ub = swap ? q1 : q2;
                const NativeInteger Q = m_lweParams->GetQ();
                const auto Q2p = Q / (input->GetptModulus() * 2) + 1;
                const auto Q2pNeg = Q - Q2p;
                const auto lv = swap ? Q2p : Q2pNeg;
                const auto uv = swap ? Q2pNeg : Q2p;
                const uint32_t N = m_lweParams->GetN();
                const uint32_t factor = N / qHalf;
                NativeVector m(N, Q);
                NativeInteger b = input->GetB();
                for (uint32_t i = 0; i < N; i += factor) {
                    m[i] = (b >= lb && b < ub) ? lv : uv;
                    b.ModSubFastEq(1, q);
                }
                std::vector<NativePoly> res(2);
                const auto polyParams = m_ringParams->GetPolyParams();
                res[0] = NativePoly(polyParams, Format::EVALUATION, true);
                res[1] = NativePoly(polyParams, Format::COEFFICIENT, false);
                res[1].SetValues(std::move(m), Format::COEFFICIENT);
                return std::make_shared<RLWECiphertextImpl>(std::move(res));
            },
            [](const RLWECiphertext& lhs, const RLWECiphertext& rhs) {
                return lhs != nullptr && rhs != nullptr && PolyVectorsEqual(lhs->GetElements(), rhs->GetElements());
            },
            detected, detector);
    }

    void EvalAccChecked(RLWECiphertext& accumulator, const NativeVector& a, bool& detected,
                        std::string& detector) const {
        const auto modulus = a.GetModulus();
        const auto MbyMod = NativeInteger(2 * m_ringParams->GetN()) / modulus;
        for (uint32_t i = 0; i < a.GetLength() && !detected; ++i) {
            const NativeInteger index = RunDMR<NativeInteger>(
                "accumulator-index",
                [&]() { return NativeInteger(0).ModSubFast(a[i], modulus) * MbyMod; },
                [](const NativeInteger& lhs, const NativeInteger& rhs) { return lhs == rhs; }, detected, detector);
            if (!detected)
                AddToAccChecked((*m_refreshKey)[0][0][i], (*m_refreshKey)[0][1][i], index, m_vcmmPre[i],
                                accumulator, detected, detector);
        }
    }

    void AddToAccChecked(ConstRingGSWEvalKey& ek1, ConstRingGSWEvalKey& ek2, const NativeInteger& index,
                         const vcmm::ChecksumPrecomp& pre, RLWECiphertext& accumulator, bool& detected,
                         std::string& detector) const {
        std::vector<NativePoly> ct(accumulator->GetElements());
        if (InverseNTTChecked(ct[0]) || InverseNTTChecked(ct[1])) {
            detected = true;
            detector = "INTT";
            return;
        }

        auto dct = decompose::MakeZeroOutput(m_ringParams);
        bool decomposeDetected;
        {
#if defined(WITH_PBS_OPERATOR_BREAKDOWN)
            ScopedPBSOperatorBreakdownTimer timer(BreakdownStage::DECOMPOSE);
#endif
            decomposeDetected = decompose::SignedDigitDecomposeOnlineChecked(m_ringParams, ct, dct);
        }
        if (decomposeDetected) {
            detected = true;
            detector = "Decompose";
            return;
        }
        for (auto& digit : dct) {
            if (ForwardNTTChecked(digit)) {
                detected = true;
                detector = "NTT";
                return;
            }
        }

        struct MonomialIndexes {
            uint32_t positive;
            uint32_t negative;
        };
        const auto monomialIndexes = RunDMR<MonomialIndexes>(
            "monomial-selection",
            [&]() {
                const uint32_t MInt = 2 * m_ringParams->GetN();
                const NativeInteger M(MInt);
                const uint32_t positive = index.ConvertToInt<uint32_t>();
                const uint32_t negative = NativeInteger(0).ModSubFast(index, M).ConvertToInt<uint32_t>();
                return MonomialIndexes{positive == MInt ? 0 : positive, negative == MInt ? 0 : negative};
            },
            [](const MonomialIndexes& lhs, const MonomialIndexes& rhs) {
                return lhs.positive == rhs.positive && lhs.negative == rhs.negative;
            },
            detected, detector);
        if (detected)
            return;
        const auto& monomial = m_ringParams->GetMonomial(monomialIndexes.positive);
        const auto& monomialNeg = m_ringParams->GetMonomial(monomialIndexes.negative);

        vcmm::AccumulatingCheckResult vcmmCheck;
        {
#if defined(WITH_PBS_OPERATOR_BREAKDOWN)
            ScopedPBSOperatorBreakdownTimer timer(BreakdownStage::VCMM);
#endif
            vcmmCheck =
                m_level == STD128
                    ? vcmm::VCMMABFTAndAccumulate<uint64_t, true>(dct, ek1, ek2, monomial, monomialNeg, pre,
                                                                  accumulator)
                    : vcmm::VCMMABFTAndAccumulate<uint128_t, false>(dct, ek1, ek2, monomial, monomialNeg, pre,
                                                                    accumulator);
        }
        if (vcmmCheck.vcmmDetected || vcmmCheck.accumulatorDetected) {
            detected = true;
            detector = vcmmCheck.vcmmDetected ? "VCMM" : "DMR:accumulator-update";
        }
    }

    BinFHEContext& m_context;
    BINFHE_PARAMSET m_level;
    std::shared_ptr<BinFHECryptoParams> m_params;
    std::shared_ptr<LWECryptoParams> m_lweParams;
    std::shared_ptr<RingGSWCryptoParams> m_ringParams;
    RingGSWACCKey m_refreshKey;
    LWESwitchingKey m_switchKey;
    keyswitch::ABFTPrecomp m_keySwitchPre;
    NativeVector m_roots;
    NativeVector m_preconRoots;
    NativeVector m_inverseRoots;
    NativeVector m_preconInverseRoots;
    NativeInteger m_nInv;
    NativeInteger m_preconNInv;
    NativeVector m_transformWeights;
    std::vector<vcmm::ChecksumPrecomp> m_vcmmPre;
};

BINFHE_PARAMSET ParseLevel(const std::string& name) {
    if (name == "STD128")
        return STD128;
    if (name == "STD192")
        return STD192;
    if (name == "STD256")
        return STD256;
    throw std::invalid_argument("security must be STD128, STD192, STD256, or all");
}

void RunLevel(const std::string& name, uint32_t samples) {
    const auto level = ParseLevel(name);
    BinFHEContext context;
    context.GenerateBinFHEContext(level, GINX);
    auto secretKey = context.KeyGen();
    context.BTKeyGen(secretKey);
    ProtectedPBS protectedPBS(context, level);

    std::vector<LWECiphertext> inputs;
    inputs.reserve(samples);
    for (uint32_t sample = 0; sample < samples; ++sample)
        inputs.emplace_back(context.Encrypt(secretKey, sample & 1U, SMALL_DIM, 4));

    // Warm both paths without including warm-up work in either breakdown.
#if defined(WITH_PBS_OPERATOR_BREAKDOWN)
    PBSOperatorBreakdown::Enable(false);
#endif
    context.Bootstrap(inputs.front());
    auto protectedWarmup = protectedPBS.Bootstrap(inputs.front());
    if (protectedWarmup.detected)
        throw std::runtime_error("protected PBS warm-up detected an error at " + protectedWarmup.detector);

    uint64_t detectedSamples = 0;
    uint64_t baselineIncorrectSamples = 0;
    uint64_t protectedIncorrectSamples = 0;
    uint64_t outputMismatchSamples = 0;
    std::vector<LWECiphertext> baselineOutputs(samples);
    std::vector<CheckedPBSResult> protectedOutputs(samples);
    std::chrono::nanoseconds baselineElapsed{0};
    std::chrono::nanoseconds protectedElapsed{0};
#if defined(WITH_PBS_OPERATOR_BREAKDOWN)
    BreakdownData baselineBreakdown;
    BreakdownData protectedBreakdown;
#endif

    const auto runBaseline = [&](uint32_t sample) {
#if defined(WITH_PBS_OPERATOR_BREAKDOWN)
        PBSOperatorBreakdown::Reset();
        PBSOperatorBreakdown::Enable(true);
#endif
        const auto start = std::chrono::steady_clock::now();
        baselineOutputs[sample] = context.Bootstrap(inputs[sample]);
        baselineElapsed +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start);
#if defined(WITH_PBS_OPERATOR_BREAKDOWN)
        PBSOperatorBreakdown::Enable(false);
        AddBreakdown(baselineBreakdown, PBSOperatorBreakdown::Snapshot());
#endif
    };
    const auto runProtected = [&](uint32_t sample) {
#if defined(WITH_PBS_OPERATOR_BREAKDOWN)
        PBSOperatorBreakdown::Reset();
        PBSOperatorBreakdown::Enable(true);
#endif
        const auto start = std::chrono::steady_clock::now();
        protectedOutputs[sample] = protectedPBS.Bootstrap(inputs[sample]);
        protectedElapsed +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start);
#if defined(WITH_PBS_OPERATOR_BREAKDOWN)
        PBSOperatorBreakdown::Enable(false);
        AddBreakdown(protectedBreakdown, PBSOperatorBreakdown::Snapshot());
#endif
        detectedSamples += protectedOutputs[sample].detected;
        if (protectedOutputs[sample].detected)
            std::cout << "parameter=" << name << " sample=" << sample << " detected=1 detector="
                      << protectedOutputs[sample].detector << '\n';
    };

    // Alternate execution order to avoid attributing thermal/frequency drift
    // systematically to one implementation.
    for (uint32_t sample = 0; sample < samples; ++sample) {
        if ((sample & 1U) == 0) {
            runBaseline(sample);
            runProtected(sample);
        }
        else {
            runProtected(sample);
            runBaseline(sample);
        }
    }

    for (uint32_t sample = 0; sample < samples; ++sample) {
        const LWEPlaintext expected = sample & 1U;
        LWEPlaintext baselineDecoded = 0;
        context.Decrypt(secretKey, baselineOutputs[sample], &baselineDecoded, 4);
        baselineIncorrectSamples += baselineDecoded != expected;
        if (protectedOutputs[sample].detected)
            continue;
        LWEPlaintext protectedDecoded = 0;
        context.Decrypt(secretKey, protectedOutputs[sample].output, &protectedDecoded, 4);
        protectedIncorrectSamples += protectedDecoded != expected;
        outputMismatchSamples += !FullCiphertextsEqual(baselineOutputs[sample], protectedOutputs[sample].output);
    }

    const double baselineAverage = static_cast<double>(baselineElapsed.count()) / samples;
    const double protectedAverage = static_cast<double>(protectedElapsed.count()) / samples;
    const double overhead = 100.0 * (protectedAverage - baselineAverage) / baselineAverage;
    std::cout << "parameter=" << name << " samples=" << samples << " detected_samples=" << detectedSamples
              << " baseline_incorrect_samples=" << baselineIncorrectSamples
              << " protected_incorrect_undetected_samples=" << protectedIncorrectSamples
              << " output_mismatch_samples=" << outputMismatchSamples
              << " no_protection_average_time_ns=" << baselineAverage
              << " full_protection_average_time_ns=" << protectedAverage << " full_protection_overhead_pct=" << overhead
              << '\n';
#if defined(WITH_PBS_OPERATOR_BREAKDOWN)
    PrintBreakdown("No-Protection-PBS", baselineBreakdown, samples);
    PrintBreakdown("Full-Protection-PBS", protectedBreakdown, samples);
    PrintPBSComparison(baselineAverage, protectedAverage);
#else
    std::cout << "operator_breakdown=disabled rebuild_with=-DWITH_PBS_OPERATOR_BREAKDOWN=ON\n";
#endif
    if (detectedSamples != 0 || baselineIncorrectSamples != 0 || protectedIncorrectSamples != 0 ||
        outputMismatchSamples != 0)
        throw std::runtime_error("fault-free protected PBS validation failed for " + name);
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const std::string security = argc > 1 ? argv[1] : "all";
        const uint32_t samples = argc > 2 ? static_cast<uint32_t>(std::stoul(argv[2])) : 2;
        if (samples == 0)
            throw std::invalid_argument("samples must be nonzero");
        if (security == "all") {
            RunLevel("STD128", samples);
            RunLevel("STD192", samples);
            RunLevel("STD256", samples);
        }
        else {
            RunLevel(security, samples);
        }
        return EXIT_SUCCESS;
    }
    catch (const std::exception& error) {
        std::cerr << "evaluation-pbs-operator-protection: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
