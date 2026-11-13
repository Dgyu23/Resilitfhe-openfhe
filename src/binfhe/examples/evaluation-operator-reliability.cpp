//==================================================================================
// Reliability evaluation of the PBS operators protected by our ABFT checks.
//
// The copied operator bodies deliberately remain local to this example.  This
// lets faults be applied to exactly the evaluated body, while checksum state
// and reference computations remain trusted.
//==================================================================================

#include "binfhecontext.h"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

using namespace lbcrypto;

namespace {

// Keep these names and the fault model aligned with
// evaluation-pbs-operator-fault-injection.cpp.  A plan selects one site and
// one or two distinct dynamic occurrences at that site; each occurrence gets
// a random single-bit flip.  A detector still produces one boolean result per
// sample, independently of whether the sample contains one or two faults.
enum class FaultSite {
    ButterflyLoad,
    ButterflyAdd,
    ButterflySub,
    ButterflyMul,
    ButterflyStore,
    VCMMInput,
    VCMMMultiplication,
    VCMMAccumulation,
    VCMMOutput,
    DecomposeDigitExtraction,
    DecomposeIntermediateDigit,
    KeySwitchKeyLookup,
    KeySwitchAccumulation,
    KeySwitchUpdate,
};

struct FaultPlan {
    FaultSite site;
    std::set<uint64_t> occurrences;
    uint64_t seen = 0;
    uint32_t applied = 0;
    std::mt19937_64* prng = nullptr;

    bool Select(FaultSite candidate) {
        if (candidate != site)
            return false;
        return occurrences.count(seen++) != 0;
    }

    template <typename Integer>
    void Inject(Integer& value, FaultSite candidate, uint32_t bitWidth) {
        if (!Select(candidate))
            return;
        if (bitWidth == 0)
            throw std::runtime_error("fault site has no valid bits");
        std::uniform_int_distribution<uint32_t> bit(0, bitWidth - 1);
        using Raw = typename Integer::Integer;
        value = Integer(value.template ConvertToInt<Raw>() ^ (static_cast<Raw>(1) << bit(*prng)));
        ++applied;
    }

    template <typename Raw>
    void InjectRaw(Raw& value, FaultSite candidate, uint32_t bitWidth) {
        if (!Select(candidate))
            return;
        if (bitWidth == 0)
            throw std::runtime_error("fault site has no valid bits");
        std::uniform_int_distribution<uint32_t> bit(0, bitWidth - 1);
        value ^= (static_cast<Raw>(1) << bit(*prng));
        ++applied;
    }
};

FaultPlan* gFaultPlan = nullptr;

template <typename Integer>
inline void Inject(Integer& value, FaultSite site, uint32_t bitWidth) {
    if (gFaultPlan != nullptr)
        gFaultPlan->Inject(value, site, bitWidth);
}

template <typename Raw>
inline void InjectRaw(Raw& value, FaultSite site, uint32_t bitWidth) {
    if (gFaultPlan != nullptr)
        gFaultPlan->InjectRaw(value, site, bitWidth);
}

struct ReliabilityStats {
    uint64_t detectedSamples = 0;
    uint64_t maskedSamples = 0;
};

FaultPlan MakeFaultPlan(FaultSite site, uint64_t candidates, uint32_t faultsPerSample,
                        std::mt19937_64& prng) {
    if (faultsPerSample == 0 || faultsPerSample > 2 || candidates < faultsPerSample)
        throw std::invalid_argument("faults_per_sample must be 1 or 2 and fit the selected site");
    std::uniform_int_distribution<uint64_t> occurrence(0, candidates - 1);
    std::set<uint64_t> selected;
    while (selected.size() < faultsPerSample)
        selected.insert(occurrence(prng));
    return {site, std::move(selected), 0, 0, &prng};
}

const char* FaultSiteName(FaultSite site) {
    switch (site) {
        case FaultSite::ButterflyLoad: return "butterfly-load";
        case FaultSite::ButterflyAdd: return "butterfly-add";
        case FaultSite::ButterflySub: return "butterfly-sub";
        case FaultSite::ButterflyMul: return "butterfly-mul";
        case FaultSite::ButterflyStore: return "butterfly-store";
        case FaultSite::VCMMInput: return "input";
        case FaultSite::VCMMMultiplication: return "multiplication";
        case FaultSite::VCMMAccumulation: return "accumulation";
        case FaultSite::VCMMOutput: return "output";
        case FaultSite::DecomposeDigitExtraction: return "digit-extraction";
        case FaultSite::DecomposeIntermediateDigit: return "intermediate-digit";
        case FaultSite::KeySwitchKeyLookup: return "key-lookup";
        case FaultSite::KeySwitchAccumulation: return "accumulation";
        case FaultSite::KeySwitchUpdate: return "update";
    }
    return "unknown";
}

void PrintReliability(const std::string& operation, const std::string& parameterName, FaultSite site,
                      uint32_t faultsPerSample, uint64_t samples, const ReliabilityStats& stats) {
    const uint64_t undetected = samples - stats.detectedSamples;
    std::cout << std::fixed << std::setprecision(6)
              << "operator=" << operation << " parameter=" << parameterName
              << " site=" << FaultSiteName(site) << " faults_per_sample=" << faultsPerSample
              << " samples=" << samples << " detected_samples=" << stats.detectedSamples
              << " undetected_samples=" << undetected << " masked_samples=" << stats.maskedSamples
              << " detection_rate_pct=" << 100.0 * static_cast<double>(stats.detectedSamples) / samples
              << " masking_rate_pct=" << 100.0 * static_cast<double>(stats.maskedSamples) / samples << '\n';
}

bool VectorsEqual(const NativeVector& lhs, const NativeVector& rhs) {
    if (lhs.GetLength() != rhs.GetLength() || lhs.GetModulus() != rhs.GetModulus())
        return false;

    for (uint32_t i = 0; i < lhs.GetLength(); ++i) {
        if (lhs[i] != rhs[i])
            return false;
    }
    return true;
}

namespace ntt {

// Unprotected baseline copied from OpenFHE's
// NumberTheoreticTransformNat<NativeVector>::ForwardTransformToBitReverseInPlace.
// It is retained here so the benchmark evaluates the operator body directly.
void ForwardTransformToBitReverseInPlace(const NativeVector& rootOfUnityTable,
                                         const NativeVector& preconRootOfUnityTable, NativeVector* element) {
    const auto modulus{element->GetModulus()};
    const uint32_t n(element->GetLength() >> 1);
    for (uint32_t m{1}, t{n}, logt{GetMSB(t)}; m < n; m <<= 1, t >>= 1, --logt) {
        for (uint32_t i{0}; i < m; ++i) {
            auto omega{rootOfUnityTable[i + m]};
            auto preconOmega{preconRootOfUnityTable[i + m]};
            for (uint32_t j1{i << logt}, j2{j1 + t}; j1 < j2; ++j1) {
                auto omegaFactor{(*element)[j1 + t]};
                omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
                auto loVal{(*element)[j1 + 0]};
#if defined(__GNUC__) && !defined(__clang__)
                auto hiVal{loVal + omegaFactor};
                if (hiVal >= modulus)
                    hiVal -= modulus;
                if (loVal < omegaFactor)
                    loVal += modulus;
                loVal -= omegaFactor;
                (*element)[j1 + 0] = hiVal;
                (*element)[j1 + t] = loVal;
#else
                (*element)[j1 + 0] += omegaFactor - (omegaFactor >= (modulus - loVal) ? modulus : 0);
                if (omegaFactor > loVal)
                    loVal += modulus;
                (*element)[j1 + t] = loVal - omegaFactor;
#endif
            }
        }
    }

    // The PBS NTT peels off the final stage for performance.
    for (uint32_t i{0}; i < (n << 1); i += 2) {
        auto omegaFactor{(*element)[i + 1]};
        auto omega{rootOfUnityTable[(i >> 1) + n]};
        auto preconOmega{preconRootOfUnityTable[(i >> 1) + n]};
        omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
        auto loVal{(*element)[i + 0]};
#if defined(__GNUC__) && !defined(__clang__)
        auto hiVal{loVal + omegaFactor};
        if (hiVal >= modulus)
            hiVal -= modulus;
        if (loVal < omegaFactor)
            loVal += modulus;
        loVal -= omegaFactor;
        (*element)[i + 0] = hiVal;
        (*element)[i + 1] = loVal;
#else
        (*element)[i + 0] += omegaFactor - (omegaFactor >= (modulus - loVal) ? modulus : 0);
        if (omegaFactor > loVal)
            loVal += modulus;
        (*element)[i + 1] = loVal - omegaFactor;
#endif
    }
}

// Generate the bit-reversed checksum weights used by the fused ABFT NTT.
NativeVector MakeABFTBitReverseWeights(const NativeInteger& rootOfUnity, uint32_t cyclotomicOrder,
                                       const NativeInteger& modulus) {
    const uint32_t length   = cyclotomicOrder >> 1;
    const uint32_t bitWidth = GetMSB(length - 1);
    const auto nInverse     = NativeInteger(length).ModInverse(modulus);
    const auto twoNInverse  = NativeInteger(2).ModMul(nInverse, modulus);
    const auto rootInverse  = rootOfUnity.ModInverse(modulus);
    const auto step         = rootInverse.ModMul(rootInverse, modulus);
    const auto mu           = modulus.ComputeMu();

    NativeVector weights(length, modulus);
    auto rootPower = rootInverse;
    for (uint32_t index = 0; index < length; ++index) {
        const auto denominator = NativeInteger(1).ModSub(rootPower, modulus);
        const auto weight      = twoNInverse.ModMul(denominator.ModInverse(modulus), modulus);
        weights[ReverseBits(index, bitWidth)] = weight;
        rootPower.ModMulFastEq(step, modulus, mu);
    }
    return weights;
}

// The last NTT stage and output decoding are fused.  Accumulator is selected at
// the call site, so the generated STD128 and STD192/STD256 loops retain their
// respective uint64_t and __uint128_t arithmetic without a per-iteration branch.
template <typename Accumulator>
inline uint64_t FusedLastStageAndOutputChecksum(const NativeVector& rootOfUnityTable,
                                                const NativeVector& preconRootOfUnityTable, NativeVector* element,
                                                const NativeVector& aBitReverse, uint32_t n, uint64_t nativeQ) {
    const auto modulus    = element->GetModulus();
    const uint32_t length = n << 1;
    Accumulator cOutAcc   = 0;
    for (uint32_t i{0}; i < length; i += 2) {
        auto omegaFactor = (*element)[i + 1];
        Inject(omegaFactor, FaultSite::ButterflyLoad, modulus.GetMSB());
        omegaFactor.ModMulFastConstEq(rootOfUnityTable[(i >> 1) + n], modulus,
                                       preconRootOfUnityTable[(i >> 1) + n]);
        Inject(omegaFactor, FaultSite::ButterflyMul, modulus.GetMSB());
        auto loVal = (*element)[i];
        Inject(loVal, FaultSite::ButterflyLoad, modulus.GetMSB());
        auto hiVal = loVal + omegaFactor;
        if (hiVal >= modulus)
            hiVal -= modulus;
        Inject(hiVal, FaultSite::ButterflyAdd, modulus.GetMSB());
        if (loVal < omegaFactor)
            loVal += modulus;
        loVal -= omegaFactor;
        Inject(loVal, FaultSite::ButterflySub, modulus.GetMSB());
        Inject(hiVal, FaultSite::ButterflyStore, modulus.GetMSB());
        Inject(loVal, FaultSite::ButterflyStore, modulus.GetMSB());

        cOutAcc += static_cast<Accumulator>(hiVal.ConvertToInt()) *
                   static_cast<uint64_t>(aBitReverse[i].ConvertToInt());
        cOutAcc += static_cast<Accumulator>(loVal.ConvertToInt()) *
                   static_cast<uint64_t>(aBitReverse[i + 1].ConvertToInt());
        (*element)[i]     = hiVal;
        (*element)[i + 1] = loVal;
    }
    return static_cast<uint64_t>(cOutAcc % static_cast<Accumulator>(nativeQ));
}

// Performs the same forward NTT as ForwardTransformToBitReverseInPlace(), while
// fusing the ABFT output checksum into the final butterfly stage.  A true return
// value denotes a checksum mismatch.
bool ForwardTransformABFT(const NativeVector& rootOfUnityTable,
                                   const NativeVector& preconRootOfUnityTable, NativeVector* element,
                                   const NativeVector& aBitReverse, BINFHE_PARAMSET level) {
    const auto modulus    = element->GetModulus();
    const uint32_t n      = element->GetLength() >> 1;
    const uint64_t nativeQ = static_cast<uint64_t>(modulus.ConvertToInt());

    uint64_t cInAcc = 0;
    {
        const auto omega       = rootOfUnityTable[1];
        const auto preconOmega = preconRootOfUnityTable[1];
        for (uint32_t j1 = 0; j1 < n; ++j1) {
            const auto hiIn = (*element)[j1];
            const auto loIn = (*element)[j1 + n];
            cInAcc += static_cast<uint64_t>(hiIn.ConvertToInt());
            cInAcc += static_cast<uint64_t>(loIn.ConvertToInt());

            auto omegaFactor = loIn;
            Inject(omegaFactor, FaultSite::ButterflyLoad, modulus.GetMSB());
            omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
            Inject(omegaFactor, FaultSite::ButterflyMul, modulus.GetMSB());
            auto loVal = hiIn;
            Inject(loVal, FaultSite::ButterflyLoad, modulus.GetMSB());
            auto hiVal = loVal + omegaFactor;
            if (hiVal >= modulus)
                hiVal -= modulus;
            Inject(hiVal, FaultSite::ButterflyAdd, modulus.GetMSB());
            if (loVal < omegaFactor)
                loVal += modulus;
            loVal -= omegaFactor;
            Inject(loVal, FaultSite::ButterflySub, modulus.GetMSB());
            Inject(hiVal, FaultSite::ButterflyStore, modulus.GetMSB());
            Inject(loVal, FaultSite::ButterflyStore, modulus.GetMSB());
            (*element)[j1]     = hiVal;
            (*element)[j1 + n] = loVal;
        }
    }

    for (uint32_t m{2}, t{n >> 1}, logt{GetMSB(t)}; m < n; m <<= 1, t >>= 1, --logt) {
        for (uint32_t i{0}; i < m; ++i) {
            const auto omega       = rootOfUnityTable[i + m];
            const auto preconOmega = preconRootOfUnityTable[i + m];
            for (uint32_t j1{i << logt}, j2{j1 + t}; j1 < j2; ++j1) {
                auto omegaFactor = (*element)[j1 + t];
                Inject(omegaFactor, FaultSite::ButterflyLoad, modulus.GetMSB());
                omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
                Inject(omegaFactor, FaultSite::ButterflyMul, modulus.GetMSB());
                auto loVal = (*element)[j1];
                Inject(loVal, FaultSite::ButterflyLoad, modulus.GetMSB());
                auto hiVal = loVal + omegaFactor;
                if (hiVal >= modulus)
                    hiVal -= modulus;
                Inject(hiVal, FaultSite::ButterflyAdd, modulus.GetMSB());
                if (loVal < omegaFactor)
                    loVal += modulus;
                loVal -= omegaFactor;
                Inject(loVal, FaultSite::ButterflySub, modulus.GetMSB());
                Inject(hiVal, FaultSite::ButterflyStore, modulus.GetMSB());
                Inject(loVal, FaultSite::ButterflyStore, modulus.GetMSB());
                (*element)[j1]     = hiVal;
                (*element)[j1 + t] = loVal;
            }
        }
    }

    const uint64_t cOut = level == STD128
                              ? FusedLastStageAndOutputChecksum<uint64_t>(rootOfUnityTable, preconRootOfUnityTable,
                                                                          element, aBitReverse, n, nativeQ)
                              : FusedLastStageAndOutputChecksum<__uint128_t>(rootOfUnityTable,
                                                                             preconRootOfUnityTable, element,
                                                                             aBitReverse, n, nativeQ);

    return cInAcc % nativeQ != cOut;
}

}  // namespace ntt

namespace intt {

// Unprotected baseline copied from OpenFHE's GS inverse transform.  Input is
// bit-reversed and output is in standard coefficient order.
void InverseTransformFromBitReverseInPlace(const NativeVector& inverseRoots, const NativeVector& preconInverseRoots,
                                           const NativeInteger& nInv, const NativeInteger& preconNInv,
                                           NativeVector* element) {
    const auto modulus{element->GetModulus()};
    const uint32_t n(element->GetLength());
    const auto omega1Inv{inverseRoots[1].ModMulFastConst(nInv, modulus, preconNInv)};
    const auto preconOmega1Inv{omega1Inv.PrepModMulConst(modulus)};

    if (n > 2) {
        for (uint32_t i{0}; i < n; i += 2) {
            const auto omega{inverseRoots[(i + n) >> 1]};
            const auto preconOmega{preconInverseRoots[(i + n) >> 1]};
            auto loVal{(*element)[i]};
            auto hiVal{(*element)[i + 1]};
#if defined(__GNUC__) && !defined(__clang__)
            auto omegaFactor{loVal};
            if (omegaFactor < hiVal)
                omegaFactor += modulus;
            omegaFactor -= hiVal;
            loVal += hiVal;
            if (loVal >= modulus)
                loVal -= modulus;
            omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
            (*element)[i]     = loVal;
            (*element)[i + 1] = omegaFactor;
#else
            auto omegaFactor{loVal + (hiVal > loVal ? modulus : 0) - hiVal};
            loVal += hiVal - (hiVal >= (modulus - loVal) ? modulus : 0);
            (*element)[i] = loVal;
            omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
            (*element)[i + 1] = omegaFactor;
#endif
        }
    }

    for (uint32_t m{n >> 2}, t{2}, logt{2}; m > 1; m >>= 1, t <<= 1, ++logt) {
        for (uint32_t i{0}; i < m; ++i) {
            const auto omega{inverseRoots[i + m]};
            const auto preconOmega{preconInverseRoots[i + m]};
            for (uint32_t j1{i << logt}, j2{j1 + t}; j1 < j2; ++j1) {
                auto loVal{(*element)[j1]};
                auto hiVal{(*element)[j1 + t]};
#if defined(__GNUC__) && !defined(__clang__)
                auto omegaFactor{loVal};
                if (omegaFactor < hiVal)
                    omegaFactor += modulus;
                omegaFactor -= hiVal;
                loVal += hiVal;
                if (loVal >= modulus)
                    loVal -= modulus;
                omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
                (*element)[j1]     = loVal;
                (*element)[j1 + t] = omegaFactor;
#else
                (*element)[j1] += hiVal - (hiVal >= (modulus - loVal) ? modulus : 0);
                auto omegaFactor = loVal + (hiVal > loVal ? modulus : 0) - hiVal;
                omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
                (*element)[j1 + t] = omegaFactor;
#endif
            }
        }
    }

    const uint32_t half{n >> 1};
    for (uint32_t j1{0}; j1 < half; ++j1) {
        auto loVal{(*element)[j1]};
        auto hiVal{(*element)[j1 + half]};
#if defined(__GNUC__) && !defined(__clang__)
        auto omegaFactor{loVal};
        if (omegaFactor < hiVal)
            omegaFactor += modulus;
        omegaFactor -= hiVal;
        loVal += hiVal;
        if (loVal >= modulus)
            loVal -= modulus;
        omegaFactor.ModMulFastConstEq(omega1Inv, modulus, preconOmega1Inv);
        (*element)[j1]        = loVal;
        (*element)[j1 + half] = omegaFactor;
#else
        (*element)[j1] += hiVal - (hiVal >= (modulus - loVal) ? modulus : 0);
        auto omegaFactor = loVal + (hiVal > loVal ? modulus : 0) - hiVal;
        omegaFactor.ModMulFastConstEq(omega1Inv, modulus, preconOmega1Inv);
        (*element)[j1 + half] = omegaFactor;
#endif
    }
    for (uint32_t i{0}; i < half; ++i)
        (*element)[i].ModMulFastConstEq(nInv, modulus, preconNInv);
}

template <typename Accumulator>
inline uint64_t FusedFirstStageAndInputChecksum(const NativeVector& inverseRoots,
                                                const NativeVector& preconInverseRoots, NativeVector* element,
                                                const NativeVector& aBitReverse, uint64_t nativeQ) {
    const auto modulus = element->GetModulus();
    const uint32_t n = element->GetLength();
    Accumulator cInAcc = 0;
    if (n == 2) {
        cInAcc += static_cast<Accumulator>((*element)[0].ConvertToInt()) *
                  static_cast<uint64_t>(aBitReverse[0].ConvertToInt());
        cInAcc += static_cast<Accumulator>((*element)[1].ConvertToInt()) *
                  static_cast<uint64_t>(aBitReverse[1].ConvertToInt());
    }
    else {
        for (uint32_t i{0}; i < n; i += 2) {
            cInAcc += static_cast<Accumulator>((*element)[i].ConvertToInt()) *
                      static_cast<uint64_t>(aBitReverse[i].ConvertToInt());
            cInAcc += static_cast<Accumulator>((*element)[i + 1].ConvertToInt()) *
                      static_cast<uint64_t>(aBitReverse[i + 1].ConvertToInt());
            const auto omega = inverseRoots[(i + n) >> 1];
            const auto preconOmega = preconInverseRoots[(i + n) >> 1];
            auto loVal = (*element)[i];
            auto hiVal = (*element)[i + 1];
            Inject(loVal, FaultSite::ButterflyLoad, modulus.GetMSB());
            Inject(hiVal, FaultSite::ButterflyLoad, modulus.GetMSB());
            auto omegaFactor = loVal;
            if (omegaFactor < hiVal)
                omegaFactor += modulus;
            omegaFactor -= hiVal;
            Inject(omegaFactor, FaultSite::ButterflySub, modulus.GetMSB());
            loVal += hiVal;
            if (loVal >= modulus)
                loVal -= modulus;
            Inject(loVal, FaultSite::ButterflyAdd, modulus.GetMSB());
            omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
            Inject(omegaFactor, FaultSite::ButterflyMul, modulus.GetMSB());
            Inject(loVal, FaultSite::ButterflyStore, modulus.GetMSB());
            Inject(omegaFactor, FaultSite::ButterflyStore, modulus.GetMSB());
            (*element)[i]     = loVal;
            (*element)[i + 1] = omegaFactor;
        }
    }
    return static_cast<uint64_t>(cInAcc % static_cast<Accumulator>(nativeQ));
}

// Fused GS-INTT dual: a_br weighted input encoding is fused into the peeled
// first stage; the all-one output decoding is split across its two final
// normalization stages.
bool InverseTransformABFT(const NativeVector& inverseRoots, const NativeVector& preconInverseRoots,
                          const NativeInteger& nInv, const NativeInteger& preconNInv, NativeVector* element,
                          const NativeVector& aBitReverse, BINFHE_PARAMSET level) {
    const auto modulus = element->GetModulus();
    const uint32_t n = element->GetLength();
    const uint64_t nativeQ = static_cast<uint64_t>(modulus.ConvertToInt());
    const uint64_t cIn = level == STD128
                             ? FusedFirstStageAndInputChecksum<uint64_t>(inverseRoots, preconInverseRoots, element,
                                                                         aBitReverse, nativeQ)
                             : FusedFirstStageAndInputChecksum<__uint128_t>(inverseRoots, preconInverseRoots,
                                                                            element, aBitReverse, nativeQ);

    for (uint32_t m{n >> 2}, t{2}, logt{2}; m > 1; m >>= 1, t <<= 1, ++logt) {
        for (uint32_t i{0}; i < m; ++i) {
            const auto omega = inverseRoots[i + m];
            const auto preconOmega = preconInverseRoots[i + m];
            for (uint32_t j1{i << logt}, j2{j1 + t}; j1 < j2; ++j1) {
                auto loVal = (*element)[j1];
                auto hiVal = (*element)[j1 + t];
                Inject(loVal, FaultSite::ButterflyLoad, modulus.GetMSB());
                Inject(hiVal, FaultSite::ButterflyLoad, modulus.GetMSB());
                auto omegaFactor = loVal;
                if (omegaFactor < hiVal)
                    omegaFactor += modulus;
                omegaFactor -= hiVal;
                Inject(omegaFactor, FaultSite::ButterflySub, modulus.GetMSB());
                loVal += hiVal;
                if (loVal >= modulus)
                    loVal -= modulus;
                Inject(loVal, FaultSite::ButterflyAdd, modulus.GetMSB());
                omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
                Inject(omegaFactor, FaultSite::ButterflyMul, modulus.GetMSB());
                Inject(loVal, FaultSite::ButterflyStore, modulus.GetMSB());
                Inject(omegaFactor, FaultSite::ButterflyStore, modulus.GetMSB());
                (*element)[j1]     = loVal;
                (*element)[j1 + t] = omegaFactor;
            }
        }
    }

    const auto omega1Inv = inverseRoots[1].ModMulFastConst(nInv, modulus, preconNInv);
    const auto preconOmega1Inv = omega1Inv.PrepModMulConst(modulus);
    const uint32_t half = n >> 1;
    uint64_t cOutAcc = 0;
    for (uint32_t j1{0}; j1 < half; ++j1) {
        auto loVal = (*element)[j1];
        auto hiVal = (*element)[j1 + half];
        Inject(loVal, FaultSite::ButterflyLoad, modulus.GetMSB());
        Inject(hiVal, FaultSite::ButterflyLoad, modulus.GetMSB());
        auto omegaFactor = loVal;
        if (omegaFactor < hiVal)
            omegaFactor += modulus;
        omegaFactor -= hiVal;
        Inject(omegaFactor, FaultSite::ButterflySub, modulus.GetMSB());
        loVal += hiVal;
        if (loVal >= modulus)
            loVal -= modulus;
        Inject(loVal, FaultSite::ButterflyAdd, modulus.GetMSB());
        omegaFactor.ModMulFastConstEq(omega1Inv, modulus, preconOmega1Inv);
        Inject(omegaFactor, FaultSite::ButterflyMul, modulus.GetMSB());
        Inject(loVal, FaultSite::ButterflyStore, modulus.GetMSB());
        Inject(omegaFactor, FaultSite::ButterflyStore, modulus.GetMSB());
        cOutAcc += static_cast<uint64_t>(omegaFactor.ConvertToInt());
        (*element)[j1]        = loVal;
        (*element)[j1 + half] = omegaFactor;
    }
    for (uint32_t i{0}; i < half; ++i) {
        auto value = (*element)[i];
        value.ModMulFastConstEq(nInv, modulus, preconNInv);
        cOutAcc += static_cast<uint64_t>(value.ConvertToInt());
        (*element)[i] = value;
    }
    return cIn != cOutAcc % nativeQ;
}

}  // namespace intt

namespace decompose {

using Input  = std::vector<NativePoly>;
using Output = std::vector<NativePoly>;

Output MakeZeroOutput(const std::shared_ptr<RingGSWCryptoParams>& params) {
    const uint32_t digitsG2 = (params->GetDigitsG() - 1) << 1;
    return Output(digitsG2, NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));
}

Input MakeRandomInput(const std::shared_ptr<RingGSWCryptoParams>& params, std::mt19937_64& prng) {
    Input input(2, NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));
    const uint64_t q = params->GetQ().ConvertToInt<uint64_t>();
    std::uniform_int_distribution<uint64_t> distribution(0, q - 1);
    for (auto& component : input) {
        for (uint32_t k = 0; k < params->GetN(); ++k)
            component[k] = NativeInteger(distribution(prng));
    }
    return input;
}

bool OutputsEqual(const Output& lhs, const Output& rhs) {
    if (lhs.size() != rhs.size())
        return false;
    for (size_t d = 0; d < lhs.size(); ++d) {
        if (lhs[d].GetLength() != rhs[d].GetLength())
            return false;
        for (uint32_t k = 0; k < lhs[d].GetLength(); ++k) {
            if (lhs[d][k] != rhs[d][k])
                return false;
        }
    }
    return true;
}

// Unprotected baseline copied from OpenFHE's
// RingGSWAccumulator::SignedDigitDecompose.  The output is intentionally
// accumulated with +=, so callers must provide a zeroed output buffer, exactly
// as the PBS call site does.
void SignedDigitDecompose(const std::shared_ptr<RingGSWCryptoParams>& params, const Input& input, Output& output) {
    auto QHalf{params->GetQ().ConvertToInt<BasicInteger>() >> 1};
    auto Q_int{params->GetQ().ConvertToInt<NativeInteger::SignedNativeInt>()};
    auto gBits{static_cast<NativeInteger::SignedNativeInt>(__builtin_ctz(params->GetBaseG()))};
    auto gBitsMaxBits{static_cast<NativeInteger::SignedNativeInt>(NativeInteger::MaxBits() - gBits)};
    const uint32_t digitsG2{(params->GetDigitsG() - 1) << 1};
    const uint32_t N{params->GetN()};

    for (uint32_t k{0}; k < N; ++k) {
        auto t0{input[0][k].ConvertToInt<BasicInteger>()};
        auto d0{static_cast<NativeInteger::SignedNativeInt>(t0 < QHalf ? t0 : t0 - Q_int)};
        auto t1{input[1][k].ConvertToInt<BasicInteger>()};
        auto d1{static_cast<NativeInteger::SignedNativeInt>(t1 < QHalf ? t1 : t1 - Q_int)};

        auto r0{(d0 << gBitsMaxBits) >> gBitsMaxBits};
        d0 = (d0 - r0) >> gBits;
        auto r1{(d1 << gBitsMaxBits) >> gBitsMaxBits};
        d1 = (d1 - r1) >> gBits;

        for (uint32_t d{0}; d < digitsG2; d += 2) {
            r0 = (d0 << gBitsMaxBits) >> gBitsMaxBits;
            d0 = (d0 - r0) >> gBits;
            if (r0 < 0)
                r0 += Q_int;
            output[d + 0][k] += r0;

            r1 = (d1 << gBitsMaxBits) >> gBitsMaxBits;
            d1 = (d1 - r1) >> gBits;
            if (r1 < 0)
                r1 += Q_int;
            output[d + 1][k] += r1;
        }
    }
}

// Low-overhead check over the extracted digits and quotient tail.
bool SignedDigitDecomposeOnlineChecked(const std::shared_ptr<RingGSWCryptoParams>& params, const Input& input,
                                      Output& output) {
    using CheckInt = int64_t;
    const auto qHalf = params->GetQ().ConvertToInt<BasicInteger>() >> 1;
    const auto qSigned = params->GetQ().ConvertToInt<NativeInteger::SignedNativeInt>();
    const auto gBits = static_cast<NativeInteger::SignedNativeInt>(__builtin_ctz(params->GetBaseG()));
    const auto gBitsMaxBits = static_cast<NativeInteger::SignedNativeInt>(NativeInteger::MaxBits() - gBits);
    const uint32_t digitsPerPoly = params->GetDigitsG() - 1;
    const CheckInt baseG = static_cast<CheckInt>(params->GetBaseG());
    const CheckInt halfBase = baseG >> 1;
    bool detected = false;

    for (uint32_t k = 0; k < params->GetN(); ++k) {
        const auto t0 = input[0][k].ConvertToInt<BasicInteger>();
        auto d0 = static_cast<NativeInteger::SignedNativeInt>(t0 < qHalf ? t0 : t0 - qSigned);
        const auto t1 = input[1][k].ConvertToInt<BasicInteger>();
        auto d1 = static_cast<NativeInteger::SignedNativeInt>(t1 < qHalf ? t1 : t1 - qSigned);

        auto r0 = (d0 << gBitsMaxBits) >> gBitsMaxBits;
        d0 = (d0 - r0) >> gBits;
        const CheckInt skippedQuot0 = static_cast<CheckInt>(d0);
        auto r1 = (d1 << gBitsMaxBits) >> gBitsMaxBits;
        d1 = (d1 - r1) >> gBits;
        const CheckInt skippedQuot1 = static_cast<CheckInt>(d1);
        CheckInt reconstruction0 = 0;
        CheckInt reconstruction1 = 0;
        CheckInt scale = 1;

        for (uint32_t digit = 0; digit < digitsPerPoly; ++digit) {
            r0 = (d0 << gBitsMaxBits) >> gBitsMaxBits;
            d0 = (d0 - r0) >> gBits;
            r1 = (d1 << gBitsMaxBits) >> gBitsMaxBits;
            d1 = (d1 - r1) >> gBits;
            if (r0 < 0)
                r0 += qSigned;
            if (r1 < 0)
                r1 += qSigned;
            using UnsignedNativeInt = std::make_unsigned_t<NativeInteger::SignedNativeInt>;
            auto rawR0 = static_cast<UnsignedNativeInt>(r0);
            auto rawR1 = static_cast<UnsignedNativeInt>(r1);
            InjectRaw(rawR0, FaultSite::DecomposeDigitExtraction, params->GetQ().GetMSB());
            InjectRaw(rawR1, FaultSite::DecomposeDigitExtraction, params->GetQ().GetMSB());
            r0 = static_cast<NativeInteger::SignedNativeInt>(rawR0);
            r1 = static_cast<NativeInteger::SignedNativeInt>(rawR1);
            output[digit << 1][k] += r0;
            output[(digit << 1) + 1][k] += r1;
            Inject(output[digit << 1][k], FaultSite::DecomposeIntermediateDigit, params->GetQ().GetMSB());
            Inject(output[(digit << 1) + 1][k], FaultSite::DecomposeIntermediateDigit, params->GetQ().GetMSB());
            const auto actual0 = output[digit << 1][k].ConvertToInt<NativeInteger::SignedNativeInt>();
            const auto actual1 = output[(digit << 1) + 1][k].ConvertToInt<NativeInteger::SignedNativeInt>();
            const CheckInt signed0 = static_cast<CheckInt>(actual0 > qSigned / 2 ? actual0 - qSigned : actual0);
            const CheckInt signed1 = static_cast<CheckInt>(actual1 > qSigned / 2 ? actual1 - qSigned : actual1);
            if (signed0 < -halfBase || signed0 >= halfBase || signed1 < -halfBase || signed1 >= halfBase)
                detected = true;
            reconstruction0 += scale * signed0;
            reconstruction1 += scale * signed1;
            scale <<= gBits;
        }
        if (reconstruction0 + scale * static_cast<CheckInt>(d0) != skippedQuot0 ||
            reconstruction1 + scale * static_cast<CheckInt>(d1) != skippedQuot1)
            detected = true;
    }
    return detected;
}

}  // namespace decompose

namespace keyswitch {

struct ABFTPrecomp {
    NativeVector weights;
    NativeVector weightPrecon;
    NativeInteger beta;
    NativeInteger betaPrecon;
    // chi[i][digit][j] = <weights, K_A[i][digit][j]> + beta*K_B[i][digit][j].
    std::vector<std::vector<std::vector<NativeInteger>>> chi;
};

struct Result {
    LWECiphertext output;
    bool detected;
};

uint32_t DigitCount(const std::shared_ptr<LWECryptoParams>& params) {
    const NativeInteger q(params->GetqKS());
    return static_cast<uint32_t>(std::ceil(std::log(q.ConvertToDouble()) /
                                           std::log(static_cast<double>(params->GetBaseKS()))));
}

LWECiphertext MakeRandomCiphertext(const std::shared_ptr<LWECryptoParams>& params, std::mt19937_64& prng) {
    const uint64_t q = params->GetqKS().ConvertToInt<uint64_t>();
    std::uniform_int_distribution<uint64_t> distribution(0, q - 1);
    NativeVector a(params->GetN(), params->GetqKS());
    for (uint32_t i = 0; i < params->GetN(); ++i)
        a[i] = NativeInteger(distribution(prng));
    return std::make_shared<LWECiphertextImpl>(std::move(a), NativeInteger(distribution(prng)));
}

bool CiphertextsEqual(const LWECiphertext& lhs, const LWECiphertext& rhs) {
    if (lhs->GetB() != rhs->GetB() || lhs->GetLength() != rhs->GetLength())
        return false;
    for (uint32_t k = 0; k < lhs->GetLength(); ++k) {
        if (lhs->GetA()[k] != rhs->GetA()[k])
            return false;
    }
    return true;
}

// A nonzero coefficient is not sufficient over a composite modulus: it may be
// a zero divisor and can annihilate a nonzero error.  Use Z_qKS units instead.
// For the standard parameter sets qKS is a power of two, so these are exactly
// the odd residues.
uint64_t SampleUnit(uint64_t q, std::mt19937_64& prng) {
    std::uniform_int_distribution<uint64_t> distribution(1, q - 1);
    for (;;) {
        const uint64_t candidate = distribution(prng);
        if (std::gcd(candidate, q) == 1)
            return candidate;
    }
}

ABFTPrecomp MakeABFTPrecomp(const std::shared_ptr<LWECryptoParams>& params, const LWESwitchingKey& key) {
    const uint32_t n = params->Getn();
    const uint32_t N = params->GetN();
    const uint32_t base = params->GetBaseKS();
    const uint32_t digits = DigitCount(params);
    const NativeInteger q(params->GetqKS());
    const uint64_t nativeQ = q.ConvertToInt<uint64_t>();
    std::mt19937_64 prng(20260409ULL);

    ABFTPrecomp pre{NativeVector(n, q), NativeVector(n, q), NativeInteger(SampleUnit(nativeQ, prng)), NativeInteger(0), {}};
    pre.betaPrecon = pre.beta.PrepModMulConst(q);
    for (uint32_t k = 0; k < n; ++k) {
        pre.weights[k] = NativeInteger(SampleUnit(nativeQ, prng));
        pre.weightPrecon[k] = pre.weights[k].PrepModMulConst(q);
    }

    pre.chi.resize(N);
    for (uint32_t i = 0; i < N; ++i) {
        pre.chi[i].resize(base);
        for (uint32_t digit = 0; digit < base; ++digit) {
            pre.chi[i][digit].resize(digits, NativeInteger(0));
            for (uint32_t j = 0; j < digits; ++j) {
                NativeInteger checksum(0);
                const auto& keyA = key->GetElementsA()[i][digit][j];
                for (uint32_t k = 0; k < n; ++k) {
                    auto term = keyA[k];
                    term.ModMulFastConstEq(pre.weights[k], q, pre.weightPrecon[k]);
                    checksum.ModAddFastEq(term, q);
                }
                auto scalar = key->GetElementsB()[i][digit][j];
                scalar.ModMulFastConstEq(pre.beta, q, pre.betaPrecon);
                checksum.ModAddFastEq(scalar, q);
                pre.chi[i][digit][j] = checksum;
            }
        }
    }
    return pre;
}

// Unprotected baseline copied from OpenFHE's LWEEncryptionScheme::KeySwitch.
LWECiphertext KeySwitch(const std::shared_ptr<LWECryptoParams>& params, const LWESwitchingKey& key,
                        const LWECiphertext& input) {
    const size_t n = params->Getn();
    const size_t N = params->GetN();
    NativeInteger q(params->GetqKS());
    NativeInteger::Integer base(params->GetBaseKS());
    const size_t digits = static_cast<size_t>(std::ceil(std::log(q.ConvertToDouble()) /
                                                         std::log(static_cast<double>(base))));
    NativeVector a(n, q);
    NativeInteger b(input->GetB());
    for (size_t i = 0; i < N; ++i) {
        const auto& refA = key->GetElementsA()[i];
        const auto& refB = key->GetElementsB()[i];
        NativeInteger::Integer atmp(input->GetA(i).ConvertToInt());
        for (size_t j = 0; j < digits; ++j) {
            const auto digit = atmp % base;
            atmp /= base;
            b.ModSubFastEq(refB[digit][j], q);
            const auto& selectedA = refA[digit][j];
            for (size_t k = 0; k < n; ++k)
                a[k].ModSubFastEq(selectedA[k], q);
        }
    }
    return std::make_shared<LWECiphertextImpl>(std::move(a), b);
}

// Fused KeySwitch ABFT. C_in consumes trusted, precomputed key-entry
// checksums in the original i/j order; C_out consumes the final (a,b).
Result KeySwitchABFT(const std::shared_ptr<LWECryptoParams>& params, const LWESwitchingKey& key,
                     const LWECiphertext& input, const ABFTPrecomp& pre) {
    const uint32_t n = params->Getn();
    const uint32_t N = params->GetN();
    const uint32_t digits = DigitCount(params);
    const NativeInteger q(params->GetqKS());
    const NativeInteger::Integer base(params->GetBaseKS());
    // qKS is a power of two for the evaluated parameter sets.  Values are
    // canonical residues in [0,qKS), so the modulus bit itself is not part of
    // their active representation.
    const uint32_t residueBitWidth = NativeInteger(q.ConvertToInt<uint64_t>() - 1).GetMSB();
    NativeVector a(n, q);
    NativeInteger b(input->GetB());
    NativeInteger cin(input->GetB());
    cin.ModMulFastConstEq(pre.beta, q, pre.betaPrecon);
    bool updateMismatch = false;

    for (uint32_t i = 0; i < N; ++i) {
        const auto& refA = key->GetElementsA()[i];
        const auto& refB = key->GetElementsB()[i];
        NativeInteger::Integer atmp(input->GetA(i).ConvertToInt());
        // Local DMR protects the digit-decomposition control state.  Only the
        // primary state is faulted; the replica is a trusted recomputation.
        NativeInteger::Integer atmpReplica(atmp);
        for (uint32_t j = 0; j < digits; ++j) {
            const uint32_t digit = static_cast<uint32_t>(atmp % base);
            atmp /= base;
            atmpReplica /= base;
            InjectRaw(atmp, FaultSite::KeySwitchUpdate, residueBitWidth);
            updateMismatch |= atmp != atmpReplica;
            cin.ModSubFastEq(pre.chi[i][digit][j], q);
            auto keyB = refB[digit][j];
            Inject(keyB, FaultSite::KeySwitchKeyLookup, residueBitWidth);
            b.ModSubFastEq(keyB, q);
            Inject(b, FaultSite::KeySwitchAccumulation, residueBitWidth);
            const auto& selectedA = refA[digit][j];
            for (uint32_t k = 0; k < n; ++k) {
                auto keyA = selectedA[k];
                Inject(keyA, FaultSite::KeySwitchKeyLookup, residueBitWidth);
                a[k].ModSubFastEq(keyA, q);
                Inject(a[k], FaultSite::KeySwitchAccumulation, residueBitWidth);
            }
        }
    }

    NativeInteger cout(b);
    cout.ModMulFastConstEq(pre.beta, q, pre.betaPrecon);
    for (uint32_t k = 0; k < n; ++k) {
        auto term = a[k];
        term.ModMulFastConstEq(pre.weights[k], q, pre.weightPrecon[k]);
        cout.ModAddFastEq(term, q);
    }
    return {std::make_shared<LWECiphertextImpl>(std::move(a), b), updateMismatch || cin != cout};
}

}  // namespace keyswitch

namespace vcmm {

// In this benchmark VCMM denotes the four y branches defined in the paper.
// Accumulator updates are deliberately outside this operator.
using Output = std::vector<NativePoly>;

struct Result {
    Output output;
    bool detected;
};

struct ChecksumPrecomp {
    std::vector<NativePoly> ev1sum;
    std::vector<NativePoly> ev2sum;
};

ChecksumPrecomp MakeChecksumPrecomp(const ConstRingGSWEvalKey& ek1, const ConstRingGSWEvalKey& ek2,
                                    uint32_t digitsG2) {
    ChecksumPrecomp pre;
    pre.ev1sum.reserve(digitsG2);
    pre.ev2sum.reserve(digitsG2);
    const auto& ev1 = ek1->GetElements();
    const auto& ev2 = ek2->GetElements();
    for (uint32_t i = 0; i < digitsG2; ++i) {
        pre.ev1sum.emplace_back(ev1[i][0] + ev1[i][1]);
        pre.ev2sum.emplace_back(ev2[i][0] + ev2[i][1]);
    }
    return pre;
}

bool OutputsEqual(const Output& lhs, const Output& rhs) {
    if (lhs.size() != rhs.size())
        return false;
    for (size_t i = 0; i < lhs.size(); ++i) {
        if (lhs[i] != rhs[i])
            return false;
    }
    return true;
}

// Unprotected baseline copied from OpenFHE's openfhe-development/src/binfhe/lib/rgsw-acc-cggi.cpp computation.  VCMM here is
// the four-branch y operator; accumulator updates are outside this operator.
Output VCMM(const std::vector<NativePoly>& dct, const ConstRingGSWEvalKey& ek1,
            const ConstRingGSWEvalKey& ek2, const NativePoly& monomial, const NativePoly& monomialNeg) {
    const uint32_t digitsG2 = dct.size();
    const auto& ev1 = ek1->GetElements();
    const auto& ev2 = ek2->GetElements();

    auto coreDct = dct;
    for (auto& poly : coreDct)
        for (uint32_t k = 0; k < poly.GetLength(); ++k)
            Inject(poly[k], FaultSite::VCMMInput, poly.GetModulus().GetMSB());

    const auto injectPoly = [](NativePoly& poly, FaultSite site) {
        for (uint32_t k = 0; k < poly.GetLength(); ++k)
            Inject(poly[k], site, poly.GetModulus().GetMSB());
    };

    NativePoly y0(coreDct[0] * ev1[0][0]);
    injectPoly(y0, FaultSite::VCMMMultiplication);
    for (uint32_t i = 1; i < digitsG2; ++i) {
        y0 += coreDct[i] * ev1[i][0];
        injectPoly(y0, FaultSite::VCMMAccumulation);
    }
    y0 *= monomial;
    injectPoly(y0, FaultSite::VCMMOutput);

    NativePoly y1(coreDct[0] * ev1[0][1]);
    injectPoly(y1, FaultSite::VCMMMultiplication);
    for (uint32_t i = 1; i < digitsG2; ++i) {
        y1 += coreDct[i] * ev1[i][1];
        injectPoly(y1, FaultSite::VCMMAccumulation);
    }
    y1 *= monomial;
    injectPoly(y1, FaultSite::VCMMOutput);

    NativePoly y2(coreDct[0] * ev2[0][0]);
    injectPoly(y2, FaultSite::VCMMMultiplication);
    for (uint32_t i = 1; i < digitsG2; ++i) {
        y2 += coreDct[i] * ev2[i][0];
        injectPoly(y2, FaultSite::VCMMAccumulation);
    }
    y2 *= monomialNeg;
    injectPoly(y2, FaultSite::VCMMOutput);

    NativePoly y3(coreDct[0] * ev2[0][1]);
    injectPoly(y3, FaultSite::VCMMMultiplication);
    for (uint32_t i = 1; i < digitsG2; ++i) {
        y3 += coreDct[i] * ev2[i][1];
        injectPoly(y3, FaultSite::VCMMAccumulation);
    }
    y3 *= monomialNeg;
    injectPoly(y3, FaultSite::VCMMOutput);
    return {std::move(y0), std::move(y1), std::move(y2), std::move(y3)};
}

// Preserve the branch-major VCMM body and compute both checksums as scalar
// dot products in the final coefficient pass, avoiding full Cin1/Cin2 arrays.
template <typename CheckAccumulator, bool ReduceCheckBeforeMonomial>
Result VCMMABFTBranchMajorCoefficientCheckImpl(const std::vector<NativePoly>& dct,
                                               const ConstRingGSWEvalKey& ek1,
                                               const ConstRingGSWEvalKey& ek2, const NativePoly& monomial,
                                               const NativePoly& monomialNeg, const ChecksumPrecomp& pre) {
    const uint32_t size = dct[0].GetLength();
    const uint64_t modulus = dct[0].GetModulus().ConvertToInt<uint64_t>();
    auto output = VCMM(dct, ek1, ek2, monomial, monomialNeg);
    bool detected = false;
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
        detected |= static_cast<uint64_t>(expected) != actual;
    }
    return {std::move(output), detected};
}

}  // namespace vcmm

FaultSite ParseFaultSite(const std::string& operation, const std::string& name) {
    if (operation == "ntt" || operation == "intt") {
        const std::string prefix = operation == "ntt" ? "NTT-" : "INTT-";
        if (name == "butterfly-load" || name == prefix + "butterfly-load") return FaultSite::ButterflyLoad;
        if (name == "butterfly-add" || name == prefix + "butterfly-add") return FaultSite::ButterflyAdd;
        if (name == "butterfly-sub" || name == prefix + "butterfly-sub") return FaultSite::ButterflySub;
        if (name == "butterfly-mul" || name == prefix + "butterfly-mul") return FaultSite::ButterflyMul;
        if (name == "butterfly-store" || name == prefix + "butterfly-store") return FaultSite::ButterflyStore;
    }
    else if (operation == "vcmm") {
        if (name == "input") return FaultSite::VCMMInput;
        if (name == "multiplication") return FaultSite::VCMMMultiplication;
        if (name == "accumulation") return FaultSite::VCMMAccumulation;
        if (name == "output") return FaultSite::VCMMOutput;
    }
    else if (operation == "decompose") {
        if (name == "digit-extraction") return FaultSite::DecomposeDigitExtraction;
        if (name == "intermediate-digit") return FaultSite::DecomposeIntermediateDigit;
    }
    else if (operation == "keyswitch") {
        if (name == "key-lookup") return FaultSite::KeySwitchKeyLookup;
        if (name == "accumulation") return FaultSite::KeySwitchAccumulation;
        if (name == "update") return FaultSite::KeySwitchUpdate;
    }
    throw std::invalid_argument("unsupported operator/site combination");
}

uint64_t ButterflyCandidates(uint32_t length, FaultSite site) {
    const uint64_t butterflies = static_cast<uint64_t>(length / 2) * GetMSB(length - 1);
    return (site == FaultSite::ButterflyLoad || site == FaultSite::ButterflyStore) ? 2 * butterflies : butterflies;
}

void RunTransformReliability(const std::string& operation, BINFHE_PARAMSET level,
                             const std::string& parameterName, FaultSite site, uint64_t samples,
                             uint32_t faultsPerSample, uint64_t seed) {
    BinFHEContext context;
    context.GenerateBinFHEContext(level);
    const auto poly = context.GetParams()->GetRingGSWParams()->GetPolyParams();
    const uint32_t cyclotomicOrder = poly->GetCyclotomicOrder();
    const uint32_t ringDimension = poly->GetRingDimension();
    const auto rootOfUnity = poly->GetRootOfUnity();
    const auto modulus = poly->GetModulus();
    intnat::ChineseRemainderTransformFTTNat<NativeVector>().PreCompute(rootOfUnity, cyclotomicOrder, modulus);
    using CRT = intnat::ChineseRemainderTransformFTTNat<NativeVector>;
    const auto& roots = CRT::m_rootOfUnityReverseTableByModulus[modulus];
    const auto& preconRoots = CRT::m_rootOfUnityPreconReverseTableByModulus[modulus];
    const auto& inverseRoots = CRT::m_rootOfUnityInverseReverseTableByModulus[modulus];
    const auto& preconInverseRoots = CRT::m_rootOfUnityInversePreconReverseTableByModulus[modulus];
    const uint32_t msb = GetMSB(ringDimension - 1);
    const auto& nInv = CRT::m_cycloOrderInverseTableByModulus[modulus][msb];
    const auto& preconNInv = CRT::m_cycloOrderInversePreconTableByModulus[modulus][msb];
    const auto weights = ntt::MakeABFTBitReverseWeights(rootOfUnity, cyclotomicOrder, modulus);
    const uint64_t candidates = ButterflyCandidates(ringDimension, site);
    std::mt19937_64 prng(seed);
    std::uniform_int_distribution<uint64_t> coefficient(0, modulus.ConvertToInt<uint64_t>() - 1);
    ReliabilityStats stats;

    for (uint64_t sample = 0; sample < samples; ++sample) {
        NativeVector input(ringDimension, modulus);
        for (uint32_t i = 0; i < ringDimension; ++i)
            input[i] = NativeInteger(coefficient(prng));
        NativeVector reference(input), checked(input);
        if (operation == "ntt")
            ntt::ForwardTransformToBitReverseInPlace(roots, preconRoots, &reference);
        else
            intt::InverseTransformFromBitReverseInPlace(inverseRoots, preconInverseRoots, nInv, preconNInv, &reference);

        if (sample == 0) {
            NativeVector clean(input);
            const bool falseAlarm = operation == "ntt"
                                        ? ntt::ForwardTransformABFT(roots, preconRoots, &clean, weights, level)
                                        : intt::InverseTransformABFT(inverseRoots, preconInverseRoots, nInv, preconNInv,
                                                                     &clean, weights, level);
            if (falseAlarm || !VectorsEqual(reference, clean))
                throw std::runtime_error("fault-free transform validation failed");
        }

        auto plan = MakeFaultPlan(site, candidates, faultsPerSample, prng);
        gFaultPlan = &plan;
        const bool detected = operation == "ntt"
                                  ? ntt::ForwardTransformABFT(roots, preconRoots, &checked, weights, level)
                                  : intt::InverseTransformABFT(inverseRoots, preconInverseRoots, nInv, preconNInv,
                                                               &checked, weights, level);
        gFaultPlan = nullptr;
        if (plan.applied != faultsPerSample)
            throw std::runtime_error("transform fault plan was not fully applied");
        stats.detectedSamples += detected;
        stats.maskedSamples += !detected && VectorsEqual(reference, checked);
    }
    PrintReliability(operation, parameterName, site, faultsPerSample, samples, stats);
}

void RunDecomposeReliability(BINFHE_PARAMSET level, const std::string& parameterName, FaultSite site,
                             uint64_t samples, uint32_t faultsPerSample, uint64_t seed) {
    BinFHEContext context;
    context.GenerateBinFHEContext(level);
    const auto params = context.GetParams()->GetRingGSWParams();
    const uint64_t candidates = static_cast<uint64_t>((params->GetDigitsG() - 1) << 1) * params->GetN();
    std::mt19937_64 prng(seed);
    ReliabilityStats stats;
    for (uint64_t sample = 0; sample < samples; ++sample) {
        const auto input = decompose::MakeRandomInput(params, prng);
        auto reference = decompose::MakeZeroOutput(params);
        auto checked = decompose::MakeZeroOutput(params);
        decompose::SignedDigitDecompose(params, input, reference);
        if (sample == 0) {
            auto clean = decompose::MakeZeroOutput(params);
            if (decompose::SignedDigitDecomposeOnlineChecked(params, input, clean) ||
                !decompose::OutputsEqual(reference, clean))
                throw std::runtime_error("fault-free Decompose validation failed");
        }
        auto plan = MakeFaultPlan(site, candidates, faultsPerSample, prng);
        gFaultPlan = &plan;
        const bool detected = decompose::SignedDigitDecomposeOnlineChecked(params, input, checked);
        gFaultPlan = nullptr;
        if (plan.applied != faultsPerSample)
            throw std::runtime_error("Decompose fault plan was not fully applied");
        stats.detectedSamples += detected;
        stats.maskedSamples += !detected && decompose::OutputsEqual(reference, checked);
    }
    PrintReliability("decompose", parameterName, site, faultsPerSample, samples, stats);
}

void RunKeySwitchReliability(BINFHE_PARAMSET level, const std::string& parameterName, FaultSite site,
                             uint64_t samples, uint32_t faultsPerSample, uint64_t seed) {
    BinFHEContext context;
    context.GenerateBinFHEContext(level);
    const auto secretKey = context.KeyGen();
    context.BTKeyGen(secretKey);
    const auto params = context.GetParams()->GetLWEParams();
    const auto key = context.GetSwitchKey();
    const auto pre = keyswitch::MakeABFTPrecomp(params, key);
    const uint64_t updates = static_cast<uint64_t>(params->GetN()) * keyswitch::DigitCount(params);
    const uint64_t candidates = site == FaultSite::KeySwitchUpdate ? updates : updates * (params->Getn() + 1);
    std::mt19937_64 prng(seed);
    ReliabilityStats stats;
    for (uint64_t sample = 0; sample < samples; ++sample) {
        const auto input = keyswitch::MakeRandomCiphertext(params, prng);
        const auto reference = keyswitch::KeySwitch(params, key, input);
        if (sample == 0) {
            const auto clean = keyswitch::KeySwitchABFT(params, key, input, pre);
            if (clean.detected || !keyswitch::CiphertextsEqual(reference, clean.output))
                throw std::runtime_error("fault-free KeySwitch validation failed");
        }
        auto plan = MakeFaultPlan(site, candidates, faultsPerSample, prng);
        gFaultPlan = &plan;
        const auto checked = keyswitch::KeySwitchABFT(params, key, input, pre);
        gFaultPlan = nullptr;
        if (plan.applied != faultsPerSample)
            throw std::runtime_error("KeySwitch fault plan was not fully applied");
        stats.detectedSamples += checked.detected;
        stats.maskedSamples += !checked.detected && keyswitch::CiphertextsEqual(reference, checked.output);
    }
    PrintReliability("keyswitch", parameterName, site, faultsPerSample, samples, stats);
}

void RunVCMMReliability(BINFHE_PARAMSET level, const std::string& parameterName, FaultSite site,
                        uint64_t samples, uint32_t faultsPerSample, uint64_t seed) {
    BinFHEContext context;
    context.GenerateBinFHEContext(level);
    const auto secretKey = context.KeyGen();
    context.BTKeyGen(secretKey);
    const auto params = context.GetParams();
    const auto lwe = params->GetLWEParams();
    const auto rgsw = params->GetRingGSWParams();
    const auto poly = rgsw->GetPolyParams();
    const auto refreshKey = context.GetRefreshKey();
    const uint32_t digitsG2 = (rgsw->GetDigitsG() - 1) << 1;
    const uint32_t n = rgsw->GetN();
    const auto protectedImpl = level == STD128
                                   ? vcmm::VCMMABFTBranchMajorCoefficientCheckImpl<uint64_t, true>
                                   : vcmm::VCMMABFTBranchMajorCoefficientCheckImpl<uint128_t, false>;
    std::mt19937_64 prng(seed);
    std::uniform_int_distribution<uint64_t> value(0, rgsw->GetQ().ConvertToInt<uint64_t>() - 1);
    std::vector<NativePoly> coefficientInput(2, NativePoly(poly, Format::COEFFICIENT, true));
    for (auto& component : coefficientInput)
        for (uint32_t k = 0; k < n; ++k)
            component[k] = NativeInteger(value(prng));
    std::vector<NativePoly> dct(digitsG2, NativePoly(poly, Format::COEFFICIENT, true));
    RingGSWAccumulatorDM accumulator;
    accumulator.SignedDigitDecompose(rgsw, coefficientInput, dct);
    for (auto& digit : dct)
        digit.SetFormat(Format::EVALUATION);

    const uint64_t candidates = site == FaultSite::VCMMInput ? static_cast<uint64_t>(digitsG2) * n
                                : site == FaultSite::VCMMAccumulation ? 4ULL * (digitsG2 - 1) * n
                                : 4ULL * n;
    std::uniform_int_distribution<uint32_t> keyIndex(0, lwe->Getn() - 1);
    std::uniform_int_distribution<uint32_t> monomialIndex(0, 2 * n - 1);
    ReliabilityStats stats;
    for (uint64_t sample = 0; sample < samples; ++sample) {
        const uint32_t index = keyIndex(prng);
        const uint32_t position = monomialIndex(prng);
        const auto& ek1 = (*refreshKey)[0][0][index];
        const auto& ek2 = (*refreshKey)[0][1][index];
        const auto& monomial = rgsw->GetMonomial(position);
        const auto& monomialNeg = rgsw->GetMonomial((2 * n - position) % (2 * n));
        const auto pre = vcmm::MakeChecksumPrecomp(ek1, ek2, digitsG2);
        const auto reference = vcmm::VCMM(dct, ek1, ek2, monomial, monomialNeg);
        if (sample == 0) {
            const auto clean = protectedImpl(dct, ek1, ek2, monomial, monomialNeg, pre);
            if (clean.detected || !vcmm::OutputsEqual(reference, clean.output))
                throw std::runtime_error("fault-free VCMM validation failed");
        }
        auto plan = MakeFaultPlan(site, candidates, faultsPerSample, prng);
        gFaultPlan = &plan;
        const auto checked = protectedImpl(dct, ek1, ek2, monomial, monomialNeg, pre);
        gFaultPlan = nullptr;
        if (plan.applied != faultsPerSample)
            throw std::runtime_error("VCMM fault plan was not fully applied");
        stats.detectedSamples += checked.detected;
        stats.maskedSamples += !checked.detected && vcmm::OutputsEqual(reference, checked.output);
    }
    PrintReliability("vcmm", parameterName, site, faultsPerSample, samples, stats);
}

BINFHE_PARAMSET ParseSecurity(const std::string& security) {
    if (security == "STD128")
        return STD128;
    if (security == "STD192")
        return STD192;
    if (security == "STD256")
        return STD256;
    throw std::invalid_argument("security must be STD128, STD192, or STD256");
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc < 6 || argc > 7) {
            std::cerr << "Usage: " << argv[0]
                      << " <ntt|intt|vcmm|decompose|keyswitch> <site> <STD128|STD192|STD256>"
                         " <samples> <faults_per_sample=1|2> [seed=7]\n";
            return EXIT_FAILURE;
        }
        const std::string operation(argv[1]);
        const std::string siteName(argv[2]);
        const std::string parameterName(argv[3]);
        const uint64_t samples = std::stoull(argv[4]);
        const uint32_t faultsPerSample = std::stoul(argv[5]);
        const uint64_t seed = argc == 7 ? std::stoull(argv[6]) : 7;
        if (samples == 0 || (faultsPerSample != 1 && faultsPerSample != 2))
            throw std::invalid_argument("samples must be nonzero and faults_per_sample must be 1 or 2");

        const auto level = ParseSecurity(parameterName);
        const auto site = ParseFaultSite(operation, siteName);
        if (operation == "ntt" || operation == "intt")
            RunTransformReliability(operation, level, parameterName, site, samples, faultsPerSample, seed);
        else if (operation == "decompose")
            RunDecomposeReliability(level, parameterName, site, samples, faultsPerSample, seed);
        else if (operation == "keyswitch")
            RunKeySwitchReliability(level, parameterName, site, samples, faultsPerSample, seed);
        else if (operation == "vcmm")
            RunVCMMReliability(level, parameterName, site, samples, faultsPerSample, seed);
        else
            throw std::invalid_argument("operator must be ntt, intt, vcmm, decompose, or keyswitch");
        return EXIT_SUCCESS;
    }
    catch (const std::exception& error) {
        gFaultPlan = nullptr;
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
