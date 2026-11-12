//==================================================================================
// Evaluation of NTT operator protection with ABFT checksums.
//
// This example evaluates the forward NTT used by PBS and its ABFT variants.
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
#include <stdexcept>
#include <string>
#include <vector>

using namespace lbcrypto;

namespace {

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

// Full NTT DMR: execute an independent replica of the complete NTT and
// compare its output with the primary execution.  A true return value denotes
// a mismatch.  This detects faults affecting only one replica; common-mode
// faults that identically affect both replicas are outside DMR's coverage.
bool ForwardTransformDMR(const NativeVector& rootOfUnityTable, const NativeVector& preconRootOfUnityTable,
                             NativeVector* element) {
    NativeVector replica(*element);
    ForwardTransformToBitReverseInPlace(rootOfUnityTable, preconRootOfUnityTable, element);
    ForwardTransformToBitReverseInPlace(rootOfUnityTable, preconRootOfUnityTable, &replica);
    return !VectorsEqual(*element, replica);
}

// Generate the output checksum weights in natural index order, as required by
// the unfused (naive) ABFT decoding pass.
NativeVector MakeABFTOriginalWeights(const NativeInteger& rootOfUnity, uint32_t cyclotomicOrder,
                                     const NativeInteger& modulus) {
    const uint32_t length  = cyclotomicOrder >> 1;
    const auto nInverse    = NativeInteger(length).ModInverse(modulus);
    const auto twoNInverse = NativeInteger(2).ModMul(nInverse, modulus);
    const auto rootInverse = rootOfUnity.ModInverse(modulus);
    const auto step        = rootInverse.ModMul(rootInverse, modulus);
    const auto mu          = modulus.ComputeMu();

    NativeVector weights(length, modulus);
    auto rootPower = rootInverse;
    for (uint32_t index = 0; index < length; ++index) {
        const auto denominator = NativeInteger(1).ModSub(rootPower, modulus);
        weights[index]         = twoNInverse.ModMul(denominator.ModInverse(modulus), modulus);
        rootPower.ModMulFastEq(step, modulus, mu);
    }
    return weights;
}

NativeVector MakeABFTPreconWeights(const NativeVector& weights) {
    const auto nativeModulus = NativeInteger(weights.GetModulus().ConvertToInt());
    NativeVector preconWeights(weights.GetLength(), nativeModulus);
    for (uint32_t index = 0; index < weights.GetLength(); ++index) {
        preconWeights[index] = NativeInteger(weights[index].ConvertToInt()).PrepModMulConst(nativeModulus);
    }
    return preconWeights;
}

// The input checksum is calculated before NTT and the output checksum is a
// separate, natural-order decoding pass after NTT.
bool ForwardTransformABFTNaive(const NativeVector& rootOfUnityTable, const NativeVector& preconRootOfUnityTable,
                               NativeVector* element, const NativeVector& aOrig, const NativeVector& preconAOrig) {
    const auto modulus = element->GetModulus();
    const uint32_t n   = element->GetLength() >> 1;
    auto cInObj                    = (*element)[0];
    const uint32_t length          = n << 1;
    for (uint32_t i = 1; i < length; ++i)
        cInObj.ModAddFastEq((*element)[i], modulus);
    const uint64_t cIn = static_cast<uint64_t>(cInObj.ConvertToInt());

    for (uint32_t m{1}, t{n}, logt{GetMSB(t)}; m < n; m <<= 1, t >>= 1, --logt) {
        for (uint32_t i{0}; i < m; ++i) {
            const auto omega       = rootOfUnityTable[i + m];
            const auto preconOmega = preconRootOfUnityTable[i + m];
            for (uint32_t j1{i << logt}, j2{j1 + t}; j1 < j2; ++j1) {
                auto omegaFactor = (*element)[j1 + t];
                omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
                auto loVal = (*element)[j1];
                auto hiVal = loVal + omegaFactor;
                if (hiVal >= modulus)
                    hiVal -= modulus;
                if (loVal < omegaFactor)
                    loVal += modulus;
                loVal -= omegaFactor;
                (*element)[j1]     = hiVal;
                (*element)[j1 + t] = loVal;
            }
        }
    }

    for (uint32_t i{0}; i < length; i += 2) {
        auto omegaFactor = (*element)[i + 1];
        const auto omega       = rootOfUnityTable[(i >> 1) + n];
        const auto preconOmega = preconRootOfUnityTable[(i >> 1) + n];
        omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
        auto loVal = (*element)[i];
        auto hiVal = loVal + omegaFactor;
        if (hiVal >= modulus)
            hiVal -= modulus;
        if (loVal < omegaFactor)
            loVal += modulus;
        loVal -= omegaFactor;
        (*element)[i]     = hiVal;
        (*element)[i + 1] = loVal;
    }

    const uint32_t bitWidth = GetMSB(element->GetLength() - 1);
    auto cOutObj            = NativeInteger{0};
    for (uint32_t i = 0; i < length; ++i) {
        const uint32_t reversedIndex = ReverseBits(i, bitWidth);
        auto term                    = (*element)[i];
        term.ModMulFastConstEq(aOrig[reversedIndex], modulus, preconAOrig[reversedIndex]);
        cOutObj.ModAddEq(term, modulus);
    }
    const uint64_t cOut = static_cast<uint64_t>(cOutObj.ConvertToInt());
    return cIn != cOut;
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
        omegaFactor.ModMulFastConstEq(rootOfUnityTable[(i >> 1) + n], modulus,
                                       preconRootOfUnityTable[(i >> 1) + n]);
        auto loVal = (*element)[i];
        auto hiVal = loVal + omegaFactor;
        if (hiVal >= modulus)
            hiVal -= modulus;
        if (loVal < omegaFactor)
            loVal += modulus;
        loVal -= omegaFactor;

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
            omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
            auto loVal = hiIn;
            auto hiVal = loVal + omegaFactor;
            if (hiVal >= modulus)
                hiVal -= modulus;
            if (loVal < omegaFactor)
                loVal += modulus;
            loVal -= omegaFactor;
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
                omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
                auto loVal = (*element)[j1];
                auto hiVal = loVal + omegaFactor;
                if (hiVal >= modulus)
                    hiVal -= modulus;
                if (loVal < omegaFactor)
                    loVal += modulus;
                loVal -= omegaFactor;
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

bool InverseTransformDMR(const NativeVector& inverseRoots, const NativeVector& preconInverseRoots,
                         const NativeInteger& nInv, const NativeInteger& preconNInv, NativeVector* element) {
    NativeVector replica(*element);
    InverseTransformFromBitReverseInPlace(inverseRoots, preconInverseRoots, nInv, preconNInv, element);
    InverseTransformFromBitReverseInPlace(inverseRoots, preconInverseRoots, nInv, preconNInv, &replica);
    return !VectorsEqual(*element, replica);
}

// Naive dual of NTT ABFT. ReverseBits is deliberately retained in the input
// encoding pass: aOrig is stored in natural, not bit-reversed, order.
bool InverseTransformABFTNaive(const NativeVector& inverseRoots, const NativeVector& preconInverseRoots,
                               const NativeInteger& nInv, const NativeInteger& preconNInv, NativeVector* element,
                               const NativeVector& aOrig, const NativeVector& preconAOrig) {
    const auto modulus = element->GetModulus();
    const uint32_t length = element->GetLength();
    const uint32_t bitWidth = GetMSB(length - 1);
    auto cInObj = NativeInteger{0};
    for (uint32_t i = 0; i < length; ++i) {
        const uint32_t reversedIndex = ReverseBits(i, bitWidth);
        auto term = (*element)[i];
        term.ModMulFastConstEq(aOrig[reversedIndex], modulus, preconAOrig[reversedIndex]);
        cInObj.ModAddEq(term, modulus);
    }

    InverseTransformFromBitReverseInPlace(inverseRoots, preconInverseRoots, nInv, preconNInv, element);

    auto cOutObj = NativeInteger{0};
    for (uint32_t i = 0; i < length; ++i)
        cOutObj.ModAddFastEq((*element)[i], modulus);
    return cInObj != cOutObj;
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
            auto omegaFactor = loVal;
            if (omegaFactor < hiVal)
                omegaFactor += modulus;
            omegaFactor -= hiVal;
            loVal += hiVal;
            if (loVal >= modulus)
                loVal -= modulus;
            omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
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
                auto omegaFactor = loVal;
                if (omegaFactor < hiVal)
                    omegaFactor += modulus;
                omegaFactor -= hiVal;
                loVal += hiVal;
                if (loVal >= modulus)
                    loVal -= modulus;
                omegaFactor.ModMulFastConstEq(omega, modulus, preconOmega);
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
        auto omegaFactor = loVal;
        if (omegaFactor < hiVal)
            omegaFactor += modulus;
        omegaFactor -= hiVal;
        loVal += hiVal;
        if (loVal >= modulus)
            loVal -= modulus;
        omegaFactor.ModMulFastConstEq(omega1Inv, modulus, preconOmega1Inv);
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

bool SignedDigitDecomposeDMR(const std::shared_ptr<RingGSWCryptoParams>& params, const Input& input, Output& output,
                             Output& replica) {
    SignedDigitDecompose(params, input, output);
    SignedDigitDecompose(params, input, replica);
    return !OutputsEqual(output, replica);
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
            const CheckInt signed0 = static_cast<CheckInt>(r0);
            const CheckInt signed1 = static_cast<CheckInt>(r1);
            if (signed0 < -halfBase || signed0 >= halfBase || signed1 < -halfBase || signed1 >= halfBase)
                detected = true;
            reconstruction0 += scale * signed0;
            reconstruction1 += scale * signed1;
            if (r0 < 0)
                r0 += qSigned;
            if (r1 < 0)
                r1 += qSigned;
            output[digit << 1][k] += r0;
            output[(digit << 1) + 1][k] += r1;
            scale <<= gBits;
        }
        if (reconstruction0 + scale * static_cast<CheckInt>(d0) != skippedQuot0 ||
            reconstruction1 + scale * static_cast<CheckInt>(d1) != skippedQuot1)
            detected = true;
    }
    return detected;
}

template <typename Transform>
void Benchmark(const char* method, const std::shared_ptr<RingGSWCryptoParams>& params, const std::string& parameterName,
               const std::vector<Input>& source, int warmup, int measure, Transform&& transform) {
    // Bounded batches keep the benchmark's memory use independent of measure.
    // Output reset is outside the timed interval because SignedDigitDecompose
    // receives a pre-zeroed dct buffer from its PBS call site.
    const Output zero = MakeZeroOutput(params);
    std::vector<Output> output(source.size(), zero);
    std::vector<Output> replica(source.size(), zero);
    const auto reset = [&](size_t count) {
        for (size_t i = 0; i < count; ++i) {
            output[i] = zero;
            replica[i] = zero;
        }
    };
    const auto execute = [&](int count, bool timed) {
        bool detected = false;
        const auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < count; ++i)
            detected |= transform(source[i], output[i], replica[i]);
        const auto elapsed = std::chrono::high_resolution_clock::now() - start;
        if (detected)
            throw std::runtime_error(std::string(method) + " detected a Decompose mismatch");
        return timed ? elapsed : std::chrono::high_resolution_clock::duration::zero();
    };

    for (int completed = 0; completed < warmup;) {
        const int count = std::min<int>(source.size(), warmup - completed);
        reset(count);
        execute(count, false);
        completed += count;
    }

    std::chrono::high_resolution_clock::duration elapsed{};
    for (int completed = 0; completed < measure;) {
        const int count = std::min<int>(source.size(), measure - completed);
        reset(count);
        elapsed += execute(count, true);
        completed += count;
    }
    const auto elapsedNs = std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count();
    std::cout << "Operator=decompose, Method=" << method << ", Params=" << parameterName
              << ", ringDim=" << params->GetN() << ", modulus=" << params->GetQ() << ", All time: " << elapsedNs
              << " ns, Average time: " << static_cast<double>(elapsedNs) / measure << " ns\n";
}

void Run(const std::shared_ptr<RingGSWCryptoParams>& params, const std::string& parameterName, int warmup, int measure) {
    constexpr size_t kBatchSize = 1024;
    const size_t batchSize = std::min(kBatchSize, static_cast<size_t>(warmup + measure));
    std::mt19937_64 prng(20260409ULL);
    std::vector<Input> source;
    source.reserve(batchSize);
    for (size_t i = 0; i < batchSize; ++i)
        source.emplace_back(MakeRandomInput(params, prng));

    const Output zero = MakeZeroOutput(params);
    // Validate every evaluated path on the whole input batch before timing.
    for (const auto& input : source) {
        Output reference = zero;
        SignedDigitDecompose(params, input, reference);
        Output dmr = zero;
        Output dmrReplica = zero;
        const bool dmrMismatch = SignedDigitDecomposeDMR(params, input, dmr, dmrReplica);
        Output online = zero;
        const bool onlineMismatch = SignedDigitDecomposeOnlineChecked(params, input, online);
        if (dmrMismatch || onlineMismatch || !OutputsEqual(reference, dmr) || !OutputsEqual(reference, online))
            throw std::runtime_error("Decompose validation failed");
    }

    Benchmark("No-Protection", params, parameterName, source, warmup, measure,
              [&](const Input& in, Output& out, Output&) {
                  SignedDigitDecompose(params, in, out);
                  return false;
              });
    Benchmark("DMR", params, parameterName, source, warmup, measure,
              [&](const Input& in, Output& out, Output& replica) {
                  return SignedDigitDecomposeDMR(params, in, out, replica);
              });
    Benchmark("Ours", params, parameterName, source, warmup, measure,
              [&](const Input& in, Output& out, Output&) {
                  return SignedDigitDecomposeOnlineChecked(params, in, out);
              });
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

uint64_t DigestCiphertext(const LWECiphertext& ciphertext) {
    uint64_t digest = ciphertext->GetB().ConvertToInt<uint64_t>();
    for (uint32_t k = 0; k < ciphertext->GetLength(); ++k)
        digest = (digest * 1315423911ULL) ^ ciphertext->GetA()[k].ConvertToInt<uint64_t>();
    return digest;
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

Result KeySwitchDMR(const std::shared_ptr<LWECryptoParams>& params, const LWESwitchingKey& key,
                    const LWECiphertext& input) {
    auto primary = KeySwitch(params, key, input);
    const auto replica = KeySwitch(params, key, input);
    const bool detected = !CiphertextsEqual(primary, replica);
    return {std::move(primary), detected};
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
    NativeVector a(n, q);
    NativeInteger b(input->GetB());
    NativeInteger cin(input->GetB());
    cin.ModMulFastConstEq(pre.beta, q, pre.betaPrecon);

    for (uint32_t i = 0; i < N; ++i) {
        const auto& refA = key->GetElementsA()[i];
        const auto& refB = key->GetElementsB()[i];
        NativeInteger::Integer atmp(input->GetA(i).ConvertToInt());
        for (uint32_t j = 0; j < digits; ++j) {
            const uint32_t digit = static_cast<uint32_t>(atmp % base);
            atmp /= base;
            cin.ModSubFastEq(pre.chi[i][digit][j], q);
            b.ModSubFastEq(refB[digit][j], q);
            const auto& selectedA = refA[digit][j];
            for (uint32_t k = 0; k < n; ++k)
                a[k].ModSubFastEq(selectedA[k], q);
        }
    }

    NativeInteger cout(b);
    cout.ModMulFastConstEq(pre.beta, q, pre.betaPrecon);
    for (uint32_t k = 0; k < n; ++k) {
        auto term = a[k];
        term.ModMulFastConstEq(pre.weights[k], q, pre.weightPrecon[k]);
        cout.ModAddFastEq(term, q);
    }
    return {std::make_shared<LWECiphertextImpl>(std::move(a), b), cin != cout};
}

template <typename Transform>
void Benchmark(const char* method, const std::shared_ptr<LWECryptoParams>& params, const std::string& parameterName,
               const std::vector<LWECiphertext>& source, int warmup, int measure, Transform&& transform) {
    volatile uint64_t sink = 0;
    const auto execute = [&](int count, bool timed) {
        bool detected = false;
        const auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < count; ++i) {
            auto result = transform(source[i]);
            sink ^= DigestCiphertext(result.output);
            detected |= result.detected;
        }
        const auto elapsed = std::chrono::high_resolution_clock::now() - start;
        if (detected)
            throw std::runtime_error(std::string(method) + " detected a KeySwitch mismatch");
        return timed ? elapsed : std::chrono::high_resolution_clock::duration::zero();
    };

    for (int completed = 0; completed < warmup;) {
        const int count = std::min<int>(source.size(), warmup - completed);
        execute(count, false);
        completed += count;
    }
    std::chrono::high_resolution_clock::duration elapsed{};
    for (int completed = 0; completed < measure;) {
        const int count = std::min<int>(source.size(), measure - completed);
        elapsed += execute(count, true);
        completed += count;
    }
    const auto elapsedNs = std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count();
    std::cout << "Operator=keyswitch, Method=" << method << ", Params=" << parameterName
              << ", n=" << params->Getn() << ", N=" << params->GetN() << ", qKS=" << params->GetqKS()
              << ", baseKS=" << params->GetBaseKS() << ", digitCount=" << DigitCount(params)
              << ", All time: " << elapsedNs << " ns, Average time: " << std::fixed << std::setprecision(3)
              << static_cast<double>(elapsedNs) / measure
              << " ns, Digest sink: " << sink << '\n';
}

void Run(BinFHEContext& context, const std::string& parameterName, int warmup, int measure) {
    const auto secretKey = context.KeyGen();
    context.BTKeyGen(secretKey);
    const auto params = context.GetParams()->GetLWEParams();
    const auto key = context.GetSwitchKey();
    const auto pre = MakeABFTPrecomp(params, key);

    constexpr size_t kBatchSize = 64;
    const size_t batchSize = std::min(kBatchSize, static_cast<size_t>(warmup + measure));
    std::mt19937_64 prng(20260409ULL);
    std::vector<LWECiphertext> source;
    source.reserve(batchSize);
    for (size_t i = 0; i < batchSize; ++i)
        source.emplace_back(MakeRandomCiphertext(params, prng));

    // No-Protection, DMR, and ABFT must return exactly the same ciphertext
    // in fault-free execution; validate every reused input before timing.
    for (const auto& input : source) {
        const auto reference = KeySwitch(params, key, input);
        const auto dmr = KeySwitchDMR(params, key, input);
        const auto ours = KeySwitchABFT(params, key, input, pre);
        if (dmr.detected || ours.detected || !CiphertextsEqual(reference, dmr.output) ||
            !CiphertextsEqual(reference, ours.output))
            throw std::runtime_error("KeySwitch validation failed");
    }

    Benchmark("No-Protection", params, parameterName, source, warmup, measure,
              [&](const LWECiphertext& input) { return Result{KeySwitch(params, key, input), false}; });
    Benchmark("DMR", params, parameterName, source, warmup, measure,
              [&](const LWECiphertext& input) { return KeySwitchDMR(params, key, input); });
    Benchmark("Ours", params, parameterName, source, warmup, measure,
              [&](const LWECiphertext& input) { return KeySwitchABFT(params, key, input, pre); });
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

uint64_t DigestOutput(const Output& output) {
    uint64_t digest = 0;
    for (const auto& poly : output) {
        for (uint32_t i = 0; i < poly.GetLength(); ++i)
            digest = (digest * 1315423911ULL) ^ poly[i].ConvertToInt<uint64_t>();
    }
    return digest;
}

// Unprotected baseline copied from OpenFHE's openfhe-development/src/binfhe/lib/rgsw-acc-cggi.cpp computation.  VCMM here is
// the four-branch y operator; accumulator updates are outside this operator.
Output VCMM(const std::vector<NativePoly>& dct, const ConstRingGSWEvalKey& ek1,
            const ConstRingGSWEvalKey& ek2, const NativePoly& monomial, const NativePoly& monomialNeg) {
    const uint32_t digitsG2 = dct.size();
    const auto& ev1 = ek1->GetElements();
    const auto& ev2 = ek2->GetElements();

    NativePoly y0(dct[0] * ev1[0][0]);
    for (uint32_t i = 1; i < digitsG2; ++i)
        y0 += dct[i] * ev1[i][0];
    y0 *= monomial;

    NativePoly y1(dct[0] * ev1[0][1]);
    for (uint32_t i = 1; i < digitsG2; ++i)
        y1 += dct[i] * ev1[i][1];
    y1 *= monomial;

    NativePoly y2(dct[0] * ev2[0][0]);
    for (uint32_t i = 1; i < digitsG2; ++i)
        y2 += dct[i] * ev2[i][0];
    y2 *= monomialNeg;

    NativePoly y3(dct[0] * ev2[0][1]);
    for (uint32_t i = 1; i < digitsG2; ++i)
        y3 += dct[i] * ev2[i][1];
    y3 *= monomialNeg;
    return {std::move(y0), std::move(y1), std::move(y2), std::move(y3)};
}

Result VCMMDMR(const std::vector<NativePoly>& dct, const ConstRingGSWEvalKey& ek1,
               const ConstRingGSWEvalKey& ek2, const NativePoly& monomial, const NativePoly& monomialNeg) {
    auto primary = VCMM(dct, ek1, ek2, monomial, monomialNeg);
    const auto replica = VCMM(dct, ek1, ek2, monomial, monomialNeg);
    const bool detected = !OutputsEqual(primary, replica);
    return {std::move(primary), detected};
}

// NaiveCheckSumVCMM: encode with u=(1,1,1,1), run the same VCMM core, and
// decode by summing the four y branches.  The checksum intentionally is not
// precomputed or fused.
Result VCMMNaiveABFT(const std::vector<NativePoly>& dct, const ConstRingGSWEvalKey& ek1,
                     const ConstRingGSWEvalKey& ek2, const NativePoly& monomial,
                     const NativePoly& monomialNeg) {
    const uint32_t digitsG2 = dct.size();
    const auto& ev1 = ek1->GetElements();
    const auto& ev2 = ek2->GetElements();
    NativePoly cin1(dct[0] * (ev1[0][0] + ev1[0][1]));
    NativePoly cin2(dct[0] * (ev2[0][0] + ev2[0][1]));
    for (uint32_t i = 1; i < digitsG2; ++i) {
        cin1 += dct[i] * (ev1[i][0] + ev1[i][1]);
        cin2 += dct[i] * (ev2[i][0] + ev2[i][1]);
    }
    NativePoly expected(cin1 * monomial);
    expected += cin2 * monomialNeg;

    auto output = VCMM(dct, ek1, ek2, monomial, monomialNeg);
    NativePoly actual(output[0]);
    actual += output[1];
    actual += output[2];
    actual += output[3];
    return {std::move(output), expected != actual};
}

template <typename CoreAccumulator, typename CheckAccumulator, bool ReduceCheckBeforeMonomial>
bool FinalizeFusedCoefficient(CoreAccumulator y0Raw, CoreAccumulator y1Raw, CoreAccumulator y2Raw,
                              CoreAccumulator y3Raw, CheckAccumulator cin1, CheckAccumulator cin2,
                              uint64_t monomial, uint64_t monomialNeg, uint64_t modulus, NativeInteger& y0,
                              NativeInteger& y1, NativeInteger& y2, NativeInteger& y3) {
    const auto coreModulus = static_cast<CoreAccumulator>(modulus);
    const uint64_t y0Value = static_cast<uint64_t>(((y0Raw % coreModulus) * monomial) % coreModulus);
    const uint64_t y1Value = static_cast<uint64_t>(((y1Raw % coreModulus) * monomial) % coreModulus);
    const uint64_t y2Value = static_cast<uint64_t>(((y2Raw % coreModulus) * monomialNeg) % coreModulus);
    const uint64_t y3Value = static_cast<uint64_t>(((y3Raw % coreModulus) * monomialNeg) % coreModulus);
    y0 = NativeInteger(y0Value);
    y1 = NativeInteger(y1Value);
    y2 = NativeInteger(y2Value);
    y3 = NativeInteger(y3Value);

    if constexpr (ReduceCheckBeforeMonomial) {
        cin1 %= modulus;
        cin2 %= modulus;
    }
    const auto expected = (cin1 * static_cast<CheckAccumulator>(monomial) +
                           cin2 * static_cast<CheckAccumulator>(monomialNeg)) %
                          modulus;
    const uint64_t actual = (y0Value + y1Value + y2Value + y3Value) % modulus;
    return static_cast<uint64_t>(expected) != actual;
}

// Experimental candidate: fuse the VCMM core and checksum computations by
// coefficient and use lazy reduction for the core itself.  It is retained for
// future comparison, but is not the Ours path used by this benchmark.
template <typename CoreAccumulator, typename CheckAccumulator, bool ReduceCheckBeforeMonomial>
Result VCMMABFTFusedCoefficientImpl(const std::vector<NativePoly>& dct, const ConstRingGSWEvalKey& ek1,
                                    const ConstRingGSWEvalKey& ek2, const NativePoly& monomial,
                                    const NativePoly& monomialNeg, const ChecksumPrecomp& pre) {
    const uint32_t size = dct[0].GetLength();
    const uint64_t modulus = dct[0].GetModulus().ConvertToInt<uint64_t>();
    const auto& ev1 = ek1->GetElements();
    const auto& ev2 = ek2->GetElements();
    Output output;
    output.reserve(4);
    for (uint32_t branch = 0; branch < 4; ++branch)
        output.emplace_back(dct[0].GetParams(), Format::EVALUATION, true);

    bool detected = false;
    for (uint32_t j = 0; j < size; ++j) {
        CoreAccumulator y0Raw = 0;
        CoreAccumulator y1Raw = 0;
        CoreAccumulator y2Raw = 0;
        CoreAccumulator y3Raw = 0;
        CheckAccumulator cin1 = 0;
        CheckAccumulator cin2 = 0;
        for (uint32_t i = 0; i < dct.size(); ++i) {
            const uint64_t digit = dct[i][j].ConvertToInt<uint64_t>();
            y0Raw += static_cast<CoreAccumulator>(digit) * ev1[i][0][j].ConvertToInt<uint64_t>();
            y1Raw += static_cast<CoreAccumulator>(digit) * ev1[i][1][j].ConvertToInt<uint64_t>();
            y2Raw += static_cast<CoreAccumulator>(digit) * ev2[i][0][j].ConvertToInt<uint64_t>();
            y3Raw += static_cast<CoreAccumulator>(digit) * ev2[i][1][j].ConvertToInt<uint64_t>();
            cin1 += static_cast<CheckAccumulator>(digit) * pre.ev1sum[i][j].ConvertToInt<uint64_t>();
            cin2 += static_cast<CheckAccumulator>(digit) * pre.ev2sum[i][j].ConvertToInt<uint64_t>();
        }
        detected |= FinalizeFusedCoefficient<CoreAccumulator, CheckAccumulator, ReduceCheckBeforeMonomial>(
            y0Raw, y1Raw, y2Raw, y3Raw, cin1, cin2, monomial[j].ConvertToInt<uint64_t>(),
            monomialNeg[j].ConvertToInt<uint64_t>(), modulus, output[0][j], output[1][j], output[2][j], output[3][j]);
    }
    return {std::move(output), detected};
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

template <typename Transform>
void Benchmark(const char* method, const std::string& parameterName, const std::vector<uint32_t>& keyIndexes,
               const std::vector<NativePoly>& monomials, const std::vector<NativePoly>& monomialNegs,
               const std::vector<ChecksumPrecomp>& pre, const std::vector<NativePoly>& dct,
               const ConstRingGSWACCKey& refreshKey, int warmup, int measure,
               Transform&& transform) {
    volatile uint64_t sink = 0;
    const auto execute = [&](int count, bool timed) {
        bool detected = false;
        const auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < count; ++i) {
            const size_t input = static_cast<size_t>(i) % keyIndexes.size();
            const auto& ek1 = (*refreshKey)[0][0][keyIndexes[input]];
            const auto& ek2 = (*refreshKey)[0][1][keyIndexes[input]];
            const auto result = transform(dct, ek1, ek2, monomials[input], monomialNegs[input], pre[input]);
            sink ^= DigestOutput(result.output);
            detected |= result.detected;
        }
        const auto elapsed = std::chrono::high_resolution_clock::now() - start;
        if (detected)
            throw std::runtime_error(std::string(method) + " detected a VCMM mismatch");
        return timed ? elapsed : std::chrono::high_resolution_clock::duration::zero();
    };

    for (int completed = 0; completed < warmup;) {
        const int count = std::min<int>(keyIndexes.size(), warmup - completed);
        execute(count, false);
        completed += count;
    }
    std::chrono::high_resolution_clock::duration elapsed{};
    for (int completed = 0; completed < measure;) {
        const int count = std::min<int>(keyIndexes.size(), measure - completed);
        elapsed += execute(count, true);
        completed += count;
    }
    const auto elapsedNs = std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count();
    std::cout << "Operator=vcmm, Method=" << method << ", Params=" << parameterName
              << ", digitsG2=" << dct.size() << ", ringDim=" << dct[0].GetLength()
              << ", All time: " << elapsedNs << " ns, Average time: " << std::fixed << std::setprecision(3)
              << static_cast<double>(elapsedNs) / measure << " ns, Digest sink: " << sink << '\n';
}

void Run(BinFHEContext& context, BINFHE_PARAMSET level, const std::string& parameterName, int warmup, int measure) {
    const auto secretKey = context.KeyGen();
    context.BTKeyGen(secretKey);
    const auto params = context.GetParams();
    const auto lwe = params->GetLWEParams();
    const auto rgsw = params->GetRingGSWParams();
    const auto poly = rgsw->GetPolyParams();
    const ConstRingGSWACCKey refreshKey = context.GetRefreshKey();
    const uint32_t digitsG2 = (rgsw->GetDigitsG() - 1) << 1;
    const uint32_t keyCount = lwe->Getn();
    const uint32_t monomialModulus = 2 * rgsw->GetN();
    // Retain the branch-major NativePoly core and avoid materializing full
    // checksum arrays by completing both checksum dot products per coefficient.
    const auto protectedImpl =
        level == STD128 ? VCMMABFTBranchMajorCoefficientCheckImpl<uint64_t, true>
                        : VCMMABFTBranchMajorCoefficientCheckImpl<uint128_t, false>;

    // Generate one decomposed input shared by all methods.  None of the VCMM
    // implementations modifies dct, so duplicating it is unnecessary.
    DiscreteUniformGeneratorImpl<NativeVector> dug(lwe->Getq());
    NativeVector m0(lwe->GetN(), lwe->GetQ());
    NativeVector m1(lwe->GetN(), lwe->GetQ());
    for (uint32_t i = 0; i < lwe->GetN(); ++i) {
        m0[i] = dug.GenerateInteger();
        m1[i] = dug.GenerateInteger();
    }
    std::vector<NativePoly> accumulator(2);
    accumulator[0] = NativePoly(poly, Format::COEFFICIENT, false);
    accumulator[0].SetValues(std::move(m0), Format::COEFFICIENT);
    accumulator[0].SetFormat(Format::EVALUATION);
    accumulator[1] = NativePoly(poly, Format::COEFFICIENT, false);
    accumulator[1].SetValues(std::move(m1), Format::COEFFICIENT);
    accumulator[1].SetFormat(Format::EVALUATION);
    std::vector<NativePoly> dct(digitsG2, NativePoly(poly, Format::COEFFICIENT, true));
    RingGSWAccumulatorDM accumulatorDM;
    accumulatorDM.SignedDigitDecompose(rgsw, accumulator, dct);
    for (auto& digit : dct)
        digit.SetFormat(Format::EVALUATION);

    constexpr size_t kBatchSize = 64;
    const size_t batchSize = std::min(kBatchSize, static_cast<size_t>(warmup + measure));
    std::vector<uint32_t> keyIndexes;
    std::vector<NativePoly> monomials;
    std::vector<NativePoly> monomialNegs;
    std::vector<ChecksumPrecomp> pre;
    keyIndexes.reserve(batchSize);
    monomials.reserve(batchSize);
    monomialNegs.reserve(batchSize);
    for (size_t i = 0; i < batchSize; ++i) {
        const NativeInteger a = dug.GenerateInteger();
        const uint32_t position = a.ConvertToInt<uint32_t>();
        const uint32_t neg = NativeInteger(0).ModSubFast(a, NativeInteger(monomialModulus)).ConvertToInt<uint32_t>();
        keyIndexes.emplace_back(i % keyCount);
        monomials.emplace_back(rgsw->GetMonomial(position == monomialModulus ? 0 : position));
        monomialNegs.emplace_back(rgsw->GetMonomial(neg == monomialModulus ? 0 : neg));
    }
    pre.reserve(batchSize);
    for (const auto keyIndex : keyIndexes)
        pre.emplace_back(MakeChecksumPrecomp((*refreshKey)[0][0][keyIndex], (*refreshKey)[0][1][keyIndex], digitsG2));

    // Validate all reused inputs before timing: protection must leave the
    // fault-free four-branch VCMM output unchanged.
    for (size_t i = 0; i < batchSize; ++i) {
        const auto& ek1 = (*refreshKey)[0][0][keyIndexes[i]];
        const auto& ek2 = (*refreshKey)[0][1][keyIndexes[i]];
        const auto reference = VCMM(dct, ek1, ek2, monomials[i], monomialNegs[i]);
        const auto dmr = VCMMDMR(dct, ek1, ek2, monomials[i], monomialNegs[i]);
        const auto naive = VCMMNaiveABFT(dct, ek1, ek2, monomials[i], monomialNegs[i]);
        const auto ours = protectedImpl(dct, ek1, ek2, monomials[i], monomialNegs[i], pre[i]);
        if (dmr.detected || naive.detected || ours.detected || !OutputsEqual(reference, dmr.output) ||
            !OutputsEqual(reference, naive.output) || !OutputsEqual(reference, ours.output))
            throw std::runtime_error("VCMM validation failed");
    }

    Benchmark("No-Protection", parameterName, keyIndexes, monomials, monomialNegs, pre, dct, refreshKey, warmup, measure,
              [](const auto& input, const auto& ek1, const auto& ek2, const auto& monomial, const auto& monomialNeg,
                 const auto&) {
                  return Result{VCMM(input, ek1, ek2, monomial, monomialNeg), false};
              });
    Benchmark("DMR", parameterName, keyIndexes, monomials, monomialNegs, pre, dct, refreshKey, warmup, measure,
              [](const auto& input, const auto& ek1, const auto& ek2, const auto& monomial, const auto& monomialNeg,
                 const auto&) {
                  return VCMMDMR(input, ek1, ek2, monomial, monomialNeg);
              });
    Benchmark("Naive-ABFT", parameterName, keyIndexes, monomials, monomialNegs, pre, dct, refreshKey, warmup, measure,
              [](const auto& input, const auto& ek1, const auto& ek2, const auto& monomial, const auto& monomialNeg,
                 const auto&) {
                  return VCMMNaiveABFT(input, ek1, ek2, monomial, monomialNeg);
              });
    Benchmark("Ours", parameterName, keyIndexes, monomials, monomialNegs, pre, dct, refreshKey, warmup, measure,
              [&](const auto& input, const auto& ek1, const auto& ek2, const auto& monomial, const auto& monomialNeg,
                  const auto& checksumPre) {
                  return protectedImpl(input, ek1, ek2, monomial, monomialNeg, checksumPre);
              });
}

}  // namespace vcmm

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
    const std::string operation = argc > 1 ? argv[1] : "ntt";
    const std::string params = argc > 2 ? argv[2] : "STD256";
    const auto level = ParseSecurity(params);
    const int warmup = argc > 3 ? std::atoi(argv[3]) : 10000;
    const int measure = argc > 4 ? std::atoi(argv[4]) : 100000;
    if (argc > 5 || (operation != "ntt" && operation != "intt" && operation != "decompose" && operation != "keyswitch" &&
                     operation != "vcmm") ||
        warmup < 0 || measure <= 0) {
        std::cerr << "Usage: " << argv[0]
                  << " [ntt|intt|decompose|keyswitch|vcmm] [STD128|STD192|STD256] [warmup] [measure]\n";
        return 1;
    }

    BinFHEContext context;
    context.GenerateBinFHEContext(level);
    if (operation == "decompose") {
        decompose::Run(context.GetParams()->GetRingGSWParams(), params, warmup, measure);
        return 0;
    }
    if (operation == "keyswitch") {
        keyswitch::Run(context, params, warmup, measure);
        return 0;
    }
    if (operation == "vcmm") {
        vcmm::Run(context, level, params, warmup, measure);
        return 0;
    }
    const auto poly = context.GetParams()->GetRingGSWParams()->GetPolyParams();
    const auto cyclotomicOrder = poly->GetCyclotomicOrder();
    const auto ringDimension   = poly->GetRingDimension();
    const auto rootOfUnity     = poly->GetRootOfUnity();
    const auto modulus         = poly->GetModulus();

    const int total = warmup + measure;
    std::vector<NativeVector> sourcePolys;
    sourcePolys.reserve(total);
    DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(modulus);
    for (int k = 0; k < total; ++k) {
        NativeVector polyInput(ringDimension, modulus);
        for (uint32_t i = 0; i < ringDimension; ++i)
            polyInput[i] = dug.GenerateInteger();
        sourcePolys.emplace_back(std::move(polyInput));
    }
    // Both transforms are in-place. Preserve this source batch and refresh
    // the reusable working batch before each method outside its timed interval.
    std::vector<NativeVector> workingPolys(sourcePolys);

    // Precompute the transform tables used by the evaluated bodies.
    intnat::ChineseRemainderTransformFTTNat<NativeVector>().PreCompute(rootOfUnity, cyclotomicOrder, modulus);
    using CRT = intnat::ChineseRemainderTransformFTTNat<NativeVector>;
    const auto& rootTable =
        CRT::m_rootOfUnityReverseTableByModulus[modulus];
    const auto& preconRootTable =
        CRT::m_rootOfUnityPreconReverseTableByModulus[modulus];
    const auto& inverseRootTable =
        CRT::m_rootOfUnityInverseReverseTableByModulus[modulus];
    const auto& preconInverseRootTable =
        CRT::m_rootOfUnityInversePreconReverseTableByModulus[modulus];
    const uint32_t msb = GetMSB(ringDimension - 1);
    const auto& nInv = CRT::m_cycloOrderInverseTableByModulus[modulus][msb];
    const auto& preconNInv = CRT::m_cycloOrderInversePreconTableByModulus[modulus][msb];
    const auto aOrig = ntt::MakeABFTOriginalWeights(rootOfUnity, cyclotomicOrder, modulus);
    const auto preconAOrig = ntt::MakeABFTPreconWeights(aOrig);
    const auto aBitReverse = ntt::MakeABFTBitReverseWeights(rootOfUnity, cyclotomicOrder, modulus);

    const auto benchmark = [&](const char* method, auto&& transform) {
        workingPolys = sourcePolys;
        bool detectedMismatch = false;
        for (int k = 0; k < warmup; ++k)
            detectedMismatch |= transform(&workingPolys[k]);

        const auto start = std::chrono::high_resolution_clock::now();
        for (int k = warmup; k < total; ++k)
            detectedMismatch |= transform(&workingPolys[k]);
        const auto end = std::chrono::high_resolution_clock::now();
        const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        if (detectedMismatch)
            throw std::runtime_error(std::string(method) + " detected a " + operation + " mismatch");

        std::cout << "Operator=" << operation << ", Method=" << method << ", Params=" << params << ", ringDim=" << ringDimension
                  << ", cycloOrder=" << cyclotomicOrder << ", modulus=" << modulus << ", All time: " << elapsed
                  << " ns, Average time: " << static_cast<double>(elapsed) / measure << " ns\n";
    };

    if (operation == "ntt") {
        NativeVector reference(sourcePolys.front());
        ntt::ForwardTransformToBitReverseInPlace(rootTable, preconRootTable, &reference);
        NativeVector dmrResult(sourcePolys.front());
        const bool dmrMismatch = ntt::ForwardTransformDMR(rootTable, preconRootTable, &dmrResult);
        NativeVector naiveResult(sourcePolys.front());
        const bool naiveMismatch = ntt::ForwardTransformABFTNaive(rootTable, preconRootTable, &naiveResult, aOrig, preconAOrig);
        NativeVector oursResult(sourcePolys.front());
        const bool oursMismatch = ntt::ForwardTransformABFT(rootTable, preconRootTable, &oursResult, aBitReverse, level);
        if (dmrMismatch || naiveMismatch || oursMismatch || !VectorsEqual(reference, dmrResult) ||
            !VectorsEqual(reference, naiveResult) || !VectorsEqual(reference, oursResult)) {
            std::cerr << "NTT validation failed\n";
            return 1;
        }
        benchmark("No-Protection", [&](NativeVector* element) {
            ntt::ForwardTransformToBitReverseInPlace(rootTable, preconRootTable, element);
            return false;
        });
        benchmark("DMR", [&](NativeVector* element) {
            return ntt::ForwardTransformDMR(rootTable, preconRootTable, element);
        });
        benchmark("Naive-ABFT", [&](NativeVector* element) {
            return ntt::ForwardTransformABFTNaive(rootTable, preconRootTable, element, aOrig, preconAOrig);
        });
        benchmark("Ours", [&](NativeVector* element) {
            return ntt::ForwardTransformABFT(rootTable, preconRootTable, element, aBitReverse, level);
        });
    }
    else {
        NativeVector reference(sourcePolys.front());
        intt::InverseTransformFromBitReverseInPlace(inverseRootTable, preconInverseRootTable, nInv, preconNInv, &reference);
        NativeVector dmrResult(sourcePolys.front());
        const bool dmrMismatch = intt::InverseTransformDMR(inverseRootTable, preconInverseRootTable, nInv, preconNInv, &dmrResult);
        NativeVector naiveResult(sourcePolys.front());
        const bool naiveMismatch = intt::InverseTransformABFTNaive(inverseRootTable, preconInverseRootTable, nInv, preconNInv,
                                                                    &naiveResult, aOrig, preconAOrig);
        NativeVector oursResult(sourcePolys.front());
        const bool oursMismatch = intt::InverseTransformABFT(inverseRootTable, preconInverseRootTable, nInv, preconNInv,
                                                              &oursResult, aBitReverse, level);
        if (dmrMismatch || naiveMismatch || oursMismatch || !VectorsEqual(reference, dmrResult) ||
            !VectorsEqual(reference, naiveResult) || !VectorsEqual(reference, oursResult)) {
            std::cerr << "INTT validation failed\n";
            return 1;
        }
        benchmark("No-Protection", [&](NativeVector* element) {
            intt::InverseTransformFromBitReverseInPlace(inverseRootTable, preconInverseRootTable, nInv, preconNInv, element);
            return false;
        });
        benchmark("DMR", [&](NativeVector* element) {
            return intt::InverseTransformDMR(inverseRootTable, preconInverseRootTable, nInv, preconNInv, element);
        });
        benchmark("Naive-ABFT", [&](NativeVector* element) {
            return intt::InverseTransformABFTNaive(inverseRootTable, preconInverseRootTable, nInv, preconNInv, element,
                                                   aOrig, preconAOrig);
        });
        benchmark("Ours", [&](NativeVector* element) {
            return intt::InverseTransformABFT(inverseRootTable, preconInverseRootTable, nInv, preconNInv, element,
                                              aBitReverse, level);
        });
    }
    return 0;
}
