#include "binfhecontext.h"
#include "pbs-operator-breakdown.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace lbcrypto;

namespace {

constexpr size_t kStageCount = static_cast<size_t>(PBSOperatorBreakdownStage::COUNT);

// Accumulates the elapsed time of every profiling stage from one trial.
void Add(PBSOperatorBreakdownData& destination, const PBSOperatorBreakdownData& source) {
    for (size_t i = 0; i < kStageCount; ++i)
        destination.nanoseconds[i] += source.nanoseconds[i];
}

// Returns the elapsed time recorded for the requested profiling stage.
uint64_t Value(const PBSOperatorBreakdownData& data, PBSOperatorBreakdownStage stage) {
    return data.nanoseconds[static_cast<size_t>(stage)];
}

// f_8 = [1, 2, 3, 4, 7, 6, 5, 4] is negacyclic modulo 8:
// f_8(x + 4) = 8 - f_8(x). It therefore uses the one-PBS EvalFunc path.
NativeInteger NegacyclicLUTP8(NativeInteger m, NativeInteger p) {
    if (p != 8)
        throw std::invalid_argument("NegacyclicLUTP8 requires p=8");
    const auto value = m.ConvertToInt<uint32_t>();
    return value < 4 ? NativeInteger(value + 1) : NativeInteger(11 - value);
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

BINFHE_PARAMSET ParseP8Security(const std::string& security) {
    if (security == "STD128")
        return STD128_4;
    if (security == "STD192")
        return STD192_4;
    if (security == "STD256")
        return STD256_4;
    throw std::invalid_argument("security must be STD128, STD192, or STD256");
}

void PrintStage(const char* name, uint64_t nanoseconds, uint64_t total, uint32_t iterations) {
    const double meanMilliseconds = static_cast<double>(nanoseconds) / iterations / 1.0e6;
    const double share = total == 0 ? 0.0 : 100.0 * static_cast<double>(nanoseconds) / total;
    std::cout << std::left << std::setw(18) << name << std::right << std::setw(12) << std::fixed
              << std::setprecision(3) << meanMilliseconds << " ms  " << std::setw(7) << std::setprecision(2) << share
              << "%\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    const std::string security = argc > 1 ? argv[1] : "STD128";
    const uint32_t plaintextModulus = argc > 2 ? static_cast<uint32_t>(std::stoul(argv[2])) : 4;
    const uint32_t warmups = argc > 3 ? static_cast<uint32_t>(std::stoul(argv[3])) : 2;
    const uint32_t iterations = argc > 4 ? static_cast<uint32_t>(std::stoul(argv[4])) : 5;

    if (plaintextModulus != 2 && plaintextModulus != 4 && plaintextModulus != 8)
        throw std::invalid_argument("plaintext modulus must be 2, 4, or 8");
    if (iterations == 0)
        throw std::invalid_argument("iterations must be nonzero");

    BinFHEContext cc;
    if (plaintextModulus == 8) {
        // EvalFunc uses beta=128. The _4 parameter sets use q >= 2048,
        // which supports p=8 while retaining the requested security level.
        cc.GenerateBinFHEContext(ParseP8Security(security));
    }
    else {
        cc.GenerateBinFHEContext(ParseSecurity(security));
    }
    auto sk = cc.KeyGen();
    cc.BTKeyGen(sk);

    // The standard Bootstrap path is a Boolean-gate refresh. For p=2 and p=4
    // benchmark only Boolean messages {0,1}; p=8 uses the functional LUT path.
    const uint32_t inputCount = plaintextModulus == 8 ? plaintextModulus : 2;
    std::vector<LWECiphertext> inputs;
    inputs.reserve(inputCount);
    for (uint32_t m = 0; m < inputCount; ++m)
        inputs.push_back(cc.Encrypt(sk, m, SMALL_DIM, plaintextModulus));

    std::vector<NativeInteger> p8LUT;
    if (plaintextModulus == 8)
        p8LUT = cc.GenerateLUTviaFunction(NegacyclicLUTP8, plaintextModulus);

    const auto Evaluate = [&](ConstLWECiphertext& input) {
        return plaintextModulus == 8 ? cc.EvalFunc(input, p8LUT) : cc.Bootstrap(input);
    };
    const auto Expected = [&](uint32_t input) -> LWEPlaintext {
        return plaintextModulus == 8 ? NegacyclicLUTP8(input, plaintextModulus).ConvertToInt<LWEPlaintext>() : input;
    };

    PBSOperatorBreakdown::Enable(false);
    for (uint32_t i = 0; i < warmups; ++i)
        Evaluate(inputs[i % inputs.size()]);

    PBSOperatorBreakdownData aggregate;
    for (uint32_t i = 0; i < iterations; ++i) {
        PBSOperatorBreakdown::Reset();
        PBSOperatorBreakdown::Enable(true);
        auto output = Evaluate(inputs[i % inputs.size()]);
        PBSOperatorBreakdown::Enable(false);
        Add(aggregate, PBSOperatorBreakdown::Snapshot());

        LWEPlaintext decoded = 0;
        cc.Decrypt(sk, output, &decoded, plaintextModulus);
        const auto expected = Expected(i % inputs.size());
        if (decoded != expected)
            throw std::runtime_error("Bootstrap correctness check failed");
    }

    const uint64_t total = Value(aggregate, PBSOperatorBreakdownStage::TOTAL);
    const uint64_t ntt = Value(aggregate, PBSOperatorBreakdownStage::NTT);
    const uint64_t intt = Value(aggregate, PBSOperatorBreakdownStage::INTT);
    const uint64_t decompose = Value(aggregate, PBSOperatorBreakdownStage::DECOMPOSE);
    const uint64_t vcmm = Value(aggregate, PBSOperatorBreakdownStage::VCMM);
    const uint64_t keySwitch = Value(aggregate, PBSOperatorBreakdownStage::KEYSWITCH);
    const uint64_t fiveOperators = ntt + intt + decompose + vcmm + keySwitch;
    const uint64_t other = total > fiveOperators ? total - fiveOperators : 0;

    const std::string parameterSet = plaintextModulus == 8 ? security + "_4" : security;
    std::cout << "PBS operator breakdown: parameter set=" << parameterSet << ", p=" << plaintextModulus
              << (plaintextModulus == 8 ? ", negacyclic EvalFunc LUT" : ", Bootstrap") << ", warmups=" << warmups
              << ", iterations=" << iterations << "\n";
    PrintStage("NTT", ntt, total, iterations);
    PrintStage("INTT", intt, total, iterations);
    PrintStage("VCMM", vcmm, total, iterations);
    PrintStage("Decompose", decompose, total, iterations);
    PrintStage("KeySwitch", keySwitch, total, iterations);
    PrintStage("Other", other, total, iterations);
    PrintStage("Five operators", fiveOperators, total, iterations);
    PrintStage("PBS total", total, total, iterations);
    return 0;
}
