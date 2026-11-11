#include "binfhecontext.h"
#include "utils/pbs-sbf-fault-injection.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace lbcrypto;

namespace {
BINFHE_PARAMSET ParseSecurity(const std::string& security) {
    if (security == "STD128") return STD128;
    if (security == "STD192") return STD192;
    if (security == "STD256") return STD256;
    throw std::invalid_argument("security must be STD128, STD192, or STD256");
}

BINFHE_PARAMSET ParseP8Security(const std::string& security) {
    if (security == "STD128") return STD128_4;
    if (security == "STD192") return STD192_4;
    if (security == "STD256") return STD256_4;
    throw std::invalid_argument("security must be STD128, STD192, or STD256");
}

PBSFaultSite Site(const std::string& op, const std::string& value) {
    if (op == "ntt" && value == "NTT-butterfly-load") return PBSFaultSite::BUTTERFLY_LOAD;
    if (op == "ntt" && value == "NTT-butterfly-add") return PBSFaultSite::BUTTERFLY_ADD;
    if (op == "ntt" && value == "NTT-butterfly-sub") return PBSFaultSite::BUTTERFLY_SUB;
    if (op == "ntt" && value == "NTT-butterfly-mul") return PBSFaultSite::BUTTERFLY_MUL;
    if (op == "ntt" && value == "NTT-butterfly-store") return PBSFaultSite::BUTTERFLY_STORE;
    if (op == "intt" && value == "INTT-butterfly-load") return PBSFaultSite::BUTTERFLY_LOAD;
    if (op == "intt" && value == "INTT-butterfly-add") return PBSFaultSite::BUTTERFLY_ADD;
    if (op == "intt" && value == "INTT-butterfly-sub") return PBSFaultSite::BUTTERFLY_SUB;
    if (op == "intt" && value == "INTT-butterfly-mul") return PBSFaultSite::BUTTERFLY_MUL;
    if (op == "intt" && value == "INTT-butterfly-store") return PBSFaultSite::BUTTERFLY_STORE;
    if (op == "vcmm" && value == "input") return PBSFaultSite::VCMM_INPUT;
    if (op == "vcmm" && value == "multiplication") return PBSFaultSite::VCMM_MULTIPLICATION;
    if (op == "vcmm" && value == "accumulation") return PBSFaultSite::VCMM_ACCUMULATION;
    if (op == "vcmm" && value == "output") return PBSFaultSite::VCMM_OUTPUT;
    if (op == "decompose" && value == "digit-extraction") return PBSFaultSite::DECOMPOSE_DIGIT_EXTRACTION;
    if (op == "decompose" && value == "intermediate-digit") return PBSFaultSite::DECOMPOSE_INTERMEDIATE_DIGIT;
    if (op == "keyswitch" && value == "key-lookup") return PBSFaultSite::KEYSWITCH_KEY_LOOKUP;
    if (op == "keyswitch" && value == "accumulation") return PBSFaultSite::KEYSWITCH_ACCUMULATION;
    if (op == "keyswitch" && value == "update") return PBSFaultSite::KEYSWITCH_UPDATE;
    throw std::invalid_argument("operator/site combination is not supported");
}

using NamedSite = std::pair<std::string, PBSFaultSite>;

struct BlockLayout {
    uint32_t evalAccCount;
    uint32_t digitsG2;
    uint32_t ringDimension;
    uint32_t bitWidth;
};

BlockLayout GetBlockLayout(const std::string& security, uint32_t p) {
    const bool p8 = p == 8;
    if (security == "STD128") return p8 ? BlockLayout{595, 6, 1024, 27} : BlockLayout{503, 4, 1024, 27};
    if (security == "STD192") return p8 ? BlockLayout{875, 4, 2048, 37} : BlockLayout{790, 2, 2048, 37};
    if (security == "STD256") return p8 ? BlockLayout{1145, 6, 2048, 29} : BlockLayout{1076, 4, 2048, 29};
    throw std::invalid_argument("security must be STD128, STD192, or STD256");
}

uint64_t ButterflyCount(uint32_t n) {
    uint32_t logn = 0;
    for (uint32_t value = n; value > 1; value >>= 1)
        ++logn;
    return static_cast<uint64_t>(n / 2) * logn;
}

uint64_t BlockCandidateCount(const std::string& op, PBSFaultSite site, const BlockLayout& layout) {
    if (op == "ntt" || op == "intt") {
        const auto butterflies = ButterflyCount(layout.ringDimension);
        return (site == PBSFaultSite::BUTTERFLY_LOAD || site == PBSFaultSite::BUTTERFLY_STORE) ?
                   2 * butterflies : butterflies;
    }
    if (op == "decompose")
        return static_cast<uint64_t>(layout.digitsG2) * layout.ringDimension;
    if (op == "vcmm") {
        if (site == PBSFaultSite::VCMM_INPUT)
            return static_cast<uint64_t>(layout.digitsG2) * layout.ringDimension;
        if (site == PBSFaultSite::VCMM_ACCUMULATION)
            return 4ULL * (layout.digitsG2 - 1) * layout.ringDimension;
        return 4ULL * layout.ringDimension;
    }
    return 0;
}

PBSFaultBlockKind BlockKind(const std::string& op) {
    if (op == "ntt") return PBSFaultBlockKind::NTT;
    if (op == "intt") return PBSFaultBlockKind::INTT;
    if (op == "decompose") return PBSFaultBlockKind::DECOMPOSE;
    if (op == "vcmm") return PBSFaultBlockKind::VCMM;
    return PBSFaultBlockKind::NONE;
}

std::vector<NamedSite> Sites(const std::string& op) {
    if (op == "ntt")
        return {{"NTT-butterfly-load", PBSFaultSite::BUTTERFLY_LOAD},
                {"NTT-butterfly-add", PBSFaultSite::BUTTERFLY_ADD},
                {"NTT-butterfly-sub", PBSFaultSite::BUTTERFLY_SUB},
                {"NTT-butterfly-mul", PBSFaultSite::BUTTERFLY_MUL},
                {"NTT-butterfly-store", PBSFaultSite::BUTTERFLY_STORE}};
    if (op == "intt")
        return {{"INTT-butterfly-load", PBSFaultSite::BUTTERFLY_LOAD},
                {"INTT-butterfly-add", PBSFaultSite::BUTTERFLY_ADD},
                {"INTT-butterfly-sub", PBSFaultSite::BUTTERFLY_SUB},
                {"INTT-butterfly-mul", PBSFaultSite::BUTTERFLY_MUL},
                {"INTT-butterfly-store", PBSFaultSite::BUTTERFLY_STORE}};
    if (op == "vcmm")
        return {{"input", PBSFaultSite::VCMM_INPUT},
                {"multiplication", PBSFaultSite::VCMM_MULTIPLICATION},
                {"accumulation", PBSFaultSite::VCMM_ACCUMULATION},
                {"output", PBSFaultSite::VCMM_OUTPUT}};
    if (op == "decompose")
        return {{"digit-extraction", PBSFaultSite::DECOMPOSE_DIGIT_EXTRACTION},
                {"intermediate-digit", PBSFaultSite::DECOMPOSE_INTERMEDIATE_DIGIT}};
    if (op == "keyswitch")
        return {{"key-lookup", PBSFaultSite::KEYSWITCH_KEY_LOOKUP},
                {"accumulation", PBSFaultSite::KEYSWITCH_ACCUMULATION},
                {"update", PBSFaultSite::KEYSWITCH_UPDATE}};
    throw std::invalid_argument("operator must be ntt, intt, vcmm, decompose, or keyswitch");
}

// Keep this LUT identical to evaluation-pbs-operator-breakdown.cpp.  It is
// negacyclic modulo 8 and therefore follows the one-PBS EvalFunc path.
NativeInteger NegacyclicLUTP8(NativeInteger m, NativeInteger p) {
    if (p != 8)
        throw std::invalid_argument("NegacyclicLUTP8 requires p=8");
    const auto value = m.ConvertToInt<uint32_t>();
    return value < 4 ? NativeInteger(value + 1) : NativeInteger(11 - value);
}

const char* SiteName(PBSFaultSite site) {
    static constexpr const char* names[] = {
        "butterfly-load", "butterfly-add", "butterfly-sub", "butterfly-mul", "butterfly-store",
        "VCMM-input", "VCMM-multiplication", "VCMM-accumulation", "VCMM-output",
        "Decompose-digit-extraction", "Decompose-intermediate-digit",
        "KeySwitch-key-lookup", "KeySwitch-accumulation", "KeySwitch-update"};
    return names[static_cast<size_t>(site)];
}

const char* TransformName(PBSFaultTransform transform) {
    switch (transform) {
        case PBSFaultTransform::ANY:
            return "other";
        case PBSFaultTransform::NTT:
            return "NTT";
        case PBSFaultTransform::INTT:
            return "INTT";
    }
    return "unknown";
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc == 4 && std::string(argv[1]) == "pre-explore") {
            const std::string security(argv[2]);
            const uint32_t p = std::stoul(argv[3]);
            if (p != 2 && p != 4 && p != 8)
                throw std::invalid_argument("p must be 2/4/8");

            BinFHEContext cc;
            if (p == 8)
                cc.GenerateBinFHEContext(ParseP8Security(security));
            else
                cc.GenerateBinFHEContext(ParseSecurity(security));
            auto sk = cc.KeyGen();
            cc.BTKeyGen(sk);
            std::vector<NativeInteger> p8LUT;
            if (p == 8)
                p8LUT = cc.GenerateLUTviaFunction(NegacyclicLUTP8, p);

            PBSFaultInjection::Disable();
            PBSFaultInjection::BeginPreExplore();
            auto probe = cc.Encrypt(sk, 0, SMALL_DIM, p);
            if (p == 8)
                cc.EvalFunc(probe, p8LUT);
            else
                cc.Bootstrap(probe);
            PBSFaultInjection::EndPreExplore();

            std::cout << "Site,Transform,SiteOffset,CandidateCount,BitWidth,Dim0,Dim1,Dim2,Dim3\n";
            for (const auto& block : PBSFaultInjection::GetPreExploreBlocks()) {
                std::cout << SiteName(block.site) << ',' << TransformName(block.transform) << ','
                          << block.siteOffset << ',' << block.candidateCount << ',' << block.bitWidth;
                for (const auto dimension : block.dimensions)
                    std::cout << ',' << dimension;
                std::cout << '\n';
            }
            return EXIT_SUCCESS;
        }
        if (argc < 5 || argc > 8) {
            std::cerr << "usage: " << argv[0]
                      << " <ntt|intt|vcmm|decompose|keyswitch> <site|all> <STD128|STD192|STD256> <p=2|4|8>"
                         " [iterations=10] [faults=1] [seed=7]\n"
                      << "   or: " << argv[0] << " pre-explore <STD128|STD192|STD256> <p=2|4|8>\n";
            return EXIT_FAILURE;
        }
        const std::string op(argv[1]), siteName(argv[2]), security(argv[3]);
        const uint32_t p = std::stoul(argv[4]);
        const uint64_t iterations = argc > 5 ? std::stoull(argv[5]) : 10;
        const uint32_t faults = argc > 6 ? std::stoul(argv[6]) : 1;
        const uint64_t seed = argc > 7 ? std::stoull(argv[7]) : 7;
        if ((p != 2 && p != 4 && p != 8) || iterations == 0 || (faults != 1 && faults != 2))
            throw std::invalid_argument("p must be 2/4/8, iterations nonzero, and faults 1/2");
        const auto selectedSites = siteName == "all" ? Sites(op) : std::vector<NamedSite>{{siteName, Site(op, siteName)}};
        const auto transform = op == "ntt" ? PBSFaultTransform::NTT :
                               op == "intt" ? PBSFaultTransform::INTT : PBSFaultTransform::ANY;

        BinFHEContext cc;
        // Match the breakdown evaluation exactly: p=8 needs the _4 parameter
        // sets (q >= 2048) and uses EvalFunc; p=2/4 uses Bootstrap.
        if (p == 8)
            cc.GenerateBinFHEContext(ParseP8Security(security));
        else
            cc.GenerateBinFHEContext(ParseSecurity(security));
        auto sk = cc.KeyGen();
        cc.BTKeyGen(sk);
        std::vector<NativeInteger> p8LUT;
        if (p == 8)
            p8LUT = cc.GenerateLUTviaFunction(NegacyclicLUTP8, p);
        const auto Evaluate = [&](ConstLWECiphertext& input) {
            return p == 8 ? cc.EvalFunc(input, p8LUT) : cc.Bootstrap(input);
        };
        const auto Expected = [&](uint32_t input) -> LWEPlaintext {
            return p == 8 ? NegacyclicLUTP8(input, p).ConvertToInt<LWEPlaintext>() : input;
        };

        std::mt19937_64 rng(seed);
        const auto blockLayout = GetBlockLayout(security, p);
        const uint64_t totalIterations = selectedSites.size() * iterations;
        uint64_t completedIterations = 0;
        for (const auto& [selectedName, faultSite] : selectedSites) {
            uint64_t candidates = 0;
            uint32_t bitWidth = 0;
            if (op == "keyswitch") {
                // KeySwitch remains a single whole-operator block.
                PBSFaultInjection::Configure(faultSite, {{UINT64_MAX, 0}}, transform);
                auto probe = cc.Encrypt(sk, 0, SMALL_DIM, p);
                Evaluate(probe);
                candidates = PBSFaultInjection::GetSeenCount();
                bitWidth = PBSFaultInjection::GetObservedBitWidth();
                PBSFaultInjection::Disable();
            }
            else {
                candidates = BlockCandidateCount(op, faultSite, blockLayout);
                bitWidth = blockLayout.bitWidth;
            }
            if (candidates < faults || bitWidth == 0)
                throw std::runtime_error("not enough dynamic fault candidates or no valid bit width");

            uint64_t errors = 0;
            for (uint64_t iteration = 0; iteration < iterations; ++iteration) {
                std::uniform_int_distribution<uint64_t> occurrence(0, candidates - 1);
                std::set<uint64_t> ids;
                while (ids.size() < faults)
                    ids.insert(occurrence(rng));
                std::vector<PBSFaultTarget> targets;
                std::uniform_int_distribution<uint32_t> bit(0, bitWidth - 1);
                for (const auto id : ids)
                    targets.push_back({id, bit(rng)});
                if (op == "keyswitch") {
                    PBSFaultInjection::Configure(faultSite, std::move(targets), transform);
                }
                else {
                    std::uniform_int_distribution<uint32_t> evalAcc(0, blockLayout.evalAccCount - 1);
                    uint32_t subBlock = 0;
                    if (op == "ntt") {
                        std::uniform_int_distribution<uint32_t> digit(0, blockLayout.digitsG2 - 1);
                        subBlock = digit(rng);
                    }
                    else if (op == "intt") {
                        std::uniform_int_distribution<uint32_t> ciphertext(0, 1);
                        subBlock = ciphertext(rng);
                    }
                    PBSFaultInjection::ConfigureBlock(faultSite, std::move(targets), transform, BlockKind(op),
                                                       evalAcc(rng), subBlock);
                }
                // The breakdown evaluation uses Boolean messages for p=2/4 and
                // the full message domain only for the p=8 functional LUT path.
                const uint32_t input = p == 8 ? rng() % p : rng() % 2;
                auto output = Evaluate(cc.Encrypt(sk, input, SMALL_DIM, p));
                if (PBSFaultInjection::GetInjectedCount() != faults)
                    throw std::runtime_error("requested faults were not all injected");
                PBSFaultInjection::Disable();
                LWEPlaintext decrypted = 0;
                cc.Decrypt(sk, output, &decrypted, p);
                const auto expected = Expected(input);
                errors += decrypted != expected;
                ++completedIterations;
                // stderr is deliberately used and flushed so that a trial log
                // can be monitored live with tail -f while stdout remains the
                // machine-readable result stream.
                std::cerr << "progress=" << completedIterations << '/' << totalIterations
                          << " site=" << selectedName << " site_iteration=" << iteration + 1 << '/' << iterations
                          << std::endl;
            }
            std::cout << "site=" << selectedName << " iterations=" << iterations << " errors=" << errors
                      << " error_rate=" << static_cast<double>(errors) / iterations << '\n';
        }
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }
}
