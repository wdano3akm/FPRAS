#include "../algorithm/algorithm.cpp"
#include "../cfg/CFGParser.hpp"

#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

void assertNear(double actual, double expected)
{
    const double tolerance = 1e-12 * std::max(1.0, std::abs(expected));
    assert(std::abs(actual - expected) <= tolerance);
}

double bruteRoundDown(std::size_t n, double kappa, int degree, double v)
{
    if (!(v > 0.0)) {
        return 0.0;
    }
    if (v > 1.0) {
        return 1.0;
    }

    const std::size_t maxEll = static_cast<std::size_t>(1) << n;
    const double factor = std::exp(kappa * degree);
    const double lowScale = 16.0 * n / factor;
    const double highScale = 16.0 * n * factor;
    double best = 0.0;

    for (std::size_t ell = 1; ell <= maxEll; ++ell) {
        const double low = lowScale / ell;
        const double high = highScale / ell;
        if (low < v) {
            best = std::max(best, low);
        }
        if (high < v) {
            best = std::max(best, high);
        }
    }

    return best;
}

PlusTimesProgram makeTwoVarAdd()
{
    PlusTimesProgram P;
    P.numVars = 2;
    P.degree = 1;
    P.nodes.push_back(PTNode{PTKind::Var, 1, 0, {}});
    P.nodes.push_back(PTNode{PTKind::Var, 1, 1, {}});
    P.nodes.push_back(PTNode{PTKind::Add, 1, 0, {0, 1}});
    P.root = 2;
    return P;
}

PlusTimesProgram makeTwoVarProduct()
{
    PlusTimesProgram P;
    P.numVars = 2;
    P.degree = 2;
    P.nodes.push_back(PTNode{PTKind::Var, 1, 0, {}});
    P.nodes.push_back(PTNode{PTKind::Var, 1, 1, {}});
    P.nodes.push_back(PTNode{PTKind::Mul, 2, 0, {0, 1}});
    P.root = 2;
    return P;
}

PlusTimesProgram makeNestedAdd()
{
    PlusTimesProgram P;
    P.numVars = 2;
    P.degree = 1;
    P.nodes.push_back(PTNode{PTKind::Var, 1, 0, {}});
    P.nodes.push_back(PTNode{PTKind::Var, 1, 1, {}});
    P.nodes.push_back(PTNode{PTKind::Add, 1, 0, {0, 1}});
    P.nodes.push_back(PTNode{PTKind::Add, 1, 0, {2, 0}});
    P.root = 3;
    return P;
}

PlusTimesProgram makeSharedSampleDag()
{
    PlusTimesProgram P;
    P.numVars = 4;
    P.degree = 2;
    P.nodes.push_back(PTNode{PTKind::Var, 1, 0, {}});
    P.nodes.push_back(PTNode{PTKind::Var, 1, 1, {}});
    P.nodes.push_back(PTNode{PTKind::Var, 1, 2, {}});
    P.nodes.push_back(PTNode{PTKind::Var, 1, 3, {}});
    P.nodes.push_back(PTNode{PTKind::Add, 1, 0, {0, 1}});
    P.nodes.push_back(PTNode{PTKind::Mul, 2, 0, {4, 2}});
    P.nodes.push_back(PTNode{PTKind::Mul, 2, 0, {4, 3}});
    P.nodes.push_back(PTNode{PTKind::Add, 2, 0, {5, 6}});
    P.root = 7;
    return P;
}

PlusTimesProgram makeLeftDeepProduct(int n)
{
    PlusTimesProgram P;
    P.numVars = static_cast<uint32_t>(n);
    P.degree = n;

    for (int i = 0; i < n; ++i) {
        P.nodes.push_back(PTNode{PTKind::Var, 1, static_cast<PolyVarId>(i), {}});
    }

    NodeId root = 0;
    int degree = 1;
    for (int i = 1; i < n; ++i) {
        ++degree;
        P.nodes.push_back(PTNode{
            PTKind::Mul,
            degree,
            0,
            {root, static_cast<NodeId>(i)}});
        root = static_cast<NodeId>(P.nodes.size() - 1);
    }

    P.root = root;
    return P;
}

NodeId appendBalancedProduct(PlusTimesProgram& P, const std::vector<NodeId>& factors, int begin, int end)
{
    if (end - begin == 1) {
        return factors[begin];
    }

    const int mid = begin + (end - begin) / 2;
    const NodeId left = appendBalancedProduct(P, factors, begin, mid);
    const NodeId right = appendBalancedProduct(P, factors, mid, end);
    P.nodes.push_back(PTNode{
        PTKind::Mul,
        P.nodes[left].degree + P.nodes[right].degree,
        0,
        {left, right}});
    return static_cast<NodeId>(P.nodes.size() - 1);
}

PlusTimesProgram makeBinaryChoiceProduct(int degree)
{
    PlusTimesProgram P;
    P.numVars = static_cast<uint32_t>(2 * degree);
    P.degree = degree;

    std::vector<NodeId> factors;
    factors.reserve(static_cast<std::size_t>(degree));

    for (int i = 0; i < degree; ++i) {
        const NodeId first = static_cast<NodeId>(P.nodes.size());
        P.nodes.push_back(PTNode{PTKind::Var, 1, static_cast<PolyVarId>(2 * i), {}});
        const NodeId second = static_cast<NodeId>(P.nodes.size());
        P.nodes.push_back(PTNode{PTKind::Var, 1, static_cast<PolyVarId>(2 * i + 1), {}});
        P.nodes.push_back(PTNode{PTKind::Add, 1, 0, {first, second}});
        factors.push_back(static_cast<NodeId>(P.nodes.size() - 1));
    }

    P.root = appendBalancedProduct(P, factors, 0, degree);
    return P;
}

PlusTimesProgram makeRandomizedAdditionProgram(int degree)
{
    constexpr int choicePositions = 9;
    assert(degree > choicePositions);

    PlusTimesProgram P;
    P.degree = degree;

    PolyVarId nextVariable = 0;
    const NodeId firstBranchVariable = static_cast<NodeId>(P.nodes.size());
    P.nodes.push_back(PTNode{PTKind::Var, 1, nextVariable++, {}});
    const NodeId secondBranchVariable = static_cast<NodeId>(P.nodes.size());
    P.nodes.push_back(PTNode{PTKind::Var, 1, nextVariable++, {}});

    std::vector<NodeId> tailFactors;
    tailFactors.reserve(static_cast<std::size_t>(degree - 1));
    for (int position = 0; position < choicePositions; ++position) {
        const NodeId first = static_cast<NodeId>(P.nodes.size());
        P.nodes.push_back(PTNode{
            PTKind::Var,
            1,
            nextVariable++,
            {}});
        const NodeId second = static_cast<NodeId>(P.nodes.size());
        P.nodes.push_back(PTNode{
            PTKind::Var,
            1,
            nextVariable++,
            {}});
        P.nodes.push_back(PTNode{PTKind::Add, 1, 0, {first, second}});
        tailFactors.push_back(static_cast<NodeId>(P.nodes.size() - 1));
    }
    for (int position = choicePositions; position < degree - 1; ++position) {
        P.nodes.push_back(PTNode{PTKind::Var, 1, nextVariable++, {}});
        tailFactors.push_back(static_cast<NodeId>(P.nodes.size() - 1));
    }
    P.numVars = static_cast<uint32_t>(nextVariable);

    const NodeId sharedTail = appendBalancedProduct(
        P,
        tailFactors,
        0,
        static_cast<int>(tailFactors.size()));
    P.nodes.push_back(PTNode{
        PTKind::Mul,
        degree,
        0,
        {firstBranchVariable, sharedTail}});
    const NodeId firstBranch = static_cast<NodeId>(P.nodes.size() - 1);
    P.nodes.push_back(PTNode{
        PTKind::Mul,
        degree,
        0,
        {secondBranchVariable, sharedTail}});
    const NodeId secondBranch = static_cast<NodeId>(P.nodes.size() - 1);
    P.nodes.push_back(PTNode{
        PTKind::Add,
        degree,
        0,
        {firstBranch, secondBranch}});
    P.root = static_cast<NodeId>(P.nodes.size() - 1);
    return P;
}

template <typename Fn>
void assertThrowsInvalidArgument(Fn fn)
{
    bool threw = false;
    try {
        fn();
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);
}

template <typename Fn>
void assertThrowsOverflow(Fn fn)
{
    bool threw = false;
    try {
        fn();
    } catch (const std::overflow_error&) {
        threw = true;
    }
    assert(threw);
}

bool sameMonomialSet(const MonomialSet& left, const MonomialSet& right)
{
    if (left.size() != right.size()) {
        return false;
    }
    for (const Monomial& monomial : left) {
        if (!right.contains(monomial)) {
            return false;
        }
    }
    return true;
}

std::vector<std::size_t> collectGeometricIndices(
    std::size_t count,
    double probability,
    std::uint64_t runSeed,
    const SampleKey& key,
    std::size_t* uniformDraws = nullptr)
{
    GeometricIndexGenerator generator(count, probability, runSeed, key);
    std::vector<std::size_t> indices;
    while (const std::optional<std::size_t> index = generator.next()) {
        indices.push_back(*index);
    }
    if (uniformDraws != nullptr) {
        *uniformDraws = generator.uniformDraws();
    }
    return indices;
}

void assertSortedUniqueAndInBounds(
    const std::vector<std::size_t>& indices,
    std::size_t count)
{
    for (std::size_t i = 0; i < indices.size(); ++i) {
        assert(indices[i] < count);
        if (i > 0) {
            assert(indices[i - 1] < indices[i]);
        }
    }
}

std::vector<Monomial> makeDegreeOneSupport(const std::vector<VarId>& variables)
{
    std::vector<Monomial> support;
    support.reserve(variables.size());
    for (VarId variable : variables) {
        support.push_back(Monomial{{variable}});
    }
    return support;
}

CountCoreState makeRegenerationState(const PlusTimesProgram& P, std::uint64_t seed)
{
    CountCoreState st;
    st.n = static_cast<std::size_t>(P.degree);
    st.kappa = 1.0;
    st.runSeed = seed;
    st.p.resize(P.nodes.size(), 0.0);
    st.pReady.resize(P.nodes.size(), false);
    st.exactSupports.resize(P.nodes.size());
    st.variableSets = computeVariableSets(P);
    st.membershipMemo.resize(P.nodes.size());
    return st;
}

void markExactSampleNode(
    CountCoreState& st,
    const PlusTimesProgram& P,
    NodeId q,
    double probability)
{
    Support support = enumerateSupportBounded(P, q, 100);
    assert(support.has_value());
    canonicalizeSupport(*support);
    st.exactSupports[q] = std::move(support);
    st.p[q] = probability;
    st.pReady[q] = true;
}

void testEnumerateSupportBounded()
{
    const PlusTimesProgram P = makeTwoVarAdd();

    Support complete = enumerateSupportBounded(P, P.root, 3);
    assert(complete.has_value());
    assert(complete->size() == 2);

    Support stoppedAtThreshold = enumerateSupportBounded(P, P.root, 2);
    assert(!stoppedAtThreshold.has_value());
}

void testGeometricIndexEdgeCases()
{
    const SampleKey key{3, 7, SampleStage::ExactSupport, 0};
    std::size_t draws = 99;

    assert(collectGeometricIndices(0, 0.5, 11, key, &draws).empty());
    assert(draws == 0);
    assert(collectGeometricIndices(10, 0.0, 11, key, &draws).empty());
    assert(draws == 0);
    assert(collectGeometricIndices(
        10,
        std::numeric_limits<double>::quiet_NaN(),
        11,
        key,
        &draws).empty());
    assert(draws == 0);
    assert(collectGeometricIndices(
        10,
        std::numeric_limits<double>::infinity(),
        11,
        key,
        &draws).empty());
    assert(draws == 0);

    const std::vector<std::size_t> all = collectGeometricIndices(5, 1.0, 11, key, &draws);
    assert((all == std::vector<std::size_t>{0, 1, 2, 3, 4}));
    assert(draws == 0);

    const std::vector<std::size_t> nearlyAll = collectGeometricIndices(
        100,
        std::nextafter(1.0, 0.0),
        11,
        key,
        &draws);
    assertSortedUniqueAndInBounds(nearlyAll, 100);
    assert(draws == nearlyAll.size() || draws == nearlyAll.size() + 1);

    const std::vector<std::size_t> hugeRange = collectGeometricIndices(
        std::numeric_limits<std::size_t>::max(),
        std::numeric_limits<double>::denorm_min(),
        11,
        key,
        &draws);
    assertSortedUniqueAndInBounds(
        hugeRange,
        std::numeric_limits<std::size_t>::max());
    assert(draws == hugeRange.size() || draws == hugeRange.size() + 1);
}

void testGeometricIndexDeterminism()
{
    const SampleKey key{5, 19, SampleStage::ExactSupport, 2};
    const std::vector<std::size_t> first =
        collectGeometricIndices(4096, 0.2, 1234567, key);
    const std::vector<std::size_t> intervening = collectGeometricIndices(
        4096,
        0.2,
        1234567,
        SampleKey{5, 20, SampleStage::ExactSupport, 2});
    const std::vector<std::size_t> repeated =
        collectGeometricIndices(4096, 0.2, 1234567, key);

    assert(first == repeated);
    assert(first != intervening);
    assert(first != collectGeometricIndices(4096, 0.2, 7654321, key));
    assertSortedUniqueAndInBounds(first, 4096);
}

void testSparseExactSupportCanonicalization()
{
    std::vector<Monomial> forward = makeDegreeOneSupport({1, 2, 3, 4, 5, 6, 7, 8});
    std::vector<Monomial> shuffled = makeDegreeOneSupport({6, 2, 8, 1, 7, 3, 5, 4});
    canonicalizeSupport(forward);
    canonicalizeSupport(shuffled);
    assert(forward == shuffled);

    CountCoreState st;
    st.runSeed = 99887766;
    const SampleKey key{4, 12, SampleStage::ExactSupport, 0};
    const MonomialSet first = reduceExactSupportSparse(forward, 0.4, st, key);
    const MonomialSet second = reduceExactSupportSparse(shuffled, 0.4, st, key);
    assert(sameMonomialSet(first, second));
}

void testGeometricSamplingDistributionAndSparseWork()
{
    constexpr std::size_t count = 32;
    constexpr std::size_t sampleCount = 4000;
    constexpr double probability = 0.25;
    std::vector<std::size_t> frequencies(count, 0);
    std::size_t total = 0;

    for (std::size_t sample = 0; sample < sampleCount; ++sample) {
        const std::vector<std::size_t> indices = collectGeometricIndices(
            count,
            probability,
            246813579,
            SampleKey{9, sample, SampleStage::ExactSupport, 0});
        assertSortedUniqueAndInBounds(indices, count);
        total += indices.size();
        for (std::size_t index : indices) {
            ++frequencies[index];
        }
    }

    assert(total > 30000);
    assert(total < 34000);
    for (std::size_t frequency : frequencies) {
        assert(frequency > 800);
        assert(frequency < 1200);
    }

    std::size_t draws = 0;
    const std::vector<std::size_t> sparse = collectGeometricIndices(
        1000000000,
        0.000001,
        135792468,
        SampleKey{10, 0, SampleStage::ExactSupport, 0},
        &draws);
    assertSortedUniqueAndInBounds(sparse, 1000000000);
    assert(sparse.size() < 5000);
    assert(draws == sparse.size() || draws == sparse.size() + 1);
}

void testRoundDown()
{
    CountCoreState st;
    st.n = 10;
    st.kappa = 0.01;
    PlusTimesProgram P;
    P.nodes = {
        PTNode{PTKind::Var, 2, 0, {}},
        PTNode{PTKind::Var, 4, 1, {}}};

    assert(roundDown(P, st, 0, 2.0) == 1.0);
    assert(roundDown(P, st, 0, 0.0) == 0.0);
    assert(roundDown(P, st, 0, -0.1) == 0.0);

    const double exactGridValue =
        (16.0 * st.n * std::exp(-st.kappa * P.nodes[0].degree)) / 200.0;
    const double rounded = roundDown(P, st, 0, exactGridValue);
    const double expected = bruteRoundDown(st.n, st.kappa, P.nodes[0].degree, exactGridValue);
    assert(rounded < exactGridValue);
    assertNear(rounded, expected);

    bool foundDegreeSensitiveValue = false;
    for (std::size_t ell = 100; ell < 400; ++ell) {
        const double v =
            (16.0 * st.n * std::exp(-st.kappa * P.nodes[0].degree)) /
            static_cast<double>(ell);
        const double d2 = roundDown(P, st, 0, v);
        const double d4 = roundDown(P, st, 1, v);
        if (std::abs(d2 - d4) > 1e-12) {
            assertNear(d2, bruteRoundDown(st.n, st.kappa, P.nodes[0].degree, v));
            assertNear(d4, bruteRoundDown(st.n, st.kappa, P.nodes[1].degree, v));
            foundDegreeSensitiveValue = true;
            break;
        }
    }
    assert(foundDegreeSensitiveValue);
}

void testCounterExactSmallPrograms()
{
    assert(counter(makeTwoVarProduct(), 1.0, 0.5) == 1.0);
}

void testPaperSupportThreshold()
{
    assert(supportThreshold(makeTwoVarAdd()) == 144);
    assert(supportThreshold(makeTwoVarProduct()) == 288);
    assert(supportThreshold(makeLeftDeepProduct(16)) == 246016);
    assertThrowsOverflow([] {
        supportThresholdFor(std::numeric_limits<std::size_t>::max(), 2);
    });
    assert(!theoremAllowsRandomizedCounting(15));
    assert(theoremAllowsRandomizedCounting(16));
}

void testAssertReadyForFPRAS()
{
    assertReadyForFPRAS(makeTwoVarAdd());
    assertReadyForFPRAS(makeTwoVarProduct());
    assertReadyForFPRAS(makeLeftDeepProduct(4));
    assertReadyForFPRAS(makeBinaryChoiceProduct(12));

    assertThrowsInvalidArgument([] {
        assertReadyForFPRAS(makeNestedAdd());
    });
    assertReadyForFPRAS(makeLeftDeepProduct(16));
    assert(counter(makeLeftDeepProduct(16), 1.0, 0.5) == 1.0);
}

void testDeterministicSampleRegeneration()
{
    {
        const PlusTimesProgram P = makeTwoVarProduct();
        CountCoreState st = makeRegenerationState(P, 123456789);
        markExactSampleNode(st, P, 0, 1.0);
        markExactSampleNode(st, P, 1, 1.0);
        st.p[2] = 0.5;
        st.pReady[2] = true;

        const MonomialSet first = generateSample(P, P.root, 3, st);
        const MonomialSet other = generateSample(P, P.root, 1, st);
        (void)other;
        const MonomialSet repeated = generateSample(P, P.root, 3, st);

        assert(sameMonomialSet(first, repeated));
    }

    {
        const PlusTimesProgram P = makeTwoVarAdd();
        CountCoreState st = makeRegenerationState(P, 987654321);
        markExactSampleNode(st, P, 0, 1.0);
        markExactSampleNode(st, P, 1, 1.0);
        st.p[2] = 0.5;
        st.pReady[2] = true;

        const MonomialSet first = generateSample(P, P.root, 4, st);
        const MonomialSet other = generateSample(P, P.root, 2, st);
        (void)other;
        const MonomialSet repeated = generateSample(P, P.root, 4, st);

        assert(sameMonomialSet(first, repeated));
    }
}

void testPerSampleMemoizationOnSharedDag()
{
    const PlusTimesProgram P = makeSharedSampleDag();
    CountCoreState st = makeRegenerationState(P, 11223344);
    for (NodeId q = 0; q <= 4; ++q) {
        markExactSampleNode(st, P, q, 1.0);
    }
    st.p[5] = 0.5;
    st.p[6] = 0.5;
    st.p[7] = 0.25;
    st.pReady[5] = true;
    st.pReady[6] = true;
    st.pReady[7] = true;

    SampleMemo memo(P.nodes.size());
    memo.beginSample(9);
    const MonomialSet& memoized = generateSampleMemoized(P, P.root, 9, st, memo);
    const MonomialSet regenerated = generateSample(P, P.root, 9, st);

    assert(memo.hitCount() > 0);
    assert(sameMonomialSet(memoized, regenerated));
}

void testExactThetaDecisions()
{
    const SampleSizeThreshold theta10(cpp_int(10));
    assert(exactThetaDecision(9, 1.0, theta10, 1) == ExactThetaDecision::Safe);
    assert(exactThetaDecision(10, 1.0, theta10, 1) == ExactThetaDecision::Fail);
    assert(exactThetaDecision(10, 0.5, theta10, 1) == ExactThetaDecision::CheckSamples);
    assert(exactThetaDecision(10, 0.0, theta10, 1) == ExactThetaDecision::Safe);
    assert(exactThetaDecision(10, 1.0, theta10, 0) == ExactThetaDecision::Safe);

    const SampleSizeThreshold theta0(cpp_int(0));
    assert(exactThetaDecision(0, 0.0, theta0, 1) == ExactThetaDecision::Fail);
}

void testAdditionThetaPassFolding()
{
    const PlusTimesProgram P = makeTwoVarAdd();

    CountCoreState highState = makeRegenerationState(P, 123);
    highState.n = 10;
    highState.kappa = 0.01;
    markExactSampleNode(highState, P, 0, 1.0);
    markExactSampleNode(highState, P, 1, 1.0);
    const PlusEstimateResult highResult = estimateSamplePlus(
        P, P.root, highState, 10, 1, 1, 1, false, SampleSizeThreshold(cpp_int(3)));
    assert(!highResult.needsFinalThetaPass);

    CountCoreState lowState = makeRegenerationState(P, 123);
    lowState.n = 10;
    lowState.kappa = 0.01;
    markExactSampleNode(lowState, P, 0, 1.0);
    markExactSampleNode(lowState, P, 1, 1.0);
    const PlusEstimateResult lowResult = estimateSamplePlus(
        P, P.root, lowState, 10, 1, 1, 1, false, SampleSizeThreshold(cpp_int(2)));
    assert(lowResult.needsFinalThetaPass);
}

void testSampleCountOverflowFailsFast()
{
    assertThrowsOverflow([] {
        countCoreWithSupportThreshold(
            makeTwoVarProduct(),
            std::numeric_limits<std::size_t>::max(),
            2,
            1,
            1.0,
            0);
    });
}

void testCountCoreSamplingSmoke()
{
    const int degree = 16;
    const double epsilon = 0.25;
    const double kappa = epsilon / (4.0 * degree);
    const std::size_t scaledNs = 8;
    const std::size_t scaledNt = 9;
    const std::size_t forcedSupportThreshold = 512;
    const PlusTimesProgram P = makeRandomizedAdditionProgram(degree);
    assert(P.nodes.size() == 52);
    assert(P.nodes[P.root].kind == PTKind::Add);
    assertReadyForFPRAS(P);

    const Support bruteForceSupport = enumerateSupportBounded(
        P,
        P.root,
        std::numeric_limits<std::size_t>::max());
    assert(bruteForceSupport.has_value());
    const double exact = static_cast<double>(bruteForceSupport->size());
    assert(exact == 1024.0);
    assert(exact <= static_cast<double>(std::size_t{1} << degree));
    assert(supportThreshold(P) > bruteForceSupport->size());

    for (NodeId child : P.nodes[P.root].children) {
        const Support childSupport = enumerateSupportBounded(
            P,
            child,
            forcedSupportThreshold + 1);
        assert(childSupport.has_value());
        assert(childSupport->size() == forcedSupportThreshold);
        assert(childSupport->size() > 16 * static_cast<std::size_t>(degree));
    }

    Support stoppedAtTestThreshold = enumerateSupportBounded(
        P,
        P.root,
        forcedSupportThreshold + 1);
    assert(!stoppedAtTestThreshold.has_value());

    const std::size_t productionNs = ceilToSize(
        48.0 / (kappa * kappa),
        "production ns overflows size_t");
    const std::size_t productionNt = checkedMultiply(
        checkedMultiply(8, static_cast<std::size_t>(degree), "production nt overflows size_t"),
        P.nodes.size(),
        "production nt overflows size_t");
    //assert(productionNs == 786432);
    //assert(productionNt == 6656);

    const cpp_int theta =
        cpp_int(512) * scaledNs * scaledNt * degree * P.nodes.size();

    std::vector<double> estimates;
    for (uint32_t seed : {12345u, 23456u, 34567u, 45678u, 56789u}) {
        seedSamplingForTesting(seed);
        estimates.push_back(countCoreWithSupportThreshold(
            P,
            scaledNs,
            scaledNt,
            theta,
            kappa,
            forcedSupportThreshold));
    }

    const double estimate = median(estimates);
    const double relativeError = std::abs(estimate - exact) / exact;
    std::cout << "randomized addition smoke test: runs=[";
    for (std::size_t i = 0; i < estimates.size(); ++i) {
        if (i > 0) {
            std::cout << ", ";
        }
        std::cout << estimates[i];
    }
    std::cout << "], median=" << estimate
              << ", brute-force exact=" << exact
              << ", n=" << degree
              << ", kappa=" << kappa
              << ", ns=" << scaledNs
              << ", nt=" << scaledNt
              << ", theta=" << theta
              << ", forced cutoff=" << forcedSupportThreshold
              << ", relative error=" << relativeError << std::endl;
    assert(estimate > 0.0);
    assert(estimate >= (1.0 - epsilon) * exact);
    assert(estimate <= (1.0 + epsilon) * exact);
    assert(std::adjacent_find(estimates.begin(), estimates.end(), std::not_equal_to<double>())
        != estimates.end());
}

void testRandomizedFixtures()
{
    {
        const CFG cfg = parseCFG("tests/fixtures/randomized_tiny.cfg");
        const PlusTimesProgram P = compileCFGToPlusTimes(cfg, 1);
        assert(P.nodes.size() == 18);
        assert(supportThreshold(P) == 5184);
        assert(counter(P, 0.5, 0.25) == 17.0);

        seedSamplingForTesting(4242);
        const double estimate = countCoreWithSupportThreshold(
            P,
            4,
            8,
            1000000,
            0.1,
            16);
        assert(std::isfinite(estimate));
        assert(estimate >= 0.0);
    }

    {
        const CFG cfg = parseCFG("tests/fixtures/randomized_theorem.cfg");
        const PlusTimesProgram P = compileCFGToPlusTimes(cfg, 22);
        const std::size_t exact = 4194304;
        assert(P.nodes.size() == 108);
        assert(supportThreshold(P) == 4105728);
        assert(exact > supportThreshold(P));
        assert(theoremAllowsRandomizedCounting(P.degree));

        std::vector<double> estimates;
        for (uint32_t seed : {101u, 202u, 303u}) {
            seedSamplingForTesting(seed);
            estimates.push_back(countCoreWithSupportThreshold(
                P,
                2,
                4,
                1000000,
                0.1,
                64));
        }
        const double estimate = median(estimates);
        assert(estimate > static_cast<double>(exact) / 100.0);
        assert(estimate < static_cast<double>(exact) * 100.0);
    }
}

} // namespace

int main()
{
    testEnumerateSupportBounded();
    testGeometricIndexEdgeCases();
    testGeometricIndexDeterminism();
    testSparseExactSupportCanonicalization();
    testGeometricSamplingDistributionAndSparseWork();
    testRoundDown();
    testCounterExactSmallPrograms();
    testPaperSupportThreshold();
    testAssertReadyForFPRAS();
    testDeterministicSampleRegeneration();
    testPerSampleMemoizationOnSharedDag();
    testExactThetaDecisions();
    testAdditionThetaPassFolding();
    testSampleCountOverflowFailsFast();
    testCountCoreSamplingSmoke();
    testRandomizedFixtures();
    return 0;
}
