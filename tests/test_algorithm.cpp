#include "../algorithm/algorithm.cpp"

#include <cassert>
#include <cmath>
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
    assert(counter(makeTwoVarAdd(), 1.0, 0.5) == 2.0);
    assert(counter(makeTwoVarProduct(), 1.0, 0.5) == 1.0);
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
    const int degree = 12;
    const PlusTimesProgram P = makeBinaryChoiceProduct(degree);
    const double exact = static_cast<double>(std::size_t{1} << degree);

    Support stoppedAtTestThreshold = enumerateSupportBounded(P, P.root, 65);
    assert(!stoppedAtTestThreshold.has_value());

    std::vector<double> estimates;
    for (uint32_t seed : {12345u, 23456u, 34567u, 45678u, 56789u}) {
        seedSamplingForTesting(seed);
        estimates.push_back(countCoreWithSupportThreshold(
            P,
            4,
            8,
            1000000,
            0.1,
            64));
    }

    const double estimate = median(estimates);
    assert(estimate > 0.0);
    assert(estimate > exact / 4.0);
    assert(estimate < exact * 4.0);
}

} // namespace

int main()
{
    testEnumerateSupportBounded();
    testRoundDown();
    testCounterExactSmallPrograms();
    testAssertReadyForFPRAS();
    testDeterministicSampleRegeneration();
    testSampleCountOverflowFailsFast();
    testCountCoreSamplingSmoke();
    return 0;
}
