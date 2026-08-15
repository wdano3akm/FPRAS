#include "algorithm.hpp"
#include <algorithm>
#include <boost/multiprecision/fwd.hpp>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>
#include <optional>
#include <unordered_set>
#include <boost/container_hash/hash.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>

using boost::multiprecision::cpp_int;

constexpr std::size_t supportThresholdFactor = 16;
constexpr std::size_t minimumRandomizedDegree = 16;

bool theoremAllowsRandomizedCounting(std::size_t degree)
{
    return degree >= minimumRandomizedDegree;
}

// A monomial is a vector of variables 
// whose order does not matter
// Every monomial is assumed to be homogeneous. Therefore
// Degree(Monomial) == length(CFG_word)
struct Monomial {
    std::vector<VarId> vars;

    bool operator==(const Monomial& other) const {
        return vars == other.vars;
    }
};

bool monomialLess(const Monomial& left, const Monomial& right)
{
    return left.vars < right.vars;
}

void canonicalizeSupport(std::vector<Monomial>& support)
{
    std::sort(support.begin(), support.end(), monomialLess);
}

// Compute an hash for every Monomial using
// the variables contained inside it
struct MonomialHash {
    std::size_t operator()(const Monomial& m) const {
        return boost::hash_range(m.vars.begin(), m.vars.end());
    }
};

class MonomialSet {
public:
    bool insert(const Monomial& monomial) {
        return data_.insert(monomial).second;
    }

    bool contains(const Monomial& monomial) const {
        return data_.find(monomial) != data_.end();
    }

    std::size_t size() const {
        return data_.size();
    }

    auto begin() const {
        return data_.begin();
    }

    auto end() const {
        return data_.end();
    }

private:
    std::unordered_set<Monomial, MonomialHash> data_;
};

// if has_value(), then base case. Otherwise, assume 
// we enumerated Q^0 entirely and there is still more to go
using Support = std::optional<std::vector<Monomial>>; 
using SupportMembershipMemo = std::vector<std::unordered_map<Monomial, bool, MonomialHash>>;

// Struct used to store the state of the current
// counter iteration.
struct CountCoreState {
    std::size_t n = 0;
    double kappa = 0.0;
    std::uint64_t runSeed = 0;
    std::vector<double> p;
    std::vector<char> pReady;
    std::vector<Support> exactSupports;
    VariableSets variableSets;
    SupportMembershipMemo membershipMemo;
};

class SampleMemo {
public:
    explicit SampleMemo(std::size_t nodeCount)
        : entries_(nodeCount)
    {
        touched_.reserve(nodeCount);
    }

    void beginSample(std::size_t sample)
    {
        for (NodeId node : touched_) {
            entries_[node].reset();
        }
        touched_.clear();
        sample_ = sample;
        active_ = true;
    }

    const MonomialSet* find(NodeId node, std::size_t sample)
    {
        requireCurrentSample(sample);
        if (!entries_[node].has_value()) {
            return nullptr;
        }
        ++hitCount_;
        return &*entries_[node];
    }

    const MonomialSet& store(NodeId node, std::size_t sample, MonomialSet value)
    {
        requireCurrentSample(sample);
        if (entries_[node].has_value()) {
            throw std::logic_error("sample memo entry stored twice");
        }
        entries_[node].emplace(std::move(value));
        touched_.push_back(node);
        return *entries_[node];
    }

    std::size_t hitCount() const
    {
        return hitCount_;
    }

private:
    void requireCurrentSample(std::size_t sample) const
    {
        if (!active_ || sample != sample_) {
            throw std::logic_error("sample memo used with the wrong sample index");
        }
    }

    std::vector<std::optional<MonomialSet>> entries_;
    std::vector<NodeId> touched_;
    std::size_t sample_ = 0;
    std::size_t hitCount_ = 0;
    bool active_ = false;
};

class SampleSizeThreshold {
public:
    explicit SampleSizeThreshold(const cpp_int& theta)
    {
        if (theta <= 0) {
            zeroAlwaysReaches_ = true;
            reachable_ = true;
            value_ = 0;
            return;
        }

        const cpp_int maxSize = std::numeric_limits<std::size_t>::max();
        if (theta > maxSize) {
            return;
        }
        reachable_ = true;
        value_ = theta.convert_to<std::size_t>();
    }

    bool reached(std::size_t size) const
    {
        return reachable_ && size >= value_;
    }

    bool zeroAlwaysReaches() const
    {
        return zeroAlwaysReaches_;
    }

private:
    std::size_t value_ = 0;
    bool reachable_ = false;
    bool zeroAlwaysReaches_ = false;
};

enum class ExactThetaDecision {
    Safe,
    Fail,
    CheckSamples,
};

ExactThetaDecision exactThetaDecision(
    std::size_t supportSize,
    double probability,
    const SampleSizeThreshold& theta,
    std::size_t sampleCount)
{
    if (sampleCount == 0) {
        return ExactThetaDecision::Safe;
    }
    if (probability <= 0.0) {
        return theta.zeroAlwaysReaches()
            ? ExactThetaDecision::Fail
            : ExactThetaDecision::Safe;
    }
    if (!theta.reached(supportSize)) {
        return ExactThetaDecision::Safe;
    }
    if (probability >= 1.0) {
        return ExactThetaDecision::Fail;
    }
    return ExactThetaDecision::CheckSamples;
}

// Description:
//  Checked power of two function
// Input:
//  size_t exponent 
// Output:
//  double 2**exponent
long double powerOfTwo(std::size_t exponent)
{
    if (exponent > static_cast<std::size_t>(std::numeric_limits<long double>::max_exponent)) {
        return std::numeric_limits<long double>::infinity();
    }
    return std::ldexp(1.0L, static_cast<int>(exponent));
}

double acceptableValueBelow(double scale, double v, long double maxEll)
{
    if (!(scale > 0.0) || !std::isfinite(scale)) {
        return 0.0;
    }

    const long double scaleLong = scale;
    const long double vLong = v;
    long double ell = std::floor(scaleLong / vLong) + 1.0L;
    if (ell < 1.0L) {
        ell = 1.0L;
    }
    if (ell > maxEll) {
        return 0.0;
    }

    long double candidate = scaleLong / ell;
    double roundedCandidate = static_cast<double>(candidate);
    while (!(roundedCandidate < v) && ell < maxEll) {
        ell += 1.0L;
        candidate = scaleLong / ell;
        roundedCandidate = static_cast<double>(candidate);
    }

    if (candidate < vLong && roundedCandidate < v && candidate > 0.0L) {
        return roundedCandidate;
    }
    return 0.0;
}

double roundDown(const PlusTimesProgram& P, const CountCoreState& st, NodeId q, double v) {
    if (!(v > 0.0)) {
        return 0.0;
    }
    if (v > 1.0) {
        return 1.0;
    }
    if (q >= P.nodes.size() || P.nodes[q].degree <= 0) {
        throw std::logic_error("cannot round probability for an invalid node degree");
    }

    const long double maxEll = powerOfTwo(st.n);
    const long double exponent =
        static_cast<long double>(st.kappa) * P.nodes[q].degree;
    const long double base = 16.0L * static_cast<long double>(st.n);

    double best = 0.0;
    const double lowScale = static_cast<double>(base * std::exp(-exponent));
    const double highScale = static_cast<double>(base * std::exp(exponent));
    best = std::max(best, acceptableValueBelow(lowScale, v, maxEll));
    best = std::max(best, acceptableValueBelow(highScale, v, maxEll));
    return best;
}

struct BoundedSupportResult {
    bool computed = false;
    bool complete = true;
    MonomialSet support;
};

// Description:
//  Adds a monomial to the support.
//  Checks whether the size of the support is 
//  below the threshold
// Input:
//  MonomialSet& support 
//  Monomial &monomial 
//  size_t threshold
// Output
//  bool support.size() < threshold
bool insertBounded(MonomialSet& support, const Monomial& monomial, std::size_t threshold)
{
    support.insert(monomial);
    return support.size() < threshold;
}

// Description:
//  Recursively iterates through the tree from the 
//  Var nodes up. Returns a struct with the set of monomials
//  counted and a bool on whether the enumeration from the node
//  is completed
// Input:
//  PlusTimesProgram &P
//  NodeId q
//  size_t threshold
//  vector<BoundedSupportResult>& memo
// Output:
//  BoundedSupportResult struct with complete, computed 
//  and memoization
BoundedSupportResult enumerateSupportBoundedNode(
    const PlusTimesProgram& P,
    NodeId q,
    std::size_t threshold,
    std::vector<BoundedSupportResult>& memo)
{
    if (q == EMPTY_NODE) {
        BoundedSupportResult empty;
        empty.computed = true;
        return empty;
    }
    if (q >= P.nodes.size()) {
        throw std::logic_error("invalid plus-times node id");
    }

    BoundedSupportResult& cached = memo[q];
    if (cached.computed) {
        return cached;
    }

    BoundedSupportResult result;
    result.computed = true;
    const PTNode& node = P.nodes[q];

    if (node.kind == PTKind::Var) {
        Monomial monomial;
        monomial.vars.push_back(node.var);
        result.complete = insertBounded(result.support, monomial, threshold);
    } else if (node.kind == PTKind::Add) {
        for (NodeId child : node.children) {
            BoundedSupportResult childResult =
                enumerateSupportBoundedNode(P, child, threshold, memo);
            if (!childResult.complete) {
                result.complete = false;
                break;
            }
            for (const Monomial& monomial : childResult.support) {
                if (!insertBounded(result.support, monomial, threshold)) {
                    result.complete = false;
                    break;
                }
            }
            if (!result.complete) {
                break;
            }
        }
    } else {
        assert(node.children.size() == 2);
        BoundedSupportResult left =
            enumerateSupportBoundedNode(P, node.children[0], threshold, memo);
        BoundedSupportResult right =
            enumerateSupportBoundedNode(P, node.children[1], threshold, memo);

        if ((!left.complete && right.complete && right.support.size() == 0) ||
            (!right.complete && left.complete && left.support.size() == 0)) {
            result.complete = true;
        } else if (!left.complete || !right.complete) {
            result.complete = false;
        } else {
            for (const Monomial& leftMonomial : left.support) {
                for (const Monomial& rightMonomial : right.support) {
                    Monomial product = leftMonomial;
                    product.vars.insert(product.vars.end(), rightMonomial.vars.begin(), rightMonomial.vars.end());
                    std::sort(product.vars.begin(), product.vars.end());
                    if (!insertBounded(result.support, product, threshold)) {
                        result.complete = false;
                        break;
                    }
                }
                if (!result.complete) {
                    break;
                }
            }
        }
    }

    cached = result;
    return cached;
}

// Description:
//  Wrapper over enumerateSupportBoundedNode returning
//  a vector with every monomial inside the support while 
//  under the threshold
// Input:
//  PlusTimesProgram& P
//  NodeId o 
//  size_t threshold
// Output: 
//  Support, optional vector. if has_value() then
//  base case
Support enumerateSupportBounded(const PlusTimesProgram& P, NodeId o, std::size_t threshold)
{
    if (threshold == 0) {
        return std::nullopt;
    }

    std::vector<BoundedSupportResult> memo(P.nodes.size());
    BoundedSupportResult result = enumerateSupportBoundedNode(P, o, threshold, memo);
    if (!result.complete) {
        return std::nullopt;
    }

    return std::vector<Monomial>(result.support.begin(), result.support.end());
}

// Description:
//  Computes the support size until the threshold is hit.
// Input:
//  PlusTimesProgram &P
//  NodeId o 
//  size_t threshold
// Output:
//  optional<size_t> if has_value() base case and size, otherwise treshold
std::optional<std::size_t> supportSizeBounded(
    const PlusTimesProgram& P,
    NodeId o,
    std::size_t threshold)
{
    Support support = enumerateSupportBounded(P, o, threshold);
    if (!support.has_value()) {
        return std::nullopt;
    }
    return support->size();
}

double clampProbability(double p)
{
    if (p <= 0.0) {
        return 0.0;
    }
    if (p >= 1.0) {
        return 1.0;
    }
    return p;
}
// returns a mersenne twister 
std::mt19937& samplingGenerator()
{
    static thread_local std::mt19937 generator(std::random_device{}());
    return generator;
}

// helper function with set seed
void seedSamplingForTesting(uint32_t seed)
{
    samplingGenerator().seed(seed);
}

// returns true/false with probability p
bool keepWithProbability(double p)
{
    std::bernoulli_distribution keep(clampProbability(p));
    return keep(samplingGenerator());
}

// Description:
//  Keeps Monomials from the input in the returned MonomialSet 
//  with probability p
// Input:
//  vector<Monomial>& monomials
//  double p
// Output: 
//  MonomialSet returned set
MonomialSet reduce(const std::vector<Monomial>& monomials, double p){
    MonomialSet result;
    for (const Monomial& m : monomials) {
        if (keepWithProbability(p)) {
            result.insert(m);
        }
    }
    return result;
}

// Description:
//  Keeps Monomials from the input in the returned MonomialSet 
//  with probability p
// Input:
//  MonomialSet& monomials
//  double p
// Output: 
//  MonomialSet returned set
MonomialSet reduce(const MonomialSet& monomials, double p)
{
    MonomialSet result;
    for (const Monomial& m : monomials) {
        if (keepWithProbability(p)) {
            result.insert(m);
        }
    }
    return result;
}

// Description:
//  Creates the cross product with left and right 
//  MonomialSet
// Input: 
//  MonomialSet& left
//  MonomialSet& right
// Output:
//  vector<Monomial> vector with elements repr.
//    the cross product results
std::vector<Monomial> cross(const MonomialSet& left, const MonomialSet& right)
{
    MonomialSet products;

    for (const Monomial& leftMonomial : left) {
        for (const Monomial& rightMonomial : right) {
            Monomial product = leftMonomial;
            product.vars.insert(product.vars.end(), rightMonomial.vars.begin(), rightMonomial.vars.end());
            std::sort(product.vars.begin(), product.vars.end());
            products.insert(product);
        }
    }

    return std::vector<Monomial>(products.begin(), products.end());
}

enum class SampleStage : std::uint64_t {
    ExactSupport = 1,
    MulFinal = 2,
    PlusChildReduce = 3,
    PlusFinal = 4,
};

struct SampleKey {
    NodeId node = EMPTY_NODE;
    std::size_t sample = 0;
    SampleStage stage = SampleStage::ExactSupport;
    std::size_t ordinal = 0;
};

std::uint64_t mix64(std::uint64_t value)
{
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

void hashCombine64(std::uint64_t& seed, std::uint64_t value)
{
    seed ^= mix64(value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2));
}

std::uint64_t hashSampleKeyState(
    std::uint64_t runSeed,
    const SampleKey& key)
{
    std::uint64_t seed = mix64(runSeed);
    hashCombine64(seed, static_cast<std::uint64_t>(key.node));
    hashCombine64(seed, static_cast<std::uint64_t>(key.sample));
    hashCombine64(seed, static_cast<std::uint64_t>(key.stage));
    hashCombine64(seed, static_cast<std::uint64_t>(key.ordinal));
    return seed;
}

std::uint64_t hashSampleDecision(
    std::uint64_t runSeed,
    const SampleKey& key,
    const Monomial& monomial)
{
    std::uint64_t seed = hashSampleKeyState(runSeed, key);
    hashCombine64(seed, static_cast<std::uint64_t>(monomial.vars.size()));
    for (VarId var : monomial.vars) {
        hashCombine64(seed, static_cast<std::uint64_t>(var));
    }
    return mix64(seed);
}

bool keepWithKey(
    std::uint64_t runSeed,
    const SampleKey& key,
    const Monomial& monomial,
    double p)
{
    if (p <= 0.0 || !std::isfinite(p)) {
        return false;
    }
    if (p >= 1.0) {
        return true;
    }

    const std::uint64_t bits = hashSampleDecision(runSeed, key, monomial);
    const double unit =
        static_cast<double>(bits >> 11) * (1.0 / 9007199254740992.0);
    return unit < p;
}

class DeterministicUniformStream {
public:
    DeterministicUniformStream(std::uint64_t runSeed, const SampleKey& key)
        : seed_(mix64(hashSampleKeyState(runSeed, key)))
    {
    }

    long double nextOpenUnit()
    {
        if (counter_ == std::numeric_limits<std::uint64_t>::max()) {
            throw std::overflow_error("deterministic sample stream exhausted");
        }
        const std::uint64_t bits = mix64(
            seed_ + 0x9e3779b97f4a7c15ULL * (++counter_));
        const std::uint64_t mantissa = bits >> 11;
        constexpr long double unitCount = 9007199254740992.0L;
        return (static_cast<long double>(mantissa) + 1.0L) / (unitCount + 1.0L);
    }

private:
    std::uint64_t seed_ = 0;
    std::uint64_t counter_ = 0;
};

class GeometricIndexGenerator {
public:
    GeometricIndexGenerator(
        std::size_t count,
        double probability,
        std::uint64_t runSeed,
        const SampleKey& key)
        : count_(count),
          stream_(runSeed, key)
    {
        if (count == 0 || !std::isfinite(probability) || !(probability > 0.0)) {
            mode_ = Mode::Empty;
        } else if (probability >= 1.0) {
            mode_ = Mode::All;
        } else {
            mode_ = Mode::Sparse;
            logFailureProbability_ = std::log1p(-static_cast<long double>(probability));
        }
    }

    std::optional<std::size_t> next()
    {
        if (nextCandidate_ >= count_ || mode_ == Mode::Empty) {
            return std::nullopt;
        }
        if (mode_ == Mode::All) {
            return nextCandidate_++;
        }

        const long double unit = stream_.nextOpenUnit();
        ++uniformDraws_;
        const long double gap = std::floor(std::log(unit) / logFailureProbability_);
        const std::size_t remaining = count_ - nextCandidate_;
        if (!std::isfinite(gap) || gap < 0.0L || gap >= static_cast<long double>(remaining)) {
            nextCandidate_ = count_;
            return std::nullopt;
        }

        const std::size_t index = nextCandidate_ + static_cast<std::size_t>(gap);
        nextCandidate_ = index + 1;
        return index;
    }

    std::size_t uniformDraws() const
    {
        return uniformDraws_;
    }

private:
    enum class Mode {
        Empty,
        All,
        Sparse,
    };

    std::size_t count_ = 0;
    std::size_t nextCandidate_ = 0;
    std::size_t uniformDraws_ = 0;
    long double logFailureProbability_ = 0.0L;
    Mode mode_ = Mode::Empty;
    DeterministicUniformStream stream_;
};

MonomialSet reduceExactSupportSparse(
    const std::vector<Monomial>& monomials,
    double probability,
    const CountCoreState& st,
    const SampleKey& key)
{
    MonomialSet result;
    GeometricIndexGenerator generator(
        monomials.size(),
        probability,
        st.runSeed,
        key);
    while (const std::optional<std::size_t> index = generator.next()) {
        result.insert(monomials[*index]);
    }
    return result;
}

std::size_t checkedCartesianProductSize(std::size_t leftSize, std::size_t rightSize)
{
    if (rightSize != 0 &&
        leftSize > std::numeric_limits<std::size_t>::max() / rightSize) {
        throw std::overflow_error("lazy Cartesian product size overflows size_t");
    }
    return leftSize * rightSize;
}

std::vector<Monomial> canonicalMonomials(const MonomialSet& monomials)
{
    std::vector<Monomial> result(monomials.begin(), monomials.end());
    canonicalizeSupport(result);
    return result;
}

Monomial multiplyMonomials(const Monomial& left, const Monomial& right)
{
    Monomial product;
    product.vars.reserve(left.vars.size() + right.vars.size());
    product.vars.insert(product.vars.end(), left.vars.begin(), left.vars.end());
    product.vars.insert(product.vars.end(), right.vars.begin(), right.vars.end());
    std::sort(product.vars.begin(), product.vars.end());
    return product;
}

struct LazyProductSample {
    MonomialSet monomials;
    std::size_t candidateCount = 0;
    std::size_t materializedCount = 0;
    std::size_t uniformDraws = 0;
};

LazyProductSample reduceLazyProductSparse(
    const MonomialSet& left,
    const MonomialSet& right,
    double probability,
    const CountCoreState& st,
    const SampleKey& key)
{
    const std::vector<Monomial> indexedLeft = canonicalMonomials(left);
    const std::vector<Monomial> indexedRight = canonicalMonomials(right);

    LazyProductSample result;
    result.candidateCount = checkedCartesianProductSize(
        indexedLeft.size(),
        indexedRight.size());

    GeometricIndexGenerator generator(
        result.candidateCount,
        probability,
        st.runSeed,
        key);
    while (const std::optional<std::size_t> flatIndex = generator.next()) {
        const std::size_t leftIndex = *flatIndex / indexedRight.size();
        const std::size_t rightIndex = *flatIndex % indexedRight.size();
        result.monomials.insert(multiplyMonomials(
            indexedLeft[leftIndex],
            indexedRight[rightIndex]));
        ++result.materializedCount;
    }
    result.uniformDraws = generator.uniformDraws();
    return result;
}

MonomialSet reduceDeterministic(
    const std::vector<Monomial>& monomials,
    double p,
    const CountCoreState& st,
    const SampleKey& key)
{
    MonomialSet result;
    for (const Monomial& monomial : monomials) {
        if (keepWithKey(st.runSeed, key, monomial, p)) {
            result.insert(monomial);
        }
    }
    return result;
}

MonomialSet reduceDeterministic(
    const MonomialSet& monomials,
    double p,
    const CountCoreState& st,
    const SampleKey& key)
{
    MonomialSet result;
    for (const Monomial& monomial : monomials) {
        if (keepWithKey(st.runSeed, key, monomial, p)) {
            result.insert(monomial);
        }
    }
    return result;
}

std::uint64_t chooseRunSeed()
{
    std::uniform_int_distribution<std::uint64_t> distribution;
    return distribution(samplingGenerator());
}

class MissingProgress {
public:
    MissingProgress(const char* label, NodeId q, std::size_t total, bool enabled)
        : total_(total),
          label_(label),
          enabled_(enabled && total >= 100),
          currentNodeLabel_(q),
          started_(Clock::now()),
          lastPrint_(started_)
    {
        if (enabled_) {
            print(0, true);
        }
    }

    ~MissingProgress()
    {
        if (enabled_ && !finished_) {
            std::cerr << std::endl;
        }
    }

    void update(std::size_t completed)
    {
        if (!enabled_ || total_ == 0) {
            return;
        }

        const std::size_t doneOutOf100 = std::min<std::size_t>(
            100,
            static_cast<std::size_t>(
                (static_cast<long double>(completed) * 100.0L) /
                static_cast<long double>(total_)));
        const std::size_t missingOutOf100 = 100 - doneOutOf100;
        const auto now = Clock::now();
        const bool heartbeat =
            std::chrono::duration_cast<std::chrono::seconds>(now - lastPrint_).count() >= 1;
        if (missingOutOf100 != lastMissing_ || heartbeat) {
            lastMissing_ = missingOutOf100;
            lastPrint_ = now;
            print(completed, false);
        }
    }

    void finish()
    {
        if (enabled_) {
            print(total_, true);
            std::cerr << std::endl;
            finished_ = true;
        }
    }

private:
    using Clock = std::chrono::steady_clock;

    void print(std::size_t completed, bool final)
    {
        const long double completedLong = static_cast<long double>(completed);
        const long double totalLong = static_cast<long double>(total_);
        const long double percentDone =
            totalLong > 0.0L ? (completedLong * 100.0L) / totalLong : 100.0L;
        const long double missingOutOf100 =
            totalLong > 0.0L ? 100.0L - percentDone : 0.0L;

        const auto now = Clock::now();
        const long double elapsedSeconds =
            std::max<long double>(
                1.0e-9L,
                std::chrono::duration<long double>(now - started_).count());
        const long double rate = completedLong / elapsedSeconds;

        std::cerr << "\r" << label_ << " node " << currentNodeLabel_
                  << ": " << std::fixed << std::setprecision(9)
                  << missingOutOf100 << "/100 missing, "
                  << completed << "/" << total_ << " done ("
                  << percentDone << "%), "
                  << std::setprecision(2) << static_cast<double>(rate) << "/s";

        if (completed > 0 && completed < total_ && rate > 0.0L) {
            const long double remainingSeconds = (totalLong - completedLong) / rate;
            std::cerr << ", eta " << std::setprecision(0)
                      << static_cast<double>(remainingSeconds) << "s";
        } else if (final) {
            std::cerr << ", complete";
        }
        std::cerr << std::defaultfloat << std::flush;
    }

    std::size_t total_ = 0;
    const char* label_ = "";
    bool enabled_ = false;
    bool finished_ = false;
    std::size_t lastMissing_ = 101;
    NodeId currentNodeLabel_ = EMPTY_NODE;
    Clock::time_point started_;
    Clock::time_point lastPrint_;
};

const char* nodeKindName(const PTNode& node)
{
    if (node.kind == PTKind::Var) {
        return "var";
    }
    if (node.kind == PTKind::Add) {
        return "add";
    }
    return "mul";
}

// Description:
//  Populates a vector with the the nodes in the plus,times 
//  program via bottom up dfs
// Input
//  PlusTimesProgram& P the pt program
//  NodeId q node from which iteration starts
//  vector<char>& seen already encountered nodes
//  vector<NodeId>& order vector to be populated
void nodesBottomUpDfs(const PlusTimesProgram& P, NodeId q, std::vector<char>& seen, std::vector<NodeId>& order)
{
    if (q == EMPTY_NODE || q >= P.nodes.size() || seen[q]) {
        return;
    }

    seen[q] = true;
    const std::vector<NodeId>& children = P.nodes[q].children;
    for (auto it = children.rbegin(); it != children.rend(); ++it) {
        nodesBottomUpDfs(P, *it, seen, order);
    }
    order.push_back(q);
}

// Description:
//  Wrapper over nodesBottomUpDfs.
//  Creates a vector with elements sorted in dfs fashion
// Input:
//  PlusTimesProgram& P the plus,times program
// Output: 
//  vector<NodeId> the populated vector
std::vector<NodeId> nodesBottomUp(const PlusTimesProgram &P)
{
    std::vector<NodeId> order;
    std::vector<char> seen(P.nodes.size(), false);
    nodesBottomUpDfs(P, P.root, seen, order);
    return order;
}

// Description:
//  Creates a vector of sets containining all children sets (if node not var)
//  or its value (if node is var).
//  The output is sorted as bottom up dfs would order the nodes.
// Input:
//  PlusTimesProgram& P
// Output:
//  VariableSets (alias of vec<un_set<uint>>)
VariableSets computeVariableSets(const PlusTimesProgram& P)
{
    VariableSets variableSets(P.nodes.size());

    for (NodeId q : nodesBottomUp(P)) {
        const PTNode& node = P.nodes[q];
        if (node.kind == PTKind::Var) {
            variableSets[q].insert(static_cast<VarId>(node.var));
        } else {
            for (NodeId child : node.children) {
                variableSets[q].insert(variableSets[child].begin(), variableSets[child].end());
            }
        }
    }

    return variableSets;
}

// Description:
//  Recursively descends down the tree to check whether the monomial
//  is in the support of the specified node
// Input:
//  PlusTimesProgram& P 
//  NodeId q 
//  Monomial& monomial
//  VariableSets& variableSets
//  SupportMembershipMemo& memo
// Output:
//  bool
bool isInSupport(
    const PlusTimesProgram& P,
    NodeId q,
    const Monomial& monomial,
    const VariableSets& variableSets,
    SupportMembershipMemo& memo)
{
    if (q == EMPTY_NODE || q >= P.nodes.size()) {
        return false;
    }

    auto& nodeMemo = memo[q];
    const auto memoIt = nodeMemo.find(monomial);
    if (memoIt != nodeMemo.end()) {
        return memoIt->second;
    }

    const PTNode& node = P.nodes[q];
    bool result = false;

    if (node.kind == PTKind::Var) {
        result = monomial.vars.size() == 1 &&
                 monomial.vars.front() == static_cast<VarId>(node.var);
    } else if (node.kind == PTKind::Add) {
        for (NodeId child : node.children) {
            if (isInSupport(P, child, monomial, variableSets, memo)) {
                result = true;
                break;
            }
        }
    } else {
        assert(node.children.size() == 2);
        const NodeId left = node.children[0];
        const NodeId right = node.children[1];
        const auto& leftVars = variableSets[left];
        const auto& rightVars = variableSets[right];
        Monomial leftMonomial;
        Monomial rightMonomial;

        result = true;
        for (VarId var : monomial.vars) {
            const bool inLeft = leftVars.find(var) != leftVars.end();
            const bool inRight = rightVars.find(var) != rightVars.end();
            if (inLeft == inRight) {
                result = false;
                break;
            }
            if (inLeft) {
                leftMonomial.vars.push_back(var);
            } else {
                rightMonomial.vars.push_back(var);
            }
        }

        result = result &&
                 isInSupport(P, left, leftMonomial, variableSets, memo) &&
                 isInSupport(P, right, rightMonomial, variableSets, memo);
    }

    nodeMemo.emplace(monomial, result);
    return result;
}

// Description:
//  For each of the children in the vector, 
//  checks whether the monomials contained appear 
//  in any of the support of nodes that have priority over it.
//  If not, it is added to the resulting union. Ensuring that 
//  each monomial can contribute once only 
// Input: 
//  PlusTimesProgram& P
//  vector<NodeId>& children
//  vector<MonomialSet>& childSamples
//  VariableSets& variableSets
//  SupportMembershipMemo& memo
// Output:
//  MonomialSet set of monomial result of union
MonomialSet sampleUnion(
    const PlusTimesProgram& P,
    const std::vector<NodeId>& children,
    const std::vector<MonomialSet>& childSamples,
    const VariableSets& variableSets,
    SupportMembershipMemo& memo)
{
    MonomialSet result;

    for (std::size_t i = 0; i < children.size(); ++i) {
        for (const Monomial& monomial : childSamples[i]) {
            bool appearsEarlier = false;
            for (std::size_t j = 0; j < i; ++j) {
                if (isInSupport(P, children[j], monomial, variableSets, memo)) {
                    appearsEarlier = true;
                    break;
                }
            }

            if (!appearsEarlier) {
                result.insert(monomial);
            }
        }
    }

    return result;
}

double median(std::vector<double> values)
{
    if (values.empty()) {
        return 0.0;
    }

    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
}

bool isTimes(PTNode n)
{
  return n.kind == PTKind::Mul;
}

bool isPlus(PTNode n)
{
  return n.kind == PTKind::Add;
}

void requireSampleReady(const CountCoreState& st, NodeId q)
{
    if (q == EMPTY_NODE || q >= st.pReady.size() || !st.pReady[q]) {
        throw std::logic_error("sample regeneration requested before node estimate was ready");
    }
}

double childMinimumProbability(const PlusTimesProgram& P, NodeId q, const CountCoreState& st)
{
    const auto& children = P.nodes[q].children;
    if (children.empty()) {
        throw std::logic_error("cannot compute child probability minimum for a leaf node");
    }

    double rho = std::numeric_limits<double>::infinity();
    for (NodeId child : children) {
        requireSampleReady(st, child);
        rho = std::min(rho, st.p[child]);
    }
    return rho;
}

const MonomialSet& generateSampleMemoized(
    const PlusTimesProgram& P,
    NodeId q,
    std::size_t r,
    CountCoreState& st,
    SampleMemo& memo);

MonomialSet generatePlusHat(
    const PlusTimesProgram& P,
    NodeId q,
    std::size_t r,
    double rho,
    CountCoreState& st,
    SampleMemo& memo)
{
    if (q == EMPTY_NODE || q >= P.nodes.size()) {
        throw std::logic_error("invalid plus-times node id during sample regeneration");
    }

    const auto& children = P.nodes[q].children;
    std::vector<MonomialSet> reducedChildren;
    reducedChildren.reserve(children.size());

    for (std::size_t i = 0; i < children.size(); ++i) {
        const NodeId child = children[i];
        requireSampleReady(st, child);
        const double reduceProbability = st.p[child] > 0.0 ? rho / st.p[child] : 0.0;
        reducedChildren.push_back(reduceDeterministic(
            generateSampleMemoized(P, child, r, st, memo),
            reduceProbability,
            st,
            SampleKey{q, r, SampleStage::PlusChildReduce, i}));
    }

    return sampleUnion(P, children, reducedChildren, st.variableSets, st.membershipMemo);
}

const MonomialSet& generateSampleMemoized(
    const PlusTimesProgram& P,
    NodeId q,
    std::size_t r,
    CountCoreState& st,
    SampleMemo& memo)
{
    if (q == EMPTY_NODE || q >= P.nodes.size()) {
        throw std::logic_error("invalid plus-times node id during sample regeneration");
    }
    requireSampleReady(st, q);
    if (const MonomialSet* cached = memo.find(q, r)) {
        return *cached;
    }

    if (st.exactSupports[q].has_value()) {
        return memo.store(q, r, reduceExactSupportSparse(
            *st.exactSupports[q],
            st.p[q],
            st,
            SampleKey{q, r, SampleStage::ExactSupport, 0}));
    }

    const PTNode& node = P.nodes[q];
    if (node.kind == PTKind::Mul) {
        assert(node.children.size() == 2);
        const NodeId left = node.children[0];
        const NodeId right = node.children[1];
        requireSampleReady(st, left);
        requireSampleReady(st, right);

        const MonomialSet& leftSample = generateSampleMemoized(P, left, r, st, memo);
        const MonomialSet& rightSample = generateSampleMemoized(P, right, r, st, memo);
        const double denominator = st.p[left] * st.p[right];
        const double reduceProbability =
            denominator > 0.0 ? st.p[q] / denominator : 0.0;
        LazyProductSample productSample = reduceLazyProductSparse(
            leftSample,
            rightSample,
            reduceProbability,
            st,
            SampleKey{q, r, SampleStage::MulFinal, 0});
        return memo.store(q, r, std::move(productSample.monomials));
    }

    if (node.kind == PTKind::Add) {
        const double rho = childMinimumProbability(P, q, st);
        MonomialSet shat = generatePlusHat(P, q, r, rho, st, memo);
        const double reduceProbability = rho > 0.0 ? st.p[q] / rho : 0.0;
        return memo.store(q, r, reduceDeterministic(
            shat,
            reduceProbability,
            st,
            SampleKey{q, r, SampleStage::PlusFinal, 0}));
    }

    throw std::logic_error("large leaf support should be impossible");
}

MonomialSet generateSample(
    const PlusTimesProgram& P,
    NodeId q,
    std::size_t r,
    CountCoreState& st)
{
    SampleMemo memo(P.nodes.size());
    memo.beginSample(r);
    return generateSampleMemoized(P, q, r, st, memo);
}

// Description
//  Sample times implementation.
//  Verbatim from paper
// Input:
//  PlusTimesProgram &P 
//  NodeId q 
//  CountCoreState &st
//  size_t n
//  size_t ns
//  size_t nt
void estimateSampleTimes(
    const PlusTimesProgram &P,
    NodeId q,
    CountCoreState& st,
    std::size_t n)
{
  double p = 0.0;
  auto children = P.nodes[q].children;
  assert(children.size() == 2);
  auto left = children[0]; //q1
  auto right = children[1]; // q2
  auto p1 = st.p[left];
  auto p2 = st.p[right];
  if (std::abs(p1 - 1.0) < 1e-7){
    const std::size_t sampleSize = generateSample(P, left, 0, st).size();
    p = sampleSize > 0 ? roundDown(P, st, q, p2 / sampleSize) : 0.0;
  } else if (std::abs(p2 - 1.0) < 1e-7) {
    const std::size_t sampleSize = generateSample(P, right, 0, st).size();
    p = sampleSize > 0 ? roundDown(P, st, q, p1 / sampleSize) : 0.0;
  } else {
    p = roundDown(P, st, q, (p1 * p2) / (16 * n));
  }
  st.p[q] = p;
}

// Description
//  esimator of the samples when union is performed 
struct PlusEstimateResult {
    bool needsFinalThetaPass = true;
};

PlusEstimateResult estimateSamplePlus(
    const PlusTimesProgram &P,
    NodeId q,
    CountCoreState& st,
    std::size_t n,
    std::size_t ns,
    std::size_t nt,
    std::size_t m,
    bool verbose,
    const SampleSizeThreshold& theta)
{
  const auto& children = P.nodes[q].children;
  assert(!children.empty());

  double rho = childMinimumProbability(P, q, st);

  std::vector<double> blockMeans;
  blockMeans.reserve(nt);
  MissingProgress progress("plus estimate", q, m, verbose);
  SampleMemo sampleMemo(P.nodes.size());
  bool sawPotentialThetaViolation = false;
  std::size_t processed = 0;
  for (std::size_t j = 0; j < nt; ++j) {
    double sampleSum = 0.0;
    for (std::size_t offset = 0; offset < ns; ++offset) {
        const std::size_t r = j * ns + offset;
        sampleMemo.beginSample(r);
        MonomialSet shat = generatePlusHat(P, q, r, rho, st, sampleMemo);
        sampleSum += static_cast<double>(shat.size());
        sawPotentialThetaViolation =
            sawPotentialThetaViolation || theta.reached(shat.size());
        ++processed;
        progress.update(processed);
    }
    blockMeans.push_back(rho > 0.0 ? sampleSum / (rho * ns) : 0.0);
  }
  progress.finish();

  const double medianMean = median(blockMeans);
  const double rhoHat = medianMean > 0.0 ? (16.0 * n) / medianMean : rho;
  st.p[q] = roundDown(P, st, q, std::min(rho, rhoHat));
  const bool finalSamplesAreEmpty = st.p[q] <= 0.0 && !theta.zeroAlwaysReaches();
  return PlusEstimateResult{sawPotentialThetaViolation && !finalSamplesAreEmpty};
}

std::size_t checkedMultiply(std::size_t left, std::size_t right, const char* name);

std::size_t supportThresholdFor(std::size_t n, std::size_t programSize)
{
    return checkedMultiply(
        checkedMultiply(
            checkedMultiply(supportThresholdFactor, n, "support threshold overflows size_t"),
            programSize,
            "support threshold overflows size_t"),
        programSize,
        "support threshold overflows size_t");
}

std::size_t supportThreshold(const PlusTimesProgram& P)
{
    return supportThresholdFor(
        static_cast<std::size_t>(P.degree),
        P.nodes.size());
}

double countCoreWithSupportThreshold(
    const PlusTimesProgram& P,
    std::size_t ns,
    std::size_t nt,
    cpp_int theta,
    double kappa,
    std::size_t supportThreshold
) {
    const NodeId o = P.root;
    const std::size_t n = static_cast<std::size_t>(P.degree);
    const std::size_t m = checkedMultiply(ns, nt, "sample count overflows size_t");
    const SampleSizeThreshold sampleSizeThreshold(theta);
    const std::size_t enumerationThreshold =
        supportThreshold == std::numeric_limits<std::size_t>::max()
            ? supportThreshold
            : supportThreshold + 1;
    const bool verbose = m >= 100;

    if (verbose) {
        std::cerr << "countCore: checking root support up to "
                  << enumerationThreshold << " monomials" << std::endl;
    }
    Support rootSupp = enumerateSupportBounded(P, o, enumerationThreshold);

    if (rootSupp.has_value()) {
        std::cout << "enumerating" << std::endl;
        return static_cast<double>(rootSupp->size());
    }
    if (verbose) {
        std::cerr << "countCore: root support exceeds exact threshold; starting regenerated sampling"
                  << std::endl;
    }

    const std::vector<NodeId> order = nodesBottomUp(P);

    CountCoreState st;
    st.n = n;
    st.kappa = kappa;
    st.runSeed = chooseRunSeed();
    st.p.resize(P.nodes.size());
    st.pReady.resize(P.nodes.size(), false);
    st.exactSupports.resize(P.nodes.size());
    st.variableSets = computeVariableSets(P);
    st.membershipMemo.resize(P.nodes.size());

    if (verbose) {
        std::cerr << "countCore: " << order.size() << " nodes, "
                  << m << " regenerated samples per node" << std::endl;
    }

    for (std::size_t qi = 0; qi < order.size(); ++qi) {
        const NodeId q = order[qi];
        const PTNode& node = P.nodes[q];
        bool needsFinalThetaPass = true;
        if (verbose) {
            std::cerr << "node " << (qi + 1) << "/" << order.size()
                      << " id " << q << " (" << nodeKindName(node)
                      << "): enumerating support up to "
                      << enumerationThreshold << std::endl;
        }

        auto supp = enumerateSupportBounded(P, q, enumerationThreshold);

        if (supp.has_value()) {
            canonicalizeSupport(*supp);
            st.p[q] = supp->empty() ? 0.0 : std::min(1.0, (16.0 * n) / supp->size());
            if (verbose) {
                std::cerr << "node " << q << ": exact support size "
                          << supp->size() << ", p=" << st.p[q] << std::endl;
            }
            st.exactSupports[q] = std::move(supp);
            st.pReady[q] = true;
            const ExactThetaDecision thetaDecision = exactThetaDecision(
                st.exactSupports[q]->size(),
                st.p[q],
                sampleSizeThreshold,
                m);
            if (thetaDecision == ExactThetaDecision::Fail) {
                return 0.0;
            }
            needsFinalThetaPass = thetaDecision == ExactThetaDecision::CheckSamples;
        } else {
            if (verbose) {
                std::cerr << "node " << q << ": support is larger than threshold; estimating "
                          << nodeKindName(node) << " node" << std::endl;
            }
            if (isTimes(P.nodes[q])) {
                estimateSampleTimes(P, q, st, n);
                st.pReady[q] = true;
            } else if (isPlus(P.nodes[q])) {
                const PlusEstimateResult plusResult = estimateSamplePlus(
                    P,
                    q,
                    st,
                    n,
                    ns,
                    nt,
                    m,
                    verbose,
                    sampleSizeThreshold);
                st.pReady[q] = true;
                needsFinalThetaPass = plusResult.needsFinalThetaPass;
            } else {
                throw std::runtime_error("large leaf support should be impossible");
            }
            if (verbose) {
                std::cerr << "node " << q << ": estimated p=" << st.p[q]
                          << ", degree=" << node.degree << std::endl;
            }
        }

        if (!needsFinalThetaPass) {
            if (verbose) {
                std::cerr << "node " << q
                          << ": final samples are provably below theta; skipping regeneration"
                          << std::endl;
            }
            continue;
        }

        if (verbose) {
            std::cerr << "node " << q << ": checking final samples against theta" << std::endl;
        }
        MissingProgress progress("theta check", q, m, verbose);
        SampleMemo sampleMemo(P.nodes.size());
        for (std::size_t r = 0; r < m; ++r) {
            sampleMemo.beginSample(r);
            if (sampleSizeThreshold.reached(
                    generateSampleMemoized(P, q, r, st, sampleMemo).size())) {
                return 0.0;
            }
            progress.update(r + 1);
        }
        progress.finish();
    }

    return st.p[o] > 0.0 ? (16.0 * n) / st.p[o] : 0.0;
}

double countCore(
    const PlusTimesProgram& P,
    std::size_t ns,
    std::size_t nt,
    cpp_int theta,
    double kappa
) {
    return countCoreWithSupportThreshold(P, ns, nt, theta, kappa, supportThreshold(P));
}

std::size_t ceilToSize(double value, const char* name)
{
    if (!std::isfinite(value) || value < 0.0 ||
        value > static_cast<double>(std::numeric_limits<std::size_t>::max())) {
        throw std::overflow_error(name);
    }
    return static_cast<std::size_t>(std::ceil(value));
}

std::size_t checkedMultiply(std::size_t left, std::size_t right, const char* name)
{
    if (right != 0 && left > std::numeric_limits<std::size_t>::max() / right) {
        throw std::overflow_error(name);
    }
    return left * right;
}

double counter(const PlusTimesProgram& P, double epsilon, double delta)
{
    if (!(epsilon > 0.0)) {
        throw std::invalid_argument("epsilon must be positive");
    }
    if (!(delta > 0.0 && delta < 1.0)) {
        throw std::invalid_argument("delta must be between 0 and 1");
    }
    if (P.root == EMPTY_NODE) {
        return 0.0;
    }
    assertReadyForFPRAS(P);

    const std::size_t n = static_cast<std::size_t>(P.degree);
    const std::size_t Psize = P.nodes.size();
    const std::size_t exactThreshold = supportThreshold(P);
    if (exactThreshold < std::numeric_limits<std::size_t>::max()) {
        if (exactThreshold >= 1000) {
            std::cerr << "counter: checking exact root support up to "
                      << (exactThreshold + 1) << " monomials" << std::endl;
        }
        std::optional<std::size_t> exactSupportSize =
            supportSizeBounded(P, P.root, exactThreshold + 1);

        if (exactSupportSize.has_value()) {
            if (exactThreshold >= 1000) {
                std::cerr << "counter: exact root support found: "
                          << *exactSupportSize << std::endl;
            }
            return static_cast<double>(*exactSupportSize);
        }
    }
    if (!theoremAllowsRandomizedCounting(n)) {
        std::optional<std::size_t> exactSupportSize = supportSizeBounded(
            P,
            P.root,
            std::numeric_limits<std::size_t>::max());
        if (!exactSupportSize.has_value()) {
            throw std::overflow_error("exact support size overflows size_t");
        }
        return static_cast<double>(*exactSupportSize);
    }
    if (exactThreshold >= 1000) {
        std::cerr << "counter: exact root support exceeds threshold; starting countCore repetitions"
                  << std::endl;
    }
#ifdef DEBUG
    const double epsilonPrime = std::min(epsilon, 1.);
#else
    const double epsilonPrime = std::min(epsilon, 0.25);
#endif
    const double kappa = epsilonPrime / (4.0 * static_cast<double>(n));

#ifdef TESTING 
    const std::size_t ns = 4;
    const std::size_t nt = 8;
#else
    const std::size_t ns = ceilToSize(48.0 / (kappa * kappa), "ns overflows size_t");
    const std::size_t nt = checkedMultiply(checkedMultiply(8, n, "nt overflows size_t"), Psize, "nt overflows size_t");
#endif
    /*
    const std::size_t theta = checkedMultiply(
        checkedMultiply(checkedMultiply(checkedMultiply(512, ns, "theta overflows size_t"), nt, "theta overflows size_t"), n, "theta overflows size_t"),
        Psize,
        "theta overflows size_t");
        */
    const cpp_int theta = cpp_int(512) * ns * nt * n * Psize;
    const std::size_t repetitions = std::max<std::size_t>(
        1,
        ceilToSize(16.0 * std::log(1.0 / delta), "repetitions overflow size_t"));

    std::vector<double> estimates;
    estimates.reserve(repetitions);
    for (std::size_t i = 0; i < repetitions; ++i) {
        estimates.push_back(countCore(P, ns, nt, theta, kappa));
    }

    return median(estimates);
}
