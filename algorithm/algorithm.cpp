#include "algorithm.hpp"
#include <algorithm>
#include <boost/multiprecision/fwd.hpp>
#include <cassert>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <optional>
#include <unordered_set>
#include <boost/container_hash/hash.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <iostream>

#define DEBUG

using boost::multiprecision::cpp_int;

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

// Struct used to store the state of the current
// counter iteration.
struct CountCoreState {
    std::size_t n = 0;
    double kappa = 0.0;
    std::vector<double> p;
    std::vector<std::vector<MonomialSet>> samples;
    std::vector<int> effectiveHeights;
};

using SupportMembershipMemo = std::vector<std::unordered_map<Monomial, bool, MonomialHash>>;

#ifndef CFGFPRAS_KEEP_ALL_SAMPLES
#define CFGFPRAS_RELEASE_UNUSED_SAMPLES
#endif

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

double roundDown(const CountCoreState& st, NodeId q, double v) {
    if (!(v > 0.0)) {
        return 0.0;
    }
    if (v > 1.0) {
        return 1.0;
    }
    if (q >= st.effectiveHeights.size() || st.effectiveHeights[q] <= 0) {
        return v;
    }

    const int height = st.effectiveHeights[q];
    const long double twoToHeight = powerOfTwo(static_cast<std::size_t>(height));
    const long double maxEll = powerOfTwo(st.n);
    const long double exponent = static_cast<long double>(st.kappa) * twoToHeight;
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

#ifdef CFGFPRAS_RELEASE_UNUSED_SAMPLES
std::vector<std::size_t> countRemainingSampleUses(
    const PlusTimesProgram& P,
    const std::vector<NodeId>& order)
{
    std::vector<std::size_t> remainingUses(P.nodes.size(), 0);
    for (NodeId q : order) {
        for (NodeId child : P.nodes[q].children) {
            if (child < remainingUses.size()) {
                ++remainingUses[child];
            }
        }
    }
    return remainingUses;
}

void allocateSampleSlots(CountCoreState& st, NodeId q, std::size_t count)
{
    st.samples[q].clear();
    st.samples[q].resize(count);
}

void releaseSampleSlots(std::vector<MonomialSet>& samples)
{
    std::vector<MonomialSet>().swap(samples);
}
#endif

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

int effectiveHeightFromChildren(const PlusTimesProgram& P, NodeId q, const CountCoreState& st)
{
    int height = 0;
    for (NodeId child : P.nodes[q].children) {
        assert(child < st.effectiveHeights.size());
        assert(st.effectiveHeights[child] >= 0);
        height = std::max(height, st.effectiveHeights[child] + 1);
    }
    return height;
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
    std::size_t n,
    std::size_t ns,
    std::size_t nt)
{
  double p = 0.0;
  auto children = P.nodes[q].children;
  assert(children.size() == 2);
  auto left = children[0]; //q1
  auto right = children[1]; // q2
  auto p1 = st.p[left];
  auto p2 = st.p[right];
  if (std::abs(p1 - 1.0) < 1e-7){
    const std::size_t sampleSize = st.samples[left][0].size();
    p = sampleSize > 0 ? roundDown(st, q, p2 / sampleSize) : 0.0;
  } else if (std::abs(p2 - 1.0) < 1e-7) {
    const std::size_t sampleSize = st.samples[right][0].size();
    p = sampleSize > 0 ? roundDown(st, q, p1 / sampleSize) : 0.0;
  } else {
    p = roundDown(st, q, (p1 * p2) / (16 * n));
  }
  st.p[q] = p;

  const std::size_t m = ns * nt;
  for (std::size_t r = 0; r < m; r++){
    const double reduceProbability = (p1 > 0.0 && p2 > 0.0) ? p / (p1 * p2) : 0.0;
    st.samples[q][r] = reduce(cross(st.samples[left][r], st.samples[right][r]), reduceProbability);
  }
}

// Description
//  esimator of the samples when union is performed 
void estimateSamplePlus(
    const PlusTimesProgram &P,
    NodeId q,
    CountCoreState& st,
    std::size_t n,
    std::size_t ns,
    std::size_t nt)
{
  const auto& children = P.nodes[q].children;
  assert(!children.empty());

  double rho = std::numeric_limits<double>::infinity();
  for (NodeId child : children) {
    rho = std::min(rho, st.p[child]);
  }

  const std::size_t m = ns * nt;
  const VariableSets variableSets = computeVariableSets(P);
  SupportMembershipMemo membershipMemo(P.nodes.size());
  std::vector<MonomialSet> shat(m);

  for (std::size_t r = 0; r < m; r++){
    std::vector<MonomialSet> reducedChildren;
    reducedChildren.reserve(children.size());
    for (NodeId child : children) {
        const double reduceProbability = st.p[child] > 0.0 ? rho / st.p[child] : 0.0;
        reducedChildren.push_back(reduce(st.samples[child][r], reduceProbability));
    }

    shat[r] = sampleUnion(P, children, reducedChildren, variableSets, membershipMemo);
  }

  std::vector<double> blockMeans;
  blockMeans.reserve(nt);
  for (std::size_t j = 0; j < nt; ++j) {
    double sampleSum = 0.0;
    for (std::size_t offset = 0; offset < ns; ++offset) {
        const std::size_t r = j * ns + offset;
        sampleSum += static_cast<double>(shat[r].size());
    }
    blockMeans.push_back(rho > 0.0 ? sampleSum / (rho * ns) : 0.0);
  }

  const double medianMean = median(blockMeans);
  const double rhoHat = medianMean > 0.0 ? (16.0 * n) / medianMean : rho;
  st.p[q] = roundDown(st, q, std::min(rho, rhoHat));

  for (std::size_t r = 0; r < m; r++){
    const double reduceProbability = rho > 0.0 ? st.p[q] / rho : 0.0;
    st.samples[q][r] = reduce(shat[r], reduceProbability);
  }
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
    const std::size_t m = ns * nt;

    Support rootSupp = enumerateSupportBounded(P, o, supportThreshold + 1);

    if (rootSupp.has_value()) {
        std::cout << "enumerating" << std::endl;
        return static_cast<double>(rootSupp->size());
    }

    const std::vector<NodeId> order = nodesBottomUp(P);

    CountCoreState st;
    st.n = n;
    st.kappa = kappa;
    st.p.resize(P.nodes.size());
#ifdef CFGFPRAS_RELEASE_UNUSED_SAMPLES
    st.samples.resize(P.nodes.size());
    std::vector<std::size_t> remainingSampleUses =
        countRemainingSampleUses(P, order);
#else
    st.samples.resize(P.nodes.size(), std::vector<MonomialSet>(m));
#endif
    st.effectiveHeights.resize(P.nodes.size(), -1);

    for (NodeId q : order) {
#ifdef CFGFPRAS_RELEASE_UNUSED_SAMPLES
        allocateSampleSlots(st, q, m);
#endif
        auto supp = enumerateSupportBounded(P, q, supportThreshold + 1);

        if (supp.has_value()) {
            st.effectiveHeights[q] = 0;
            st.p[q] = supp->empty() ? 0.0 : std::min(1.0, (16.0 * n) / supp->size());

            for (std::size_t r = 0; r < m; ++r) {
                st.samples[q][r] = reduce(*supp, st.p[q]);
            }
        } else {
            st.effectiveHeights[q] = effectiveHeightFromChildren(P, q, st);
            if (isTimes(P.nodes[q])) {
                estimateSampleTimes(P, q, st, n, ns, nt);
            } else if (isPlus(P.nodes[q])) {
                estimateSamplePlus(P, q, st, n, ns, nt);
            } else {
                return -1; // FIXME: throw
                // throw std::runtime_error("large leaf support should be impossible");
            }
        }

        for (std::size_t r = 0; r < m; ++r) {
            if (st.samples[q][r].size() >= theta) {
                return 0.0;
            }
        }

#ifdef CFGFPRAS_RELEASE_UNUSED_SAMPLES
        for (NodeId child : P.nodes[q].children) {
            assert(child < remainingSampleUses.size());
            assert(remainingSampleUses[child] > 0);
            --remainingSampleUses[child];
            if (remainingSampleUses[child] == 0) {
                releaseSampleSlots(st.samples[child]);
            }
        }
#endif
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
    const std::size_t n = static_cast<std::size_t>(P.degree);
    const std::size_t Psize = P.nodes.size();
    const std::size_t threshold = 16ull * n * Psize * Psize;
    return countCoreWithSupportThreshold(P, ns, nt, theta, kappa, threshold);
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
    const std::size_t exactThreshold = checkedMultiply(
        checkedMultiply(checkedMultiply(16, n, "support threshold overflows size_t"), Psize, "support threshold overflows size_t"),
        Psize,
        "support threshold overflows size_t");
    if (exactThreshold < std::numeric_limits<std::size_t>::max()) {
        std::optional<std::size_t> exactSupportSize =
            supportSizeBounded(P, P.root, exactThreshold + 1);

        std::cout << exactThreshold << std::endl;
        if (exactSupportSize.has_value()) {
            return static_cast<double>(*exactSupportSize);
        }
    }
#ifdef DEBUG
    const double epsilonPrime = std::min(epsilon, 1.);
#else
    const double epsilonPrime = std::min(epsilon, 0.25);
#endif
    const double kappa = epsilonPrime / (4.0 * std::pow(static_cast<double>(n + 1), 3.0));
    const std::size_t ns = ceilToSize(12.0 / (kappa * kappa), "ns overflows size_t");
    const std::size_t nt = checkedMultiply(checkedMultiply(8, n, "nt overflows size_t"), Psize, "nt overflows size_t");
    /*
    const std::size_t theta = checkedMultiply(
        checkedMultiply(checkedMultiply(checkedMultiply(512, ns, "theta overflows size_t"), nt, "theta overflows size_t"), n, "theta overflows size_t"),
        Psize,
        "theta overflows size_t");
        */
    const cpp_int theta = 512 * ns * nt * n * Psize;
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
