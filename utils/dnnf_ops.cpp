#include "dnnf_ops.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using Atom = uint32_t;
using AtomSet = std::set<Atom>;

struct GraphInfo {
    std::unordered_map<DNNFNodeId, std::size_t> indexes;
    std::vector<AtomSet> atoms;
};

Atom atomOf(DNNFLiteral literal)
{
    if (literal == 0) {
        throw std::invalid_argument("DNNF literals cannot be 0");
    }
    if (literal == std::numeric_limits<DNNFLiteral>::min()) {
        throw std::invalid_argument("DNNF literal magnitude cannot exceed INT32_MAX");
    }
    return static_cast<Atom>(literal < 0 ? -static_cast<int64_t>(literal) : literal);
}

std::unordered_map<DNNFNodeId, std::size_t> makeNodeIndexes(const DNNF& dnnf)
{
    std::unordered_map<DNNFNodeId, std::size_t> indexes;
    indexes.reserve(dnnf.nodes.size());

    for (std::size_t i = 0; i < dnnf.nodes.size(); ++i) {
        const DNNFNodeId id = dnnf.nodes[i].id;
        if (id == 0 || id == DNNF_EMPTY_NODE) {
            throw std::invalid_argument("DNNF node ids must be nonzero and distinct from DNNF_EMPTY_NODE");
        }
        if (!indexes.emplace(id, i).second) {
            throw std::invalid_argument("duplicate DNNF node id " + std::to_string(id));
        }
    }

    if (dnnf.root != DNNF_EMPTY_NODE && indexes.find(dnnf.root) == indexes.end()) {
        throw std::invalid_argument("DNNF root refers to an undeclared node");
    }
    return indexes;
}

AtomSet computeAtoms(
    std::size_t nodeIndex,
    const DNNF& dnnf,
    const std::unordered_map<DNNFNodeId, std::size_t>& indexes,
    std::vector<AtomSet>& memo,
    std::vector<uint8_t>& state)
{
    if (state[nodeIndex] == 2) {
        return memo[nodeIndex];
    }
    if (state[nodeIndex] == 1) {
        throw std::invalid_argument("DNNF must be acyclic");
    }

    state[nodeIndex] = 1;
    AtomSet result;
    const DNNFNode& node = dnnf.nodes[nodeIndex];
    if ((node.kind == TrueVar || node.kind == FalseVar) && !node.children.empty()) {
        throw std::invalid_argument("DNNF constant nodes cannot have children");
    }

    for (const DNNFArc& arc : node.children) {
        const auto target = indexes.find(arc.target);
        if (target == indexes.end()) {
            throw std::invalid_argument(
                "DNNF arc refers to undeclared node " + std::to_string(arc.target));
        }

        const AtomSet childAtoms = computeAtoms(target->second, dnnf, indexes, memo, state);
        result.insert(childAtoms.begin(), childAtoms.end());
        for (DNNFLiteral literal : arc.literals) {
            result.insert(atomOf(literal));
        }
    }

    memo[nodeIndex] = std::move(result);
    state[nodeIndex] = 2;
    return memo[nodeIndex];
}

GraphInfo inspectGraph(const DNNF& dnnf)
{
    GraphInfo info;
    info.indexes = makeNodeIndexes(dnnf);
    info.atoms.resize(dnnf.nodes.size());
    std::vector<uint8_t> state(dnnf.nodes.size(), 0);

    for (std::size_t i = 0; i < dnnf.nodes.size(); ++i) {
        computeAtoms(i, dnnf, info.indexes, info.atoms, state);
    }
    return info;
}

AtomSet branchAtoms(const DNNFArc& arc, const GraphInfo& info)
{
    const auto target = info.indexes.find(arc.target);
    if (target == info.indexes.end()) {
        throw std::invalid_argument("DNNF arc refers to an undeclared node");
    }

    AtomSet atoms = info.atoms[target->second];
    for (DNNFLiteral literal : arc.literals) {
        atoms.insert(atomOf(literal));
    }
    return atoms;
}

struct Polarity {
    bool positive = false;
    bool negative = false;
};

std::unordered_map<Atom, Polarity> collectPolarities(const DNNF& dnnf)
{
    std::unordered_map<Atom, Polarity> polarities;
    for (const DNNFNode& node : dnnf.nodes) {
        for (const DNNFArc& arc : node.children) {
            for (DNNFLiteral literal : arc.literals) {
                Polarity& polarity = polarities[atomOf(literal)];
                if (literal > 0) {
                    polarity.positive = true;
                } else {
                    polarity.negative = true;
                }
            }
        }
    }
    return polarities;
}

class DNNFBuilder {
public:
    explicit DNNFBuilder(DNNF& dnnf)
        : dnnf_(dnnf)
    {
        for (const DNNFNode& node : dnnf_.nodes) {
            largestId_ = std::max(largestId_, node.id);
            if (node.kind == TrueVar && trueNode_ == DNNF_EMPTY_NODE) {
                trueNode_ = node.id;
            }
            if (node.kind == FalseVar && falseNode_ == DNNF_EMPTY_NODE) {
                falseNode_ = node.id;
            }
        }
    }

    DNNFNodeId trueNode()
    {
        if (trueNode_ == DNNF_EMPTY_NODE) {
            trueNode_ = addNode(TrueVar, {});
        }
        return trueNode_;
    }

    DNNFNodeId falseNode()
    {
        if (falseNode_ == DNNF_EMPTY_NODE) {
            falseNode_ = addNode(FalseVar, {});
        }
        return falseNode_;
    }

    DNNFNodeId tautology(Atom atom)
    {
        const auto cached = tautologies_.find(atom);
        if (cached != tautologies_.end()) {
            return cached->second;
        }

        if (atom > static_cast<Atom>(std::numeric_limits<DNNFLiteral>::max())) {
            throw std::overflow_error("DNNF atom cannot be represented as an int32 literal");
        }
        const DNNFLiteral literal = static_cast<DNNFLiteral>(atom);
        const DNNFNodeId id = addNode(
            Or,
            {{trueNode(), {literal}}, {trueNode(), {-literal}}});
        tautologies_.emplace(atom, id);
        return id;
    }

    DNNFNodeId literalWithMissingNegation(DNNFLiteral literal)
    {
        // Darwiche's first operation: l becomes l OR (not l AND false).
        return addNode(
            Or,
            {{trueNode(), {literal}}, {falseNode(), {static_cast<DNNFLiteral>(-literal)}}});
    }

    DNNFNodeId conjunction(std::vector<DNNFArc> children)
    {
        return addNode(And, std::move(children));
    }

private:
    DNNFNodeId freshId()
    {
        if (largestId_ >= DNNF_EMPTY_NODE - 1) {
            throw std::overflow_error("no DNNF node ids remain for smoothing");
        }
        return ++largestId_;
    }

    DNNFNodeId addNode(DNNFKind kind, std::vector<DNNFArc> children)
    {
        const DNNFNodeId id = freshId();
        dnnf_.nodes.push_back(
            {kind, static_cast<uint32_t>(children.size()), id, std::move(children)});
        return id;
    }

    DNNF& dnnf_;
    DNNFNodeId largestId_ = 0;
    DNNFNodeId trueNode_ = DNNF_EMPTY_NODE;
    DNNFNodeId falseNode_ = DNNF_EMPTY_NODE;
    std::unordered_map<Atom, DNNFNodeId> tautologies_;
};

bool wrapOneLiteralOccurrence(DNNF& dnnf, DNNFBuilder& builder, DNNFLiteral literal)
{
    for (std::size_t nodeIndex = 0; nodeIndex < dnnf.nodes.size(); ++nodeIndex) {
        for (std::size_t arcIndex = 0; arcIndex < dnnf.nodes[nodeIndex].children.size(); ++arcIndex) {
            const DNNFArc original = dnnf.nodes[nodeIndex].children[arcIndex];
            const auto occurrence = std::find(original.literals.begin(), original.literals.end(), literal);
            if (occurrence == original.literals.end()) {
                continue;
            }

            DNNFArc remainder = original;
            remainder.literals.erase(remainder.literals.begin() + (occurrence - original.literals.begin()));
            const DNNFNodeId literalNode = builder.literalWithMissingNegation(literal);
            const DNNFNodeId conjunction = builder.conjunction(
                {std::move(remainder), {literalNode, {}}});

            // Reacquire the node after builder additions, which may reallocate nodes.
            dnnf.nodes[nodeIndex].children[arcIndex] = {conjunction, {}};
            return true;
        }
    }
    return false;
}

void normalizeDegrees(DNNF& dnnf)
{
    for (DNNFNode& node : dnnf.nodes) {
        node.degree = static_cast<uint32_t>(node.children.size());
    }
}

} // namespace

bool isSmoothDNNF(const DNNF& dnnf)
{
    if (dnnf.root == DNNF_EMPTY_NODE) {
        return dnnf.nodes.empty();
    }

    const GraphInfo info = inspectGraph(dnnf);
    const auto polarities = collectPolarities(dnnf);
    for (const auto& entry : polarities) {
        if (!entry.second.positive || !entry.second.negative) {
            return false;
        }
    }

    for (std::size_t i = 0; i < dnnf.nodes.size(); ++i) {
        const DNNFNode& node = dnnf.nodes[i];
        if (node.kind != Or) {
            continue;
        }
        for (const DNNFArc& arc : node.children) {
            if (branchAtoms(arc, info) != info.atoms[i]) {
                return false;
            }
        }
    }
    return true;
}

DNNF smoothDNNF(const DNNF& input)
{
    if (input.root == DNNF_EMPTY_NODE) {
        if (!input.nodes.empty()) {
            throw std::invalid_argument("nonempty DNNF must have a root");
        }
        return input;
    }

    // Validate the input before adding nodes and edges.
    inspectGraph(input);
    DNNF result = input;
    DNNFBuilder builder(result);

    /*
    // Definition 6(1): introduce a logically dead occurrence of every missing
    // polarity using l OR (not l AND false).
    const auto originalPolarities = collectPolarities(result);
    for (const auto& entry : originalPolarities) {
        const Atom atom = entry.first;
        const Polarity polarity = entry.second;
        if (polarity.positive && polarity.negative) {
            continue;
        }

        const DNNFLiteral occurrence = polarity.positive
            ? static_cast<DNNFLiteral>(atom)
            : -static_cast<DNNFLiteral>(atom);
        if (!wrapOneLiteralOccurrence(result, builder, occurrence)) {
            throw std::logic_error("failed to locate DNNF literal selected for smoothing");
        }
    }
    */

    // Definition 6(2): each OR branch receives A OR not A for every atom that
    // occurs in the OR but is absent from that branch.
    const GraphInfo info = inspectGraph(result);
    const std::size_t nodesBeforeTautologies = result.nodes.size();
    for (std::size_t nodeIndex = 0; nodeIndex < nodesBeforeTautologies; ++nodeIndex) {
        if (result.nodes[nodeIndex].kind != Or) {
            continue;
        }

        const std::size_t branchCount = result.nodes[nodeIndex].children.size();
        for (std::size_t arcIndex = 0; arcIndex < branchCount; ++arcIndex) {
            const DNNFArc original = result.nodes[nodeIndex].children[arcIndex];
            const AtomSet present = branchAtoms(original, info);
            std::vector<DNNFArc> conjunctionChildren;
            conjunctionChildren.push_back(original);

            for (Atom atom : info.atoms[nodeIndex]) {
                if (present.find(atom) == present.end()) {
                    conjunctionChildren.push_back({builder.tautology(atom), {}});
                }
            }

            if (conjunctionChildren.size() > 1) {
                const DNNFNodeId conjunction = builder.conjunction(std::move(conjunctionChildren));
                result.nodes[nodeIndex].children[arcIndex] = {conjunction, {}};
            }
        }
    }

    normalizeDegrees(result);
    if (!isSmoothDNNF(result)) {
        throw std::logic_error("internal error: DNNF smoothing did not produce a smooth DNNF");
    }
    return result;
}
