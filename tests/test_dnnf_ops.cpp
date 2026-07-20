#include "../dnnf/DNNFParser.hpp"
#include "../utils/dnnf_ops.hpp"

#include <cassert>
#include <cstddef>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace {

bool literalValue(DNNFLiteral literal, const std::vector<bool>& assignment)
{
    const std::size_t atom = static_cast<std::size_t>(literal < 0 ? -literal : literal);
    const bool value = assignment.at(atom);
    return literal > 0 ? value : !value;
}

bool evaluateNode(
    const DNNF& dnnf,
    DNNFNodeId id,
    const std::vector<bool>& assignment,
    const std::unordered_map<DNNFNodeId, std::size_t>& indexes)
{
    const DNNFNode& node = dnnf.nodes.at(indexes.at(id));
    if (node.kind == TrueVar) {
        return true;
    }
    if (node.kind == FalseVar) {
        return false;
    }

    const bool identity = node.kind == And;
    bool result = identity;
    for (const DNNFArc& arc : node.children) {
        bool branch = evaluateNode(dnnf, arc.target, assignment, indexes);
        for (DNNFLiteral literal : arc.literals) {
            branch = branch && literalValue(literal, assignment);
        }
        if (node.kind == Or) {
            result = result || branch;
        } else {
            result = result && branch;
        }
    }
    return result;
}

bool evaluate(const DNNF& dnnf, const std::vector<bool>& assignment)
{
    std::unordered_map<DNNFNodeId, std::size_t> indexes;
    for (std::size_t i = 0; i < dnnf.nodes.size(); ++i) {
        indexes.emplace(dnnf.nodes[i].id, i);
    }
    return evaluateNode(dnnf, dnnf.root, assignment, indexes);
}

void assertEquivalent(const DNNF& first, const DNNF& second, uint32_t variables)
{
    const uint32_t assignmentCount = 1U << variables;
    for (uint32_t mask = 0; mask < assignmentCount; ++mask) {
        std::vector<bool> assignment(variables + 1, false);
        for (uint32_t atom = 1; atom <= variables; ++atom) {
            assignment[atom] = (mask & (1U << (atom - 1))) != 0;
        }
        assert(evaluate(first, assignment) == evaluate(second, assignment));
    }
}

DNNF parse(const char* text)
{
    std::istringstream input(text);
    return parseDNNF(input);
}

} // namespace

int main()
{
    // x1 OR (not x1 AND x2): x2 is absent from the first disjunct and its
    // negative polarity is absent from the complete DNNF.
    const DNNF nonsmooth = parse(
        "o 1 0\n"
        "t 2 0\n"
        "1 2 1 0\n"
        "1 2 -1 2 0\n");
    assert(!isSmoothDNNF(nonsmooth));

    const DNNF smoothed = smoothDNNF(nonsmooth);
    assert(isSmoothDNNF(smoothed));
    assertEquivalent(nonsmooth, smoothed, 2);
    assert(smoothed.root == nonsmooth.root);
    assert(smoothed.numVars == nonsmooth.numVars);

    // Exercise Definition 6(1) when there is no OR node in the input.
    DNNF positiveConjunction = parse(
        "a 1 0\n"
        "t 2 0\n"
        "1 2 1 2 0\n");
    const DNNF originalPositiveConjunction = positiveConjunction;
    const DNNF smoothConjunction = smooth_DNNF(positiveConjunction);
    assert(isSmoothDNNF(positiveConjunction));
    assert(isSmoothDNNF(smoothConjunction));
    assertEquivalent(originalPositiveConjunction, smoothConjunction, 2);

    // An already smooth DNNF needs no additional nodes.
    const DNNF alreadySmooth = parse(
        "o 1 0\n"
        "t 2 0\n"
        "1 2 1 0\n"
        "1 2 -1 0\n");
    const DNNF unchanged = smoothDNNF(alreadySmooth);
    assert(isSmoothDNNF(unchanged));
    assert(unchanged.nodes.size() == alreadySmooth.nodes.size());
    assertEquivalent(alreadySmooth, unchanged, 1);

    const DNNF sample = parseDNNF("dnnf/input.nnf");
    const DNNF smoothSample = smoothDNNF(sample);
    assert(isSmoothDNNF(smoothSample));
    assertEquivalent(sample, smoothSample, 3);

    return 0;
}
