#include "ptfromdnnf.hpp"

#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

struct TranslatedNode {
  NodeId node = EMPTY_NODE;
  bool isOne = false;
};

} // namespace

PlusTimesProgram compileDNNFtoPlusTimes(const DNNF& dnnf)
{
  PlusTimesProgram program;

  const std::size_t numLiteralVars = checkedProduct(
    static_cast<std::size_t>(dnnf.numVars),
    2,
    "too many DNNF literal variables");
  if (numLiteralVars > std::numeric_limits<uint32_t>::max()) {
    throw std::overflow_error("too many DNNF literal variables");
  }
  program.numVars = static_cast<uint32_t>(numLiteralVars);

  if (dnnf.root == DNNF_EMPTY_NODE) {
    if (!dnnf.nodes.empty()) {
      throw std::invalid_argument("nonempty DNNF must have a root");
    }
    return program;
  }

  std::vector<uint32_t> nodeIds;
  nodeIds.reserve(dnnf.nodes.size());
  for (const DNNFNode& node : dnnf.nodes) {
    nodeIds.push_back(node.id);
  }
  const auto nodeIndexes = makeDenseIndex(nodeIds, "DNNF node");
  const std::size_t rootIndex = lookupDense(nodeIndexes, dnnf.root, "DNNF node");

  std::vector<NodeId> literalNodes(numLiteralVars, EMPTY_NODE);
  std::vector<uint8_t> state(dnnf.nodes.size(), 0);
  std::vector<TranslatedNode> translated(dnnf.nodes.size());

  const auto translateLiteral = [&](DNNFLiteral literal) -> NodeId {
    if (literal == 0 || literal == std::numeric_limits<DNNFLiteral>::min()) {
      throw std::invalid_argument("invalid DNNF literal");
    }

    const uint32_t atom = static_cast<uint32_t>(
      literal < 0 ? -static_cast<int64_t>(literal) : literal);
    if (atom > dnnf.numVars) {
      throw std::invalid_argument("DNNF literal exceeds the declared variable count");
    }

    const std::size_t variable = variableIndex(
      literal < 0 ? 1U : 0U,
      static_cast<int>(atom - 1),
      2);
    if (literalNodes[variable] == EMPTY_NODE) {
      literalNodes[variable] = makeVar(program, static_cast<PolyVarId>(variable));
    }
    return literalNodes[variable];
  };

  std::function<TranslatedNode(std::size_t)> translateNode;
  translateNode = [&](std::size_t index) -> TranslatedNode {
    if (state[index] == 2) {
      return translated[index];
    }
    if (state[index] == 1) {
      throw std::invalid_argument("DNNF must be acyclic");
    }

    state[index] = 1;
    const DNNFNode& dnnfNode = dnnf.nodes[index];

    if (dnnfNode.kind == TrueVar) {
      translated[index].isOne = true;
    } else if (dnnfNode.kind == FalseVar) {
      translated[index] = {};
    } else {
      std::vector<TranslatedNode> branches;
      branches.reserve(dnnfNode.children.size());

      for (const DNNFArc& arc : dnnfNode.children) {
        TranslatedNode branch = translateNode(
          lookupDense(nodeIndexes, arc.target, "DNNF node"));

        for (DNNFLiteral literal : arc.literals) {
          const NodeId literalNode = translateLiteral(literal);
          if (branch.isOne) {
            branch = {literalNode, false};
          } else {
            branch.node = makeMult(program, branch.node, literalNode);
          }
        }
        branches.push_back(branch);
      }

      if (dnnfNode.kind == Or) {
        std::vector<NodeId> children;
        children.reserve(branches.size());
        for (const TranslatedNode branch : branches) {
          if (branch.isOne) {
            throw std::invalid_argument(
              "PlusTimesProgram cannot represent a constant-one sum term");
          }
          children.push_back(branch.node);
        }
        translated[index].node = makeAdd(program, children);
      } else if (dnnfNode.kind == And) {
        TranslatedNode product;
        product.isOne = true;
        for (const TranslatedNode branch : branches) {
          if (branch.isOne) {
            continue;
          }
          if (product.isOne) {
            product = branch;
          } else {
            product.node = makeMult(program, product.node, branch.node);
          }
        }
        translated[index] = product;
      } else {
        throw std::invalid_argument("unknown DNNF node kind");
      }
    }

    state[index] = 2;
    return translated[index];
  };

  const TranslatedNode root = translateNode(rootIndex);
  if (root.isOne) {
    throw std::invalid_argument(
      "PlusTimesProgram cannot represent a constant-one root");
  }

  program.root = root.node;
  program = pruneReachableTopological(program);
  if (program.root != EMPTY_NODE) {
    program.degree = program.nodes[program.root].degree;
    program = depthReduceVSBR(program);
  }
  return program;
}
