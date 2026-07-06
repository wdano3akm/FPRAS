#include "plustimes.hpp"
#include "../cfg/CFG.hpp"

#include <algorithm>
#include <set>
#include <stack>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

struct DenseBinaryRule {
  std::size_t left;
  std::size_t right;
};

NodeId appendNode(PlusTimesProgram& program, PTNode node)
{
  if (program.nodes.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max()) - 1) {
    throw std::overflow_error("too many plus-times nodes");
  }

  const NodeId id = static_cast<NodeId>(program.nodes.size());
  program.nodes.push_back(std::move(node));
  return id;
}

NodeId makeVar(PlusTimesProgram& program, PolyVarId var)
{
  PTNode node;
  node.kind = PTKind::Var;
  node.degree = 1;
  node.var = var;

  return appendNode(program, std::move(node));
}

void requireNode(const PlusTimesProgram& program, NodeId id)
{
  if (id == EMPTY_NODE || id >= program.nodes.size()) {
    throw std::logic_error("invalid plus-times node id");
  }
}

NodeId makeAdd(PlusTimesProgram& program, const std::vector<NodeId>& rawChildren)
{
  std::vector<NodeId> children;
  children.reserve(rawChildren.size());

  int degree = -1;

  const auto addChild = [&](NodeId child) {
    if (child == EMPTY_NODE) {
      return;
    }

    requireNode(program, child);
    const int childDegree = program.nodes[child].degree;
    if (degree == -1) {
      degree = childDegree;
    } else if (degree != childDegree) {
      throw std::logic_error("non-homogeneous plus-times addition");
    }

    children.push_back(child);
  };

  for (NodeId child : rawChildren) {
    if (child == EMPTY_NODE) {
      continue;
    }

    requireNode(program, child);
    const PTNode& node = program.nodes[child];
    if (node.kind == PTKind::Add) {
      for (NodeId grandchild : node.children) {
        addChild(grandchild);
      }
    } else {
      addChild(child);
    }
  }

  if (children.empty()) {
    return EMPTY_NODE;
  }

  std::sort(children.begin(), children.end());
  children.erase(std::unique(children.begin(), children.end()), children.end());

  if (children.size() == 1) {
    return children.front();
  }

  PTNode node;
  node.kind = PTKind::Add;
  node.degree = degree;
  node.children = std::move(children);

  return appendNode(program, std::move(node));
}

NodeId makeMult(PlusTimesProgram& program, NodeId left, NodeId right)
{
  if (left == EMPTY_NODE || right == EMPTY_NODE) {
    return EMPTY_NODE;
  }

  requireNode(program, left);
  requireNode(program, right);

  PTNode node;
  node.kind = PTKind::Mul;
  node.degree = program.nodes[left].degree + program.nodes[right].degree;
  node.children = {left, right};

  return appendNode(program, std::move(node));
}

void validateDisjointNonzeroIds(const CFG& cfg)
{
  std::unordered_set<uint32_t> nonterminals;
  nonterminals.reserve(cfg.nonterminalIds.size());
  for (uint32_t id : cfg.nonterminalIds) {
    if (id == 0) {
      throw std::invalid_argument("nonterminal id 0 is reserved");
    }
    if (!nonterminals.insert(id).second) {
      throw std::invalid_argument("duplicate nonterminal id " + std::to_string(id));
    }
  }

  std::unordered_set<uint32_t> terminals;
  terminals.reserve(cfg.terminalIds.size());
  for (uint32_t id : cfg.terminalIds) {
    if (id == 0) {
      throw std::invalid_argument("terminal id 0 is reserved");
    }
    if (nonterminals.find(id) != nonterminals.end()) {
      throw std::invalid_argument("terminal/nonterminal id overlap at " + std::to_string(id));
    }
    if (!terminals.insert(id).second) {
      throw std::invalid_argument("duplicate terminal id " + std::to_string(id));
    }
  }
}

std::unordered_map<uint32_t, std::size_t> makeDenseIndex(
  const std::vector<uint32_t>& ids,
  const char* kind)
{
  std::unordered_map<uint32_t, std::size_t> dense;
  dense.reserve(ids.size());

  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (!dense.emplace(ids[i], i).second) {
      throw std::invalid_argument(
        "duplicate " + std::string(kind) + " id " + std::to_string(ids[i]));
    }
  }

  return dense;
}

std::size_t lookupDense(
  const std::unordered_map<uint32_t, std::size_t>& dense,
  uint32_t id,
  const char* kind)
{
  const auto it = dense.find(id);
  if (it == dense.end()) {
    throw std::invalid_argument(
      "rule references unknown " + std::string(kind) + " id " + std::to_string(id));
  }

  return it->second;
}

std::size_t dpIndex(std::size_t nonterminal, int length, int offset, int n)
{
  const std::size_t side = static_cast<std::size_t>(n) + 1;
  return (nonterminal * side + static_cast<std::size_t>(length)) * side +
         static_cast<std::size_t>(offset);
}

std::size_t variableIndex(std::size_t terminal, int offset, std::size_t numTerminals)
{
  return static_cast<std::size_t>(offset) * numTerminals + terminal;
}

std::size_t checkedProduct(std::size_t left, std::size_t right, const char* message)
{
  if (right != 0 && left > std::numeric_limits<std::size_t>::max() / right) {
    throw std::overflow_error(message);
  }

  return left * right;
}

} // namespace

void makeMinimal(PlusTimesProgram &pt)
{
  if (pt.root == EMPTY_NODE) {
    pt.nodes.clear();
    return;
  }

  std::vector<char> marked(pt.nodes.size(), false);

  std::vector<NodeId> oldToNew(pt.nodes.size(), EMPTY_NODE);
  std::vector<PTNode> nnodes;

  for (NodeId old = 0; old < pt.nodes.size(); ++old) {
    if (!marked[old]) continue;
    oldToNew[old] = static_cast<NodeId>(nnodes.size());
    nnodes.push_back(pt.nodes[old]);
  }

  for (PTNode& node : nnodes) {
    for (NodeId& child : node.children) {
      child = oldToNew[child];
    }
  }

  pt.root = oldToNew[pt.root];
  pt.nodes = std::move(nnodes);
}

PlusTimesProgram compileCFGToPlusTimes(const CFG& cfg, int n)
{
  if (n < 0) {
    throw std::invalid_argument("n cannot be negative");
  }
  if (n == 0) {
    throw std::runtime_error("length-zero CFG compilation needs constant nodes");
  }

  validateDisjointNonzeroIds(cfg);

  const auto nonterminalIndex = makeDenseIndex(cfg.nonterminalIds, "nonterminal");
  const auto terminalIndex = makeDenseIndex(cfg.terminalIds, "terminal");
  const std::size_t numNonterminals = cfg.nonterminalIds.size();
  const std::size_t numTerminals = cfg.terminalIds.size();

  if (numNonterminals == 0) {
    throw std::invalid_argument("CFG must have at least one nonterminal");
  }

  const std::size_t numVars = checkedProduct(
    static_cast<std::size_t>(n),
    numTerminals,
    "too many plus-times variables");
  if (numVars > std::numeric_limits<PolyVarId>::max()) {
    throw std::overflow_error("too many plus-times variables");
  }

  std::vector<std::vector<std::size_t>> terminalRules(numNonterminals);
  for (const TerminalRule& rule : cfg.terminalRules) {
    const std::size_t lhs = lookupDense(nonterminalIndex, rule.lhs, "nonterminal");
    const std::size_t terminal = lookupDense(terminalIndex, rule.terminal, "terminal");
    terminalRules[lhs].push_back(terminal);
  }

  std::vector<std::vector<DenseBinaryRule>> binaryRules(numNonterminals);
  for (const BinaryRule& rule : cfg.binaryRules) {
    const std::size_t lhs = lookupDense(nonterminalIndex, rule.lhs, "nonterminal");
    binaryRules[lhs].push_back({
      lookupDense(nonterminalIndex, rule.left, "nonterminal"),
      lookupDense(nonterminalIndex, rule.right, "nonterminal")
    });
  }

  PlusTimesProgram program;
  program.numVars = static_cast<uint32_t>(numVars);
  program.degree = n;

  std::vector<NodeId> variableNodes(numVars, EMPTY_NODE);
  for (int offset = 0; offset < n; ++offset) {
    for (std::size_t terminal = 0; terminal < numTerminals; ++terminal) {
      const std::size_t var = variableIndex(terminal, offset, numTerminals);
      variableNodes[var] = makeVar(program, static_cast<PolyVarId>(var));
    }
  }

  const std::size_t side = static_cast<std::size_t>(n) + 1;
  const std::size_t dpSide = checkedProduct(side, side, "DP table too large");
  const std::size_t dpSize = checkedProduct(numNonterminals, dpSide, "DP table too large");
  std::vector<NodeId> dp(dpSize, EMPTY_NODE);

  for (std::size_t nonterminal = 0; nonterminal < numNonterminals; ++nonterminal) {
    for (int offset = 0; offset < n; ++offset) {
      std::vector<NodeId> terminalChildren;
      terminalChildren.reserve(terminalRules[nonterminal].size());

      for (std::size_t terminal : terminalRules[nonterminal]) {
        terminalChildren.push_back(variableNodes[variableIndex(terminal, offset, numTerminals)]);
      }

      dp[dpIndex(nonterminal, 1, offset, n)] = makeAdd(program, terminalChildren);
    }
  }

  for (int length = 2; length <= n; ++length) {
    for (int offset = 0; offset <= n - length; ++offset) {
      for (std::size_t nonterminal = 0; nonterminal < numNonterminals; ++nonterminal) {
        std::vector<NodeId> alternatives;
        alternatives.reserve(
          checkedProduct(
            binaryRules[nonterminal].size(),
            static_cast<std::size_t>(length - 1),
            "too many DP alternatives"));

        for (const DenseBinaryRule& rule : binaryRules[nonterminal]) {
          for (int split = 1; split < length; ++split) {
            const NodeId left = dp[dpIndex(rule.left, split, offset, n)];
            const NodeId right = dp[dpIndex(rule.right, length - split, offset + split, n)];
            const NodeId product = makeMult(program, left, right);
            if (product != EMPTY_NODE) {
              alternatives.push_back(product);
            }
          }
        }

        dp[dpIndex(nonterminal, length, offset, n)] = makeAdd(program, alternatives);
      }
    }
  }

  const std::size_t start = lookupDense(nonterminalIndex, cfg.start, "nonterminal");
  program.root = dp[dpIndex(start, n, 0, n)];

  makeMinimal(program);

  return program;

}
