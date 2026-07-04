#include "plustimes.hpp"
#include "../cfg/CFG.hpp"

#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct DenseBinaryRule {
  std::size_t left;
  std::size_t right;
};

NodeId appendNode(PlusTimesProgram& program, PTNode node)
{
  if (program.nodes.size() >= EMPTY_NODE) {
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
  for (NodeId child : rawChildren) {
    if (child == EMPTY_NODE) {
      continue;
    }

    requireNode(program, child);
    const int childDegree = program.nodes[child].degree;
    if (degree == -1) {
      degree = childDegree;
    } else if (degree != childDegree) {
      throw std::logic_error("non-homogeneous plus-times addition");
    }

    children.push_back(child);
  }

  if (children.empty()) {
    return EMPTY_NODE;
  }
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

} // namespace

PlusTimesProgram compileCFGToPlusTimes(CFG &cfg, int n)
{
  if (n < 0) {
    throw std::invalid_argument("n cannot be negative");
  }
  if (n == 0) {
    throw std::runtime_error("length-zero CFG compilation needs constant nodes");
  }

  const auto nonterminalIndex = makeDenseIndex(cfg.nonterminalIds, "nonterminal");
  const auto terminalIndex = makeDenseIndex(cfg.terminalIds, "terminal");
  const std::size_t numNonterminals = cfg.nonterminalIds.size();
  const std::size_t numTerminals = cfg.terminalIds.size();

  if (numNonterminals == 0) {
    throw std::invalid_argument("CFG must have at least one nonterminal");
  }

  const std::size_t numVars = static_cast<std::size_t>(n) * numTerminals;
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

  const std::size_t dpSize =
    numNonterminals * (static_cast<std::size_t>(n) + 1) * (static_cast<std::size_t>(n) + 1);
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
        std::vector<NodeId> ruleChildren;
        ruleChildren.reserve(binaryRules[nonterminal].size());

        for (const DenseBinaryRule& rule : binaryRules[nonterminal]) {
          std::vector<NodeId> splitChildren;
          splitChildren.reserve(static_cast<std::size_t>(length - 1));

          for (int split = 1; split < length; ++split) {
            const NodeId left = dp[dpIndex(rule.left, split, offset, n)];
            const NodeId right = dp[dpIndex(rule.right, length - split, offset + split, n)];
            splitChildren.push_back(makeMult(program, left, right));
          }

          ruleChildren.push_back(makeAdd(program, splitChildren));
        }

        dp[dpIndex(nonterminal, length, offset, n)] = makeAdd(program, ruleChildren);
      }
    }
  }

  const std::size_t start = lookupDense(nonterminalIndex, cfg.start, "nonterminal");
  program.root = dp[dpIndex(start, n, 0, n)];

  return program;
}
