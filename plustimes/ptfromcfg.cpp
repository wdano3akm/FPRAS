#include "plustimes.hpp"
#include "../cfg/CFG.hpp"

#include <limits>
#include <stdexcept>
#include <vector>

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

  program = pruneReachableTopological(program);
  if (program.root != EMPTY_NODE) {
    program = depthReduceVSBR(program);
  }

  return program;

}
