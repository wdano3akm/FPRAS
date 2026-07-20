#pragma once

#include <stdint.h>
#include <vector>
#include <unordered_map>

struct CFG;

using NodeId = uint32_t;
using PolyVarId = uint32_t;

static constexpr NodeId EMPTY_NODE = UINT32_MAX;

enum class PTKind {
    Var,
    Add,
    Mul
};

struct PTNode {
    PTKind kind;
    int degree = 0;

    // Var
    PolyVarId var = 0;

    // Add: arbitrary fan-in
    // Mul: exactly two children
    std::vector<NodeId> children;
};

struct PlusTimesProgram {
    std::vector<PTNode> nodes;
    NodeId root = EMPTY_NODE;
    uint32_t numVars = 0;
    int degree = 0;
};

struct DenseBinaryRule {
  std::size_t left;
  std::size_t right;
};

PlusTimesProgram compileCFGToPlusTimes(const CFG& cfg, int n);
void validateDisjointNonzeroIds(const CFG& cfg);
std::unordered_map<uint32_t, std::size_t> makeDenseIndex(
  const std::vector<uint32_t>& ids,
  const char* kind);
std::size_t checkedProduct(std::size_t left, std::size_t right, const char* message);
std::size_t lookupDense(
  const std::unordered_map<uint32_t, std::size_t>& dense,
  uint32_t id,
  const char* kind);
std::size_t variableIndex(std::size_t terminal, int offset, std::size_t numTerminals);
NodeId makeVar(PlusTimesProgram& program, PolyVarId var);
std::size_t dpIndex(std::size_t nonterminal, int length, int offset, int n);
NodeId makeAdd(PlusTimesProgram& program, const std::vector<NodeId>& rawChildren);
NodeId makeMult(PlusTimesProgram& program, NodeId left, NodeId right);
PlusTimesProgram pruneReachableTopological(const PlusTimesProgram& input);
PlusTimesProgram depthReduceVSBR(const PlusTimesProgram& input);


