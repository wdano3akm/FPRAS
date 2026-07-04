#include <stdint.h>
#include <vector>

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

PlusTimesProgram compileCFGToPlusTimes(const CFG& cfg, int n);
