#ifndef DNNF_HPP
#define DNNF_HPP

#include <cstdint>
#include <limits>
#include <vector>

using DNNFNodeId = uint32_t;
using NodeId = DNNFNodeId;
using DNNFLiteral = int32_t;

static constexpr DNNFNodeId DNNF_EMPTY_NODE = std::numeric_limits<DNNFNodeId>::max();

enum DNNFKind {
    Or,
    And,
    TrueVar,
    FalseVar,
};

// d4 stores propagated literals on an arc from a parent to a child.
struct DNNFArc {
    DNNFNodeId target = DNNF_EMPTY_NODE;
    std::vector<DNNFLiteral> literals;
};

struct DNNFNode {
    DNNFKind kind;
    uint32_t degree = 0;
    DNNFNodeId id = 0;
    std::vector<DNNFArc> children;
};

struct DNNF {
    std::vector<DNNFNode> nodes;
    DNNFNodeId root = DNNF_EMPTY_NODE;
    uint32_t numVars = 0;
};

#endif
