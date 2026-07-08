#include "alg_checks.hpp"

#include <algorithm>
#include <stdexcept>


void requireReady(bool condition, const char* message)
{
    if (!condition) {
        throw std::invalid_argument(message);
    }
}

namespace {

int ceilLog2Positive(int value)
{
    requireReady(value > 0, "degree must be positive");

    int x = value - 1;
    int result = 0;
    while (x > 0) {
        x >>= 1;
        ++result;
    }
    return result;
}

} // namespace

// Description: 
//  Iterates through the plus,times tree and computes 
//  the degree of all nodes inside it.
//  - Variable nodes have degree 1 
//  - Addition nodes, same degree as children (all same deg)
//  - Multiplication nodes, degree is sum of the multiplied monomials degree
//  Asserts that the root is the same degree as precomputed program deg
// Input: 
//  PlusTimesProgram &P 
// Output: 
//  vector<int> the array with degrees
std::vector<int> computeDegreesOrThrow(const PlusTimesProgram& P)
{
    requireReady(P.root != EMPTY_NODE, "empty plus-times program has no FPRAS root");
    requireReady(!P.nodes.empty(), "plus-times program has no nodes");
    requireReady(P.root < P.nodes.size(), "invalid plus-times root");

    std::vector<int> degrees(P.nodes.size(), 0);

    for (std::size_t i = 0; i < P.nodes.size(); ++i) {
        const PTNode& node = P.nodes[i];

        if (node.kind == PTKind::Var) {
            requireReady(node.children.empty(), "variable plus-times node has children");
            requireReady(P.numVars > 0 && node.var < P.numVars, "plus-times variable id is out of range");
            requireReady(node.degree == 1, "variable plus-times node must have degree 1");
            degrees[i] = 1;
        } else if (node.kind == PTKind::Add) {
            requireReady(!node.children.empty(), "addition plus-times node has no children");

            int degree = -1;
            for (NodeId child : node.children) {
                requireReady(child < i, "plus-times nodes are not topologically ordered");
                if (degree == -1) {
                    degree = degrees[child];
                } else {
                    requireReady(degrees[child] == degree, "non-homogeneous plus-times addition");
                }
            }

            requireReady(node.degree == degree, "stored plus-times addition degree is incorrect");
            degrees[i] = degree;
        } else if (node.kind == PTKind::Mul) {
            requireReady(node.children.size() == 2, "multiplication plus-times node must have fan-in exactly 2");

            const NodeId left = node.children[0];
            const NodeId right = node.children[1];
            requireReady(left < i && right < i, "plus-times nodes are not topologically ordered");

            const int degree = degrees[left] + degrees[right];
            requireReady(node.degree == degree, "stored plus-times multiplication degree is incorrect");
            degrees[i] = degree;
        } else {
            throw std::invalid_argument("unknown plus-times node kind");
        }
    }

    requireReady(P.degree == degrees[P.root], "program degree does not match root degree");
    return degrees;
}

// Description:
//  Ensures that each child has higher NodeId than parent 
// Input: 
//  PlusTimesProgram &P
// Output:
//  bool is tree topologically sorted?T/F
bool isTopologicallyOrdered(const PlusTimesProgram& P)
{
    if (P.root == EMPTY_NODE || P.root >= P.nodes.size()) {
        return false;
    }

    for (std::size_t i = 0; i < P.nodes.size(); ++i) {
        for (NodeId child : P.nodes[i].children) {
            if (child >= i) {
                return false;
            }
        }
    }
    return true;
}
// Description: 
//  Checks that the plus,times tree does not contain 
//  Nodes that are not Var, Add or Mul and that every Add node has
//  at least one child 
// Input:
//  PlusTimesProgram &P 
// Output:
//  bool is tree compliant?T/F
bool noConstants(const PlusTimesProgram& P)
{
    for (const PTNode& node : P.nodes) {
        if (node.kind != PTKind::Var && node.kind != PTKind::Add && node.kind != PTKind::Mul) {
            return false;
        }
        if (node.kind == PTKind::Add && node.children.empty()) {
            return false;
        }
    }
    return true;
}

// Description:
//  Checks that all Mul nodes have exactly 2 children
// Input:
//  PlusTimesProgram &P
// Output: 
//  bool is tree compliant?T/F
bool allProductsBinary(const PlusTimesProgram& P)
{
    for (const PTNode& node : P.nodes) {
        if (node.kind == PTKind::Mul && node.children.size() != 2) {
            return false;
        }
    }
    return true;
}

// Description: 
//  Check that all nodes of type Sum have been normalized.
//  I.e. No node of type sum should have a children of type sum
// Input: 
//  PlusTimesProgram &P 
// Output: 
//  bool are all sum nodes compliant?T/F
bool noSumHasSumChild(const PlusTimesProgram& P)
{
    for (const PTNode& node : P.nodes) {
        if (node.kind != PTKind::Add) {
            continue;
        }
        for (NodeId child : node.children) {
            if (child >= P.nodes.size() || P.nodes[child].kind == PTKind::Add) {
                return false;
            }
        }
    }
    return true;
}

bool isHomogeneous(const PlusTimesProgram& P, const std::vector<int>& degrees)
{
    try {
        return degrees == computeDegreesOrThrow(P);
    } catch (const std::invalid_argument&) {
        return false;
    }
}

// Description: 
//  Checks that two unordered_set are disjoint.
// Input:
//  unordered_set<VarId> &a
//  unordered_set<VarId> &b
// Output:
//  bool true if disjoint
bool disjointVarSets(const std::unordered_set<VarId>& left, const std::unordered_set<VarId>& right)
{
    const auto& small = left.size() <= right.size() ? left : right;
    const auto& large = left.size() <= right.size() ? right : left;
    for (VarId var : small) {
        if (large.find(var) != large.end()) {
            return false;
        }
    }
    return true;
}

// Description:
//  Asserts that right and left subprogram of multiplication 
//  node reference disjoint subsets
// Input:
//  PlusTimesProgram &P
// Output: 
//  bool are all the set disjoint?T/F
bool isStructurallyMultilinear(const PlusTimesProgram& P)
{
    VariableSets variableSets(P.nodes.size());

    for (std::size_t i = 0; i < P.nodes.size(); ++i) {
        const PTNode& node = P.nodes[i];
        if (node.kind == PTKind::Var) {
            variableSets[i].insert(static_cast<VarId>(node.var));
        } else if (node.kind == PTKind::Add) {
            for (NodeId child : node.children) {
                variableSets[i].insert(variableSets[child].begin(), variableSets[child].end());
            }
        } else if (node.kind == PTKind::Mul) {
            const NodeId left = node.children[0];
            const NodeId right = node.children[1];
            if (!disjointVarSets(variableSets[left], variableSets[right])) {
                return false;
            }
            variableSets[i] = variableSets[left];
            variableSets[i].insert(variableSets[right].begin(), variableSets[right].end());
        }
    }

    return true;
}

// Description:
//  Computes overrall height of tree.
//  In addition, check that the node's height is 0 if variable node.
//  Otherwise, computes height as highest child plus one
// Input:
//  PlusTimesProgram &P
// Return: 
//  int height of root
int rootHeightOrThrow(const PlusTimesProgram& P)
{
    requireReady(P.root != EMPTY_NODE, "empty plus-times program has no FPRAS root");
    requireReady(P.root < P.nodes.size(), "invalid plus-times root");

    std::vector<int> heights(P.nodes.size(), 0);
    for (std::size_t i = 0; i < P.nodes.size(); ++i) {
        const PTNode& node = P.nodes[i];
        if (node.kind == PTKind::Var) {
            requireReady(node.children.empty(), "variable plus-times node has children");
            heights[i] = 0;
            continue;
        }

        int height = 0;
        for (NodeId child : node.children) {
            requireReady(child < i, "plus-times nodes are not topologically ordered");
            height = std::max(height, heights[child] + 1);
        }
        heights[i] = height;
    }

    return heights[P.root];
}

// Description:
//  Performs the following checks:
//  - toposort of tree 
//  - only var/mul/add nodes 
//  - all mul have binary children 
//  - no sum has sum child
//  - homogeneous tree 
//  - multilinear 
//  - height <= 3 * ceil(log(rootDeg))
// Input: 
//  PlusTimesProgram &P
// Output:
//  void
void assertReadyForFPRAS(const PlusTimesProgram& P)
{
    const std::vector<int> degrees = computeDegreesOrThrow(P);
    const int rootDegree = degrees[P.root];

    requireReady(isTopologicallyOrdered(P), "plus-times nodes are not topologically ordered");
    requireReady(noConstants(P), "plus-times program contains a constant-like node");
    requireReady(allProductsBinary(P), "multiplication plus-times node must have fan-in exactly 2");
    requireReady(noSumHasSumChild(P), "addition plus-times nodes must be flattened");
    requireReady(isHomogeneous(P, degrees), "plus-times program is not homogeneous");
    requireReady(isStructurallyMultilinear(P), "plus-times program is not structurally multilinear");

    const int allowedHeight = rootDegree <= 1 ? 1 : 3 * ceilLog2Positive(rootDegree);
    requireReady(rootHeightOrThrow(P) <= allowedHeight, "plus-times root exceeds the FPRAS height bound");
}

