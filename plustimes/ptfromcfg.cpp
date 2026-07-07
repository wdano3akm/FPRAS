#include "plustimes.hpp"
#include "../cfg/CFG.hpp"

#include <algorithm>
#include <functional>
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

PlusTimesProgram pruneReachableTopological(const PlusTimesProgram& input)
{
  PlusTimesProgram out;
  out.numVars = input.numVars;
  out.degree = input.degree;

  if (input.root == EMPTY_NODE) {
    return out;
  }
  if (input.root >= input.nodes.size()) {
    throw std::logic_error("invalid plus-times root");
  }

  std::vector<char> seen(input.nodes.size(), false);
  std::vector<NodeId> order;

  std::function<void(NodeId)> dfs = [&](NodeId node) {
    if (node == EMPTY_NODE) {
      return;
    }
    if (node >= input.nodes.size()) {
      throw std::logic_error("invalid plus-times node id");
    }
    if (seen[node]) {
      return;
    }

    seen[node] = true;
    for (NodeId child : input.nodes[node].children) {
      dfs(child);
    }
    order.push_back(node);
  };

  dfs(input.root);

  std::vector<NodeId> oldToNew(input.nodes.size(), EMPTY_NODE);
  out.nodes.reserve(order.size());

  for (NodeId old : order) {
    PTNode node = input.nodes[old];
    for (NodeId& child : node.children) {
      if (child == EMPTY_NODE || child >= oldToNew.size() || oldToNew[child] == EMPTY_NODE) {
        throw std::logic_error("invalid plus-times child during pruning");
      }
      child = oldToNew[child];
    }

    oldToNew[old] = appendNode(out, std::move(node));
  }

  out.root = oldToNew[input.root];
  if (out.root != EMPTY_NODE) {
    out.degree = out.nodes[out.root].degree;
  }

  return out;
}

enum class ExprKind {
  Unset,
  Zero,
  One,
  Node
};

struct Expr {
  ExprKind kind = ExprKind::Unset;
  NodeId id = 0;

  static Expr unset() { return {ExprKind::Unset, 0}; }
  static Expr zero() { return {ExprKind::Zero, 0}; }
  static Expr one() { return {ExprKind::One, 0}; }
  static Expr node(NodeId id) { return {ExprKind::Node, id}; }
};

bool isSet(Expr expr)
{
  return expr.kind != ExprKind::Unset;
}

bool isZero(Expr expr)
{
  return expr.kind == ExprKind::Zero;
}

struct NodeKey {
  PTKind kind = PTKind::Var;
  PolyVarId var = 0;
  std::vector<NodeId> children;

  bool operator==(const NodeKey& other) const
  {
    return kind == other.kind && var == other.var && children == other.children;
  }
};

std::size_t hashCombine(std::size_t seed, std::size_t value)
{
  return seed ^ (value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2));
}

struct NodeKeyHash {
  std::size_t operator()(const NodeKey& key) const
  {
    std::size_t seed = std::hash<int>{}(static_cast<int>(key.kind));
    seed = hashCombine(seed, std::hash<PolyVarId>{}(key.var));
    for (NodeId child : key.children) {
      seed = hashCombine(seed, std::hash<NodeId>{}(child));
    }
    return seed;
  }
};

bool disjointSets(
  const std::unordered_set<PolyVarId>& left,
  const std::unordered_set<PolyVarId>& right)
{
  const auto& small = left.size() <= right.size() ? left : right;
  const auto& large = left.size() <= right.size() ? right : left;

  for (PolyVarId var : small) {
    if (large.find(var) != large.end()) {
      return false;
    }
  }
  return true;
}

void appendVarsetUnion(
  std::unordered_set<PolyVarId>& out,
  const std::unordered_set<PolyVarId>& in)
{
  out.insert(in.begin(), in.end());
}

PlusTimesProgram normalizeForDepthReduction(const PlusTimesProgram& input)
{
  PlusTimesProgram normalized;
  normalized.numVars = input.numVars;
  normalized.degree = input.degree;

  if (input.root == EMPTY_NODE) {
    return normalized;
  }
  if (input.root >= input.nodes.size()) {
    throw std::logic_error("invalid plus-times root");
  }

  enum class Mark {
    Unseen,
    Visiting,
    Done
  };

  std::vector<Mark> marks(input.nodes.size(), Mark::Unseen);
  std::vector<NodeId> order;

  std::function<void(NodeId)> dfs = [&](NodeId node) {
    if (node == EMPTY_NODE) {
      return;
    }
    if (node >= input.nodes.size()) {
      throw std::logic_error("invalid plus-times node id");
    }
    if (marks[node] == Mark::Visiting) {
      throw std::logic_error("cyclic plus-times program");
    }
    if (marks[node] == Mark::Done) {
      return;
    }

    marks[node] = Mark::Visiting;
    const PTNode& ptNode = input.nodes[node];
    if (ptNode.kind == PTKind::Var) {
      if (!ptNode.children.empty()) {
        throw std::logic_error("variable plus-times node has children");
      }
    } else if (ptNode.kind == PTKind::Mul) {
      if (ptNode.children.size() != 2) {
        throw std::logic_error("multiplication plus-times node must have two children");
      }
    }

    for (NodeId child : ptNode.children) {
      if (child == EMPTY_NODE) {
        if (ptNode.kind == PTKind::Mul) {
          throw std::logic_error("multiplication plus-times node has an empty child");
        }
        continue;
      }
      dfs(child);
    }

    marks[node] = Mark::Done;
    order.push_back(node);
  };

  dfs(input.root);

  std::vector<NodeId> oldToNew(input.nodes.size(), EMPTY_NODE);
  std::vector<std::unordered_set<PolyVarId>> varsets;
  normalized.nodes.reserve(order.size());
  varsets.reserve(order.size());

  const auto addNode =
    [&](PTNode node, std::unordered_set<PolyVarId> varset) -> NodeId {
      const NodeId id = appendNode(normalized, std::move(node));
      varsets.push_back(std::move(varset));
      return id;
    };

  for (NodeId old : order) {
    const PTNode& oldNode = input.nodes[old];

    if (oldNode.kind == PTKind::Var) {
      if (input.numVars == 0 || oldNode.var >= input.numVars) {
        throw std::logic_error("plus-times variable id is out of range");
      }

      PTNode node;
      node.kind = PTKind::Var;
      node.degree = 1;
      node.var = oldNode.var;

      oldToNew[old] = addNode(std::move(node), {oldNode.var});
      continue;
    }

    if (oldNode.kind == PTKind::Add) {
      std::vector<NodeId> children;
      std::unordered_set<PolyVarId> varset;
      int degree = -1;

      const auto addChild = [&](NodeId child) {
        if (child == EMPTY_NODE) {
          return;
        }
        if (child >= normalized.nodes.size()) {
          throw std::logic_error("invalid normalized plus-times child");
        }

        const int childDegree = normalized.nodes[child].degree;
        if (degree == -1) {
          degree = childDegree;
        } else if (degree != childDegree) {
          throw std::logic_error("non-homogeneous plus-times addition");
        }

        children.push_back(child);
        appendVarsetUnion(varset, varsets[child]);
      };

      for (NodeId oldChild : oldNode.children) {
        if (oldChild == EMPTY_NODE) {
          continue;
        }

        const NodeId child = oldToNew[oldChild];
        if (child == EMPTY_NODE) {
          continue;
        }

        if (normalized.nodes[child].kind == PTKind::Add) {
          for (NodeId grandchild : normalized.nodes[child].children) {
            addChild(grandchild);
          }
        } else {
          addChild(child);
        }
      }

      if (children.empty()) {
        oldToNew[old] = EMPTY_NODE;
      } else if (children.size() == 1) {
        oldToNew[old] = children.front();
      } else {
        PTNode node;
        node.kind = PTKind::Add;
        node.degree = degree;
        node.children = std::move(children);
        oldToNew[old] = addNode(std::move(node), std::move(varset));
      }
      continue;
    }

    if (oldNode.kind != PTKind::Mul) {
      throw std::logic_error("unknown plus-times node kind");
    }

    NodeId left = oldToNew[oldNode.children[0]];
    NodeId right = oldToNew[oldNode.children[1]];
    if (left == EMPTY_NODE || right == EMPTY_NODE) {
      oldToNew[old] = EMPTY_NODE;
      continue;
    }

    if (normalized.nodes[left].degree < normalized.nodes[right].degree) {
      std::swap(left, right);
    }

    if (!disjointSets(varsets[left], varsets[right])) {
      throw std::logic_error("non-multilinear plus-times multiplication");
    }

    std::unordered_set<PolyVarId> varset = varsets[left];
    appendVarsetUnion(varset, varsets[right]);

    PTNode node;
    node.kind = PTKind::Mul;
    node.degree = normalized.nodes[left].degree + normalized.nodes[right].degree;
    node.children = {left, right};
    oldToNew[old] = addNode(std::move(node), std::move(varset));
  }

  normalized.root = oldToNew[input.root];
  if (normalized.root == EMPTY_NODE) {
    normalized.nodes.clear();
    return normalized;
  }

  if (input.degree > 0 && input.degree != normalized.nodes[normalized.root].degree) {
    throw std::logic_error("plus-times program degree does not match root degree");
  }
  normalized.degree = normalized.nodes[normalized.root].degree;

  return pruneReachableTopological(normalized);
}

uint32_t ceilLog2(uint32_t value)
{
  if (value <= 1) {
    return 0;
  }

  --value;
  uint32_t result = 0;
  while (value > 0) {
    value >>= 1;
    ++result;
  }
  return result;
}

int programHeight(const PlusTimesProgram& program)
{
  if (program.root == EMPTY_NODE) {
    return 0;
  }

  std::vector<int> height(program.nodes.size(), 0);
  for (std::size_t i = 0; i < program.nodes.size(); ++i) {
    const PTNode& node = program.nodes[i];
    if (node.kind == PTKind::Var) {
      height[i] = 0;
      continue;
    }

    int nodeHeight = 0;
    for (NodeId child : node.children) {
      if (child >= i) {
        throw std::logic_error("plus-times program is not topologically ordered");
      }
      nodeHeight = std::max(nodeHeight, height[child] + 1);
    }
    height[i] = nodeHeight;
  }

  return height[program.root];
}

class DepthReducer {
public:
  explicit DepthReducer(const PlusTimesProgram& input)
    : oldP(input),
      deg(input.nodes.size(), 0),
      value(input.nodes.size(), Expr::unset()),
      context(input.nodes.size(), std::vector<Expr>(input.nodes.size(), Expr::unset()))
  {
    for (std::size_t i = 0; i < oldP.nodes.size(); ++i) {
      if (oldP.nodes[i].degree <= 0) {
        throw std::logic_error("depth reduction requires positive-degree nodes");
      }
      deg[i] = static_cast<uint32_t>(oldP.nodes[i].degree);
    }

    newP.numVars = oldP.numVars;
    newP.degree = oldP.degree;
  }

  PlusTimesProgram run()
  {
    initDegreeOne();

    const uint32_t rootDegree = deg[oldP.root];
    for (uint32_t low = 1; low < rootDegree;) {
      const uint32_t high =
        low > std::numeric_limits<uint32_t>::max() / 2
          ? std::numeric_limits<uint32_t>::max()
          : low * 2;
      buildLayer(low, high);
      if (high == std::numeric_limits<uint32_t>::max()) {
        break;
      }
      low = high;
    }

    Expr root = value[oldP.root];
    if (!isSet(root)) {
      throw std::logic_error("depth reduction did not construct the root value");
    }
    return finalize(root);
  }

private:
  const PlusTimesProgram& oldP;
  PlusTimesProgram newP;
  std::vector<uint32_t> deg;
  std::vector<Expr> value;
  std::vector<std::vector<Expr>> context;
  std::unordered_map<NodeKey, NodeId, NodeKeyHash> cache;

  bool isProduct(NodeId node) const
  {
    return oldP.nodes[node].kind == PTKind::Mul;
  }

  NodeId left(NodeId node) const
  {
    return oldP.nodes[node].children[0];
  }

  NodeId right(NodeId node) const
  {
    return oldP.nodes[node].children[1];
  }

  Expr makeVar(PolyVarId var)
  {
    NodeKey key;
    key.kind = PTKind::Var;
    key.var = var;

    const auto cached = cache.find(key);
    if (cached != cache.end()) {
      return Expr::node(cached->second);
    }

    PTNode node;
    node.kind = PTKind::Var;
    node.degree = 1;
    node.var = var;

    const NodeId id = appendNode(newP, std::move(node));
    cache.emplace(std::move(key), id);
    return Expr::node(id);
  }

  NodeId internNode(NodeKey key, PTNode node)
  {
    const auto cached = cache.find(key);
    if (cached != cache.end()) {
      return cached->second;
    }

    const NodeId id = appendNode(newP, std::move(node));
    cache.emplace(std::move(key), id);
    return id;
  }

  Expr makeSum(std::vector<Expr> terms)
  {
    std::vector<NodeId> ids;
    bool sawOne = false;

    for (Expr term : terms) {
      if (term.kind == ExprKind::Unset) {
        throw std::logic_error("unset expression in depth-reduction sum");
      }
      if (term.kind == ExprKind::Zero) {
        continue;
      }
      if (term.kind == ExprKind::One) {
        sawOne = true;
        continue;
      }

      if (newP.nodes[term.id].kind == PTKind::Add) {
        ids.insert(
          ids.end(),
          newP.nodes[term.id].children.begin(),
          newP.nodes[term.id].children.end());
      } else {
        ids.push_back(term.id);
      }
    }

    if (ids.empty()) {
      return sawOne ? Expr::one() : Expr::zero();
    }
    if (sawOne) {
      throw std::logic_error("inhomogeneous constant term in depth-reduction sum");
    }

    std::sort(ids.begin(), ids.end());

    if (ids.size() == 1) {
      return Expr::node(ids.front());
    }

    int degree = newP.nodes[ids.front()].degree;
    for (NodeId child : ids) {
      if (newP.nodes[child].degree != degree) {
        throw std::logic_error("non-homogeneous depth-reduction sum");
      }
    }

    NodeKey key;
    key.kind = PTKind::Add;
    key.children = ids;

    PTNode node;
    node.kind = PTKind::Add;
    node.degree = degree;
    node.children = std::move(ids);

    return Expr::node(internNode(std::move(key), std::move(node)));
  }

  Expr makeProduct(Expr a, Expr b)
  {
    if (a.kind == ExprKind::Unset || b.kind == ExprKind::Unset) {
      throw std::logic_error("unset expression in depth-reduction product");
    }
    if (a.kind == ExprKind::Zero || b.kind == ExprKind::Zero) {
      return Expr::zero();
    }
    if (a.kind == ExprKind::One) {
      return b;
    }
    if (b.kind == ExprKind::One) {
      return a;
    }

    NodeId leftChild = a.id;
    NodeId rightChild = b.id;
    if (newP.nodes[leftChild].degree < newP.nodes[rightChild].degree) {
      std::swap(leftChild, rightChild);
    }

    NodeKey key;
    key.kind = PTKind::Mul;
    key.children = {leftChild, rightChild};

    PTNode node;
    node.kind = PTKind::Mul;
    node.degree = newP.nodes[leftChild].degree + newP.nodes[rightChild].degree;
    node.children = {leftChild, rightChild};

    return Expr::node(internNode(std::move(key), std::move(node)));
  }

  Expr requireValue(NodeId node) const
  {
    if (!isSet(value[node])) {
      throw std::logic_error("missing depth-reduction value");
    }
    return value[node];
  }

  Expr requireContext(NodeId v, NodeId w) const
  {
    if (v == w) {
      return Expr::one();
    }
    if (deg[w] < deg[v]) {
      return Expr::zero();
    }
    if (!isSet(context[v][w])) {
      throw std::logic_error("missing depth-reduction context");
    }
    return context[v][w];
  }

  Expr buildDegreeOneValue(NodeId node)
  {
    if (isSet(value[node])) {
      return value[node];
    }
    if (deg[node] != 1) {
      throw std::logic_error("attempted degree-one construction for higher-degree node");
    }

    const PTNode& ptNode = oldP.nodes[node];
    if (ptNode.kind == PTKind::Var) {
      value[node] = makeVar(ptNode.var);
      return value[node];
    }
    if (ptNode.kind == PTKind::Add) {
      std::vector<Expr> terms;
      terms.reserve(ptNode.children.size());
      for (NodeId child : ptNode.children) {
        terms.push_back(buildDegreeOneValue(child));
      }
      value[node] = makeSum(std::move(terms));
      return value[node];
    }

    throw std::logic_error("degree-one multiplication node is impossible");
  }

  Expr buildContextDiffOne(NodeId v, NodeId w)
  {
    if (v == w) {
      return Expr::one();
    }
    if (deg[w] < deg[v]) {
      return Expr::zero();
    }
    if (deg[w] - deg[v] > 1) {
      return Expr::zero();
    }
    if (isSet(context[v][w])) {
      return context[v][w];
    }

    const PTNode& wNode = oldP.nodes[w];
    if (wNode.kind == PTKind::Var) {
      context[v][w] = Expr::zero();
      return context[v][w];
    }
    if (wNode.kind == PTKind::Add) {
      std::vector<Expr> terms;
      terms.reserve(wNode.children.size());
      for (NodeId child : wNode.children) {
        terms.push_back(buildContextDiffOne(v, child));
      }
      context[v][w] = makeSum(std::move(terms));
      return context[v][w];
    }

    if (deg[w] == deg[v]) {
      context[v][w] = Expr::zero();
      return context[v][w];
    }

    const NodeId heavy = left(w);
    const NodeId light = right(w);
    Expr inner = buildContextDiffOne(v, heavy);
    if (isZero(inner)) {
      context[v][w] = Expr::zero();
    } else {
      context[v][w] = makeProduct(requireValue(light), inner);
    }
    return context[v][w];
  }

  void initDegreeOne()
  {
    for (std::size_t i = 0; i < oldP.nodes.size(); ++i) {
      if (deg[i] == 1) {
        buildDegreeOneValue(static_cast<NodeId>(i));
      }
    }

    for (std::size_t v = 0; v < oldP.nodes.size(); ++v) {
      context[v][v] = Expr::one();
    }

    for (std::size_t v = 0; v < oldP.nodes.size(); ++v) {
      for (std::size_t w = 0; w < oldP.nodes.size(); ++w) {
        if (v != w && deg[w] >= deg[v] && deg[w] - deg[v] <= 1) {
          buildContextDiffOne(static_cast<NodeId>(v), static_cast<NodeId>(w));
        }
      }
    }
  }

  std::vector<NodeId> crossingGates(uint32_t threshold) const
  {
    std::vector<NodeId> out;
    for (std::size_t i = 0; i < oldP.nodes.size(); ++i) {
      const NodeId node = static_cast<NodeId>(i);
      if (!isProduct(node)) {
        continue;
      }
      const NodeId heavy = left(node);
      if (deg[node] > threshold && deg[heavy] <= threshold) {
        out.push_back(node);
      }
    }
    return out;
  }

  Expr valueOfCrossingGate(NodeId node)
  {
    return makeProduct(requireValue(left(node)), requireValue(right(node)));
  }

  Expr buildValueFromCuts(NodeId w, uint32_t threshold)
  {
    std::vector<Expr> terms;

    for (NodeId cut : crossingGates(threshold)) {
      Expr ctx = requireContext(cut, w);
      if (isZero(ctx)) {
        continue;
      }

      Expr term = makeProduct(valueOfCrossingGate(cut), ctx);
      if (!isZero(term)) {
        terms.push_back(term);
      }
    }

    return makeSum(std::move(terms));
  }

  Expr buildContextFromCuts(NodeId v, NodeId w, uint32_t threshold)
  {
    if (v == w) {
      return Expr::one();
    }
    if (deg[w] < deg[v]) {
      return Expr::zero();
    }

    const uint32_t shiftedThreshold = deg[v] + threshold;
    std::vector<Expr> terms;

    for (NodeId cut : crossingGates(shiftedThreshold)) {
      Expr outer = requireContext(cut, w);
      if (isZero(outer)) {
        continue;
      }

      const NodeId heavy = left(cut);
      const NodeId light = right(cut);
      Expr inner = requireContext(v, heavy);
      if (isZero(inner)) {
        continue;
      }

      Expr term = makeProduct(makeProduct(requireValue(light), inner), outer);
      if (!isZero(term)) {
        terms.push_back(term);
      }
    }

    return makeSum(std::move(terms));
  }

  void buildLayer(uint32_t low, uint32_t high)
  {
    for (std::size_t w = 0; w < oldP.nodes.size(); ++w) {
      if (low < deg[w] && deg[w] <= high) {
        value[w] = buildValueFromCuts(static_cast<NodeId>(w), low);
      }
    }

    for (std::size_t v = 0; v < oldP.nodes.size(); ++v) {
      for (std::size_t w = 0; w < oldP.nodes.size(); ++w) {
        if (deg[w] <= deg[v]) {
          continue;
        }

        const uint32_t diff = deg[w] - deg[v];
        if (low < diff && diff <= high) {
          context[v][w] =
            buildContextFromCuts(static_cast<NodeId>(v), static_cast<NodeId>(w), low);
        }
      }
    }
  }

  PlusTimesProgram finalize(Expr root)
  {
    if (root.kind == ExprKind::Zero) {
      PlusTimesProgram empty;
      empty.numVars = oldP.numVars;
      empty.degree = oldP.degree;
      return empty;
    }
    if (root.kind != ExprKind::Node) {
      throw std::logic_error("depth reduction produced a non-node root");
    }

    newP.root = root.id;
    newP.degree = oldP.degree;
    PlusTimesProgram out = normalizeForDepthReduction(newP);

    const uint32_t degree = static_cast<uint32_t>(out.degree);
    const int allowedHeight =
      degree <= 1 ? 1 : static_cast<int>(3 * ceilLog2(degree));
    if (programHeight(out) > allowedHeight) {
      throw std::logic_error("depth reduction did not meet the height bound");
    }

    return out;
  }
};

PlusTimesProgram depthReduceVSBR(const PlusTimesProgram& input)
{
  PlusTimesProgram normalized = normalizeForDepthReduction(input);
  if (normalized.root == EMPTY_NODE) {
    return normalized;
  }

  DepthReducer reducer(normalized);
  return reducer.run();
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
