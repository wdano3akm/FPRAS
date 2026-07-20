#include "DNNFParser.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

struct ParsedLine {
    std::size_t number;
    std::vector<std::string> tokens;
};

struct ParsedArc {
    std::size_t line;
    DNNFNodeId source;
    DNNFArc arc;
};

bool isCommentOrBlank(const std::string& line)
{
    const auto first = line.find_first_not_of(" \t\r\n");
    return first == std::string::npos || line[first] == 'c';
}

std::vector<std::string> tokenize(const std::string& line)
{
    std::istringstream stream(line);
    std::vector<std::string> tokens;
    std::string token;
    while (stream >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<ParsedLine> readMeaningfulLines(std::istream& input)
{
    std::vector<ParsedLine> lines;
    std::string line;
    std::size_t lineNumber = 0;

    while (std::getline(input, line)) {
        ++lineNumber;
        if (isCommentOrBlank(line)) {
            continue;
        }

        lines.push_back({lineNumber, tokenize(line)});
    }

    return lines;
}

std::runtime_error parseError(std::size_t line, const std::string& message)
{
    return std::runtime_error("DNNF parse error on line " + std::to_string(line) + ": " + message);
}

uint32_t parseUint32(const std::string& token, std::size_t line)
{
    if (token.empty()) {
        throw parseError(line, "empty numeric token");
    }

    std::size_t parsed = 0;
    unsigned long value = 0;
    try {
        value = std::stoul(token, &parsed, 10);
    } catch (const std::exception&) {
        throw parseError(line, "expected uint32 token, found '" + token + "'");
    }

    if (parsed != token.size() || value > std::numeric_limits<uint32_t>::max()) {
        throw parseError(line, "expected uint32 token, found '" + token + "'");
    }

    return static_cast<uint32_t>(value);
}

DNNFLiteral parseLiteral(const std::string& token, std::size_t line)
{
    std::size_t parsed = 0;
    long long value = 0;
    try {
        value = std::stoll(token, &parsed, 10);
    } catch (const std::exception&) {
        throw parseError(line, "expected signed literal, found '" + token + "'");
    }

    if (parsed != token.size() || value == 0 ||
        value < -static_cast<long long>(std::numeric_limits<int32_t>::max()) ||
        value > std::numeric_limits<int32_t>::max()) {
        throw parseError(line, "expected nonzero int32 literal, found '" + token + "'");
    }

    return static_cast<DNNFLiteral>(value);
}

DNNFKind parseKind(const std::string& token, std::size_t line)
{
    if (token == "o") {
        return Or;
    }
    if (token == "a") {
        return And;
    }
    if (token == "t") {
        return TrueVar;
    }
    if (token == "f") {
        return FalseVar;
    }

    throw parseError(line, "unknown node type '" + token + "'");
}

bool isNodeDeclaration(const std::string& token)
{
    return token == "o" || token == "a" || token == "t" || token == "f";
}

uint32_t variableOf(DNNFLiteral literal)
{
    return static_cast<uint32_t>(literal < 0 ? -static_cast<int64_t>(literal) : literal);
}

} // namespace

// Description:
//   Opens a .nnf file from disk and parses it into the DNNF structure declared
//   in DNNF.hpp.
// Variables:
//   path - Path to the .nnf file to parse.
// Return:
//   A DNNF populated with nodes, labelled arcs, its inferred root, and the
//   number of variables referenced by arc literals.
DNNF parseDNNF(const std::string& path)
{
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("failed to open DNNF file '" + path + "': " + std::strerror(errno));
    }

    return parseDNNF(input);
}

// Description:
//   Parses d4 NNF text from an input stream. Comment and blank lines are
//   ignored. Node lines have the form "type id 0", where type is o, a, t, or
//   f. Arc lines have the form "source target [literal ...] 0".
// Variables:
//   input - Stream containing the NNF data to parse.
// Return:
//   A DNNF whose root is the unique node with no incoming arc.
DNNF parseDNNF(std::istream& input)
{
    const auto lines = readMeaningfulLines(input);
    if (lines.empty()) {
        throw std::runtime_error("DNNF input must contain at least one node declaration");
    }

    DNNF dnnf;
    std::unordered_map<DNNFNodeId, std::size_t> nodeIndexes;
    std::vector<ParsedArc> arcs;

    for (const ParsedLine& line : lines) {
        if (line.tokens.empty()) {
            continue;
        }

        if (isNodeDeclaration(line.tokens[0])) {
            if (line.tokens.size() != 3 || line.tokens.back() != "0") {
                throw parseError(line.number, "node declaration must have the form 'type id 0'");
            }

            const DNNFNodeId id = parseUint32(line.tokens[1], line.number);
            if (id == 0) {
                throw parseError(line.number, "node id cannot be 0");
            }
            if (nodeIndexes.find(id) != nodeIndexes.end()) {
                throw parseError(line.number, "duplicate node id " + std::to_string(id));
            }

            nodeIndexes.emplace(id, dnnf.nodes.size());
            dnnf.nodes.push_back({parseKind(line.tokens[0], line.number), 0, id, {}});
            continue;
        }

        if (line.tokens.size() < 3 || line.tokens.back() != "0") {
            throw parseError(line.number, "arc must have the form 'source target [literal ...] 0'");
        }

        ParsedArc parsedArc;
        parsedArc.line = line.number;
        parsedArc.source = parseUint32(line.tokens[0], line.number);
        parsedArc.arc.target = parseUint32(line.tokens[1], line.number);
        if (parsedArc.source == 0 || parsedArc.arc.target == 0) {
            throw parseError(line.number, "arc node ids cannot be 0");
        }

        for (std::size_t i = 2; i + 1 < line.tokens.size(); ++i) {
            const DNNFLiteral literal = parseLiteral(line.tokens[i], line.number);
            parsedArc.arc.literals.push_back(literal);
            dnnf.numVars = std::max(dnnf.numVars, variableOf(literal));
        }
        arcs.push_back(std::move(parsedArc));
    }

    if (dnnf.nodes.empty()) {
        throw std::runtime_error("DNNF input must contain at least one node declaration");
    }

    std::unordered_set<DNNFNodeId> hasParent;
    for (ParsedArc& parsedArc : arcs) {
        const auto source = nodeIndexes.find(parsedArc.source);
        if (source == nodeIndexes.end()) {
            throw parseError(parsedArc.line, "arc uses undeclared source node " + std::to_string(parsedArc.source));
        }
        if (nodeIndexes.find(parsedArc.arc.target) == nodeIndexes.end()) {
            throw parseError(parsedArc.line, "arc uses undeclared target node " + std::to_string(parsedArc.arc.target));
        }

        DNNFNode& sourceNode = dnnf.nodes[source->second];
        if (sourceNode.kind == TrueVar || sourceNode.kind == FalseVar) {
            throw parseError(parsedArc.line, "constant leaf node " + std::to_string(parsedArc.source) + " cannot have outgoing arcs");
        }

        hasParent.insert(parsedArc.arc.target);
        sourceNode.children.push_back(std::move(parsedArc.arc));
        sourceNode.degree = static_cast<uint32_t>(sourceNode.children.size());
    }

    for (const DNNFNode& node : dnnf.nodes) {
        if (hasParent.find(node.id) == hasParent.end()) {
            if (dnnf.root != EMPTY_NODE) {
                throw std::runtime_error(
                    "DNNF must have exactly one root; found root candidates " +
                    std::to_string(dnnf.root) + " and " + std::to_string(node.id));
            }
            dnnf.root = node.id;
        }
    }

    if (dnnf.root == EMPTY_NODE) {
        throw std::runtime_error("DNNF has no root");
    }

    return dnnf;
}
