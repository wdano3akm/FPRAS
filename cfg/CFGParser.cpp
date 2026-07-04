#include "CFGParser.hpp"

#include <cerrno>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

struct ParsedLine {
    std::size_t number;
    std::vector<std::string> tokens;
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
    return std::runtime_error("CFG parse error on line " + std::to_string(line) + ": " + message);
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

void requireNoTrailingTokens(const ParsedLine& line, std::size_t expected)
{
    if (line.tokens.size() != expected) {
        throw parseError(line.number, "unexpected token count");
    }
}

bool contains(const std::unordered_set<uint32_t>& ids, uint32_t id)
{
    return ids.find(id) != ids.end();
}

std::unordered_set<uint32_t> makeIdSet(const std::vector<uint32_t>& ids, const char* kind)
{
    std::unordered_set<uint32_t> seen;
    for (uint32_t id : ids) {
        if (id == 0) {
            throw std::runtime_error(std::string(kind) + " id cannot be 0");
        }

        if (!seen.insert(id).second) {
            throw std::runtime_error("duplicate " + std::string(kind) + " id " + std::to_string(id));
        }
    }

    return seen;
}

std::vector<uint32_t> parseSymbolLine(const ParsedLine& line, const std::string& prefix, std::size_t expectedCount)
{
    if (line.tokens.empty() || line.tokens[0] != prefix) {
        throw parseError(line.number, "expected '" + prefix + "' declaration");
    }

    if (line.tokens.size() != expectedCount + 1) {
        throw parseError(line.number, "expected " + std::to_string(expectedCount) + " ids after '" + prefix + "'");
    }

    std::vector<uint32_t> ids;
    ids.reserve(expectedCount);
    for (std::size_t i = 1; i < line.tokens.size(); ++i) {
        ids.push_back(parseUint32(line.tokens[i], line.number));
    }

    return ids;
}

void addRuleAlternative(
    const std::vector<uint32_t>& rhs,
    uint32_t lhs,
    const std::unordered_set<uint32_t>& terminals,
    const std::unordered_set<uint32_t>& nonterminals,
    CFG& cfg,
    std::size_t line)
{
    if (rhs.empty()) {
        if (lhs == cfg.start) {
            cfg.startNullable = true;
            return;
        }

        throw parseError(line, "only the start nonterminal may have an empty production");
    }

    if (rhs.size() == 1) {
        const uint32_t terminal = rhs[0];
        if (!contains(terminals, terminal)) {
            throw parseError(line, "terminal production uses undeclared terminal " + std::to_string(terminal));
        }

        cfg.terminalRules.push_back({lhs, terminal});
        return;
    }

    if (rhs.size() == 2) {
        const uint32_t left = rhs[0];
        const uint32_t right = rhs[1];
        if (!contains(nonterminals, left) || !contains(nonterminals, right)) {
            throw parseError(line, "binary production must use declared nonterminals");
        }

        cfg.binaryRules.push_back({lhs, left, right});
        return;
    }

    throw parseError(line, "CNF productions must have either one terminal or two nonterminals");
}

void parseRuleLine(
    const ParsedLine& line,
    const std::unordered_set<uint32_t>& terminals,
    const std::unordered_set<uint32_t>& nonterminals,
    CFG& cfg)
{
    if (line.tokens.size() < 2) {
        throw parseError(line.number, "rule line must contain lhs and terminating 0");
    }

    const uint32_t lhs = parseUint32(line.tokens[0], line.number);
    if (!contains(nonterminals, lhs)) {
        throw parseError(line.number, "rule lhs uses undeclared nonterminal " + std::to_string(lhs));
    }

    if (line.tokens.back() != "0") {
        throw parseError(line.number, "rule line must terminate with 0");
    }

    std::vector<uint32_t> rhs;
    bool sawAlternative = false;
    for (std::size_t i = 1; i + 1 < line.tokens.size(); ++i) {
        if (line.tokens[i] == "|") {
            if (rhs.empty()) {
                throw parseError(line.number, "alternative separator cannot have an empty left side");
            }

            addRuleAlternative(rhs, lhs, terminals, nonterminals, cfg, line.number);
            rhs.clear();
            sawAlternative = true;
            continue;
        }

        rhs.push_back(parseUint32(line.tokens[i], line.number));
    }

    if (sawAlternative && rhs.empty()) {
        throw parseError(line.number, "alternative separator cannot appear immediately before terminating 0");
    }

    addRuleAlternative(rhs, lhs, terminals, nonterminals, cfg, line.number);
}

} // namespace

// Description:
//   Opens a .cnfg file from disk and parses it into the CFG structure declared
//   in CFG.hpp.
// Variables:
//   path - Path to the .cnfg file to parse.
// Return:
//   A CFG populated with terminal ids, nonterminal ids, terminal rules,
//   binary rules, the inferred start symbol, and startNullable.
CFG parseCFG(const std::string& path)
{
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("failed to open CFG file '" + path + "': " + std::strerror(errno));
    }

    return parseCFG(input);
}

// Description:
//   Parses CNFG text from an input stream. Comment and blank lines are ignored.
//   The first non-comment line must be "p cfg num_terminals num_nonterminals
//   num_rules", followed by terminal and nonterminal declarations, then CNF
//   rule lines terminated by 0. The first declared nonterminal is used as the
//   CFG start symbol.
// Variables:
//   input - Stream containing the CNFG data to parse.
// Return:
//   A CFG populated with the symbols and rules read from the stream.
CFG parseCFG(std::istream& input)
{
    const auto lines = readMeaningfulLines(input);
    if (lines.size() < 3) {
        throw std::runtime_error("CFG input must contain p, t, and n declarations");
    }

    const ParsedLine& problem = lines[0];
    if (problem.tokens.size() != 5 || problem.tokens[0] != "p" || problem.tokens[1] != "cfg") {
        throw parseError(problem.number, "expected 'p cfg num_terminals num_nonterminals num_rules'");
    }

    const std::size_t numTerminals = parseUint32(problem.tokens[2], problem.number);
    const std::size_t numNonterminals = parseUint32(problem.tokens[3], problem.number);
    const std::size_t numRules = parseUint32(problem.tokens[4], problem.number);
    requireNoTrailingTokens(problem, 5);

    CFG cfg;
    cfg.terminalIds = parseSymbolLine(lines[1], "t", numTerminals);
    cfg.nonterminalIds = parseSymbolLine(lines[2], "n", numNonterminals);
    if (cfg.nonterminalIds.empty()) {
        throw std::runtime_error("CFG must declare at least one nonterminal");
    }
    cfg.start = cfg.nonterminalIds.front();

    const auto terminals = makeIdSet(cfg.terminalIds, "terminal");
    const auto nonterminals = makeIdSet(cfg.nonterminalIds, "nonterminal");
    for (uint32_t terminal : cfg.terminalIds) {
        if (contains(nonterminals, terminal)) {
            throw std::runtime_error("terminal and nonterminal sets must be disjoint; duplicate id " + std::to_string(terminal));
        }
    }

    if (lines.size() != numRules + 3) {
        throw std::runtime_error(
            "CFG declared " + std::to_string(numRules) + " rule lines but found " +
            std::to_string(lines.size() - 3));
    }

    for (std::size_t i = 0; i < numRules; ++i) {
        parseRuleLine(lines[i + 3], terminals, nonterminals, cfg);
    }

    return cfg;
}
