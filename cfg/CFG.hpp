#ifndef CFG_HPP
#define CFG_HPP

#include <cstdint>
#include <vector>
#include <stdint.h>
using namespace std;

using SymbolId = uint32_t;
using NonterminalId = uint32_t;
using TerminalId = uint32_t;

struct TerminalRule {
    NonterminalId lhs;
    TerminalId terminal;
};

struct BinaryRule {
    NonterminalId lhs;
    NonterminalId left;
    NonterminalId right;
};

struct CFG {
    NonterminalId start;
    vector<NonterminalId> nonterminalIds;
    vector<TerminalId> terminalIds;
    vector<TerminalRule> terminalRules;
    vector<BinaryRule> binaryRules;
    bool startNullable = false;
};

#endif
