#include "../cfg/CFGParser.hpp"

#include <cassert>
#include <sstream>

int main()
{
    std::istringstream input(
        "c comments and blank lines are ignored\n"
        "\n"
        "p cfg 2 3 4\n"
        "t 1 2\n"
        "n 10 11 12\n"
        "10 11 12 0\n"
        "10 1 0\n"
        "11 1 0\n"
        "12 2 0\n");

    const CFG cfg = parseCFG(input);

    assert(cfg.start == 10);
    assert(!cfg.startNullable);

    assert(cfg.terminalIds.size() == 2);
    assert(cfg.terminalIds[0] == 1);
    assert(cfg.terminalIds[1] == 2);

    assert(cfg.nonterminalIds.size() == 3);
    assert(cfg.nonterminalIds[0] == 10);
    assert(cfg.nonterminalIds[1] == 11);
    assert(cfg.nonterminalIds[2] == 12);

    assert(cfg.binaryRules.size() == 1);
    assert(cfg.binaryRules[0].lhs == 10);
    assert(cfg.binaryRules[0].left == 11);
    assert(cfg.binaryRules[0].right == 12);

    assert(cfg.terminalRules.size() == 3);
    assert(cfg.terminalRules[0].lhs == 10);
    assert(cfg.terminalRules[0].terminal == 1);
    assert(cfg.terminalRules[1].lhs == 11);
    assert(cfg.terminalRules[1].terminal == 1);
    assert(cfg.terminalRules[2].lhs == 12);
    assert(cfg.terminalRules[2].terminal == 2);

    return 0;
}
