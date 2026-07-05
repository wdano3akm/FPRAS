#include "cyk.hpp"

#include <set>
#include <vector>

bool cykAccepts(const CFG& cfg, const std::vector<TerminalId>& word)
{
    const int n = static_cast<int>(word.size());
    if (n == 0) {
        return cfg.startNullable;
    }

    std::vector<std::vector<std::set<NonterminalId>>> table(
        static_cast<std::size_t>(n),
        std::vector<std::set<NonterminalId>>(static_cast<std::size_t>(n)));

    for (int offset = 0; offset < n; ++offset) {
        for (const TerminalRule& rule : cfg.terminalRules) {
            if (rule.terminal == word[static_cast<std::size_t>(offset)]) {
                table[static_cast<std::size_t>(offset)][static_cast<std::size_t>(offset)].insert(rule.lhs);
            }
        }
    }

    for (int length = 2; length <= n; ++length) {
        for (int offset = 0; offset <= n - length; ++offset) {
            const int end = offset + length - 1;

            for (int split = offset; split < end; ++split) {
                const std::set<NonterminalId>& leftSet =
                    table[static_cast<std::size_t>(offset)][static_cast<std::size_t>(split)];
                const std::set<NonterminalId>& rightSet =
                    table[static_cast<std::size_t>(split + 1)][static_cast<std::size_t>(end)];

                for (const BinaryRule& rule : cfg.binaryRules) {
                    if (leftSet.find(rule.left) != leftSet.end() &&
                        rightSet.find(rule.right) != rightSet.end()) {
                        table[static_cast<std::size_t>(offset)][static_cast<std::size_t>(end)].insert(rule.lhs);
                    }
                }
            }
        }
    }

    return table[0][static_cast<std::size_t>(n - 1)].find(cfg.start) !=
           table[0][static_cast<std::size_t>(n - 1)].end();
}
