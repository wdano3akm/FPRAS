#include "../cfg/CFG.hpp"
#include "../plustimes/plustimes.hpp"
#include "cyk.hpp"

#include <algorithm>
#include <cassert>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using Monomial = std::vector<PolyVarId>;
using Support = std::set<Monomial>;
using Word = std::vector<TerminalId>;
using LanguageSlice = std::set<Word>;

Support supportOf(const PlusTimesProgram& program, NodeId node, std::vector<Support>& memo)
{
    if (node == EMPTY_NODE) {
        return {};
    }

    assert(node < program.nodes.size());
    if (!memo[node].empty()) {
        return memo[node];
    }

    const PTNode& ptNode = program.nodes[node];
    Support support;

    if (ptNode.kind == PTKind::Var) {
        support.insert({ptNode.var});
    } else if (ptNode.kind == PTKind::Add) {
        for (NodeId child : ptNode.children) {
            Support childSupport = supportOf(program, child, memo);
            support.insert(childSupport.begin(), childSupport.end());
        }
    } else {
        assert(ptNode.kind == PTKind::Mul);
        assert(ptNode.children.size() == 2);

        Support leftSupport = supportOf(program, ptNode.children[0], memo);
        Support rightSupport = supportOf(program, ptNode.children[1], memo);

        for (const Monomial& left : leftSupport) {
            for (const Monomial& right : rightSupport) {
                Monomial product = left;
                product.insert(product.end(), right.begin(), right.end());
                std::sort(product.begin(), product.end());
                support.insert(std::move(product));
            }
        }
    }

    memo[node] = support;
    return memo[node];
}

Support exactSupport(const PlusTimesProgram& program)
{
    std::vector<Support> memo(program.nodes.size());
    return supportOf(program, program.root, memo);
}

std::size_t exactSupportSize(const PlusTimesProgram& program)
{
    return exactSupport(program).size();
}

PolyVarId varId(std::size_t terminalIndex, int offset, std::size_t numTerminals)
{
    return static_cast<PolyVarId>(static_cast<std::size_t>(offset) * numTerminals + terminalIndex);
}

Word wordFromMonomial(const Monomial& monomial, const CFG& grammar, int n)
{
    Word word(static_cast<std::size_t>(n), 0);
    std::vector<bool> seen(static_cast<std::size_t>(n), false);
    const std::size_t numTerminals = grammar.terminalIds.size();

    for (PolyVarId var : monomial) {
        const std::size_t offset = var / numTerminals;
        const std::size_t terminalIndex = var % numTerminals;

        if (offset >= static_cast<std::size_t>(n) || terminalIndex >= numTerminals || seen[offset]) {
            throw std::logic_error("plus-times support contains a malformed positional monomial");
        }

        word[offset] = grammar.terminalIds[terminalIndex];
        seen[offset] = true;
    }

    if (std::find(seen.begin(), seen.end(), false) != seen.end()) {
        throw std::logic_error("plus-times support contains a monomial with missing positions");
    }

    return word;
}

LanguageSlice languageFromPlusTimes(const PlusTimesProgram& program, const CFG& grammar, int n)
{
    LanguageSlice language;
    for (const Monomial& monomial : exactSupport(program)) {
        assert(monomial.size() == static_cast<std::size_t>(n));
        language.insert(wordFromMonomial(monomial, grammar, n));
    }

    return language;
}

void enumerateCykLanguage(
    const CFG& grammar,
    int n,
    Word& prefix,
    LanguageSlice& language)
{
    if (prefix.size() == static_cast<std::size_t>(n)) {
        if (cykAccepts(grammar, prefix)) {
            language.insert(prefix);
        }
        return;
    }

    for (TerminalId terminal : grammar.terminalIds) {
        prefix.push_back(terminal);
        enumerateCykLanguage(grammar, n, prefix, language);
        prefix.pop_back();
    }
}

LanguageSlice languageFromCykBruteforce(const CFG& grammar, int n)
{
    LanguageSlice language;
    Word prefix;
    enumerateCykLanguage(grammar, n, prefix, language);
    return language;
}

void assertMatchesCykBruteforce(const CFG& grammar, int n)
{
    const PlusTimesProgram program = compileCFGToPlusTimes(grammar, n);
    assert((languageFromPlusTimes(program, grammar, n) == languageFromCykBruteforce(grammar, n)));
}

CFG cfg(
    NonterminalId start,
    std::vector<NonterminalId> nonterminals,
    std::vector<TerminalId> terminals,
    std::vector<TerminalRule> terminalRules,
    std::vector<BinaryRule> binaryRules)
{
    return CFG{start, nonterminals, terminals, terminalRules, binaryRules, false};
}

void testOneTerminal()
{
    const NonterminalId S = 10;
    const TerminalId a = 1;
    const CFG grammar = cfg(S, {S}, {a}, {{S, a}}, {});

    const PlusTimesProgram program = compileCFGToPlusTimes(grammar, 1);

    assert(program.root != EMPTY_NODE);
    assert(program.nodes[program.root].degree == 1);
    assert((exactSupport(program) == Support{{varId(0, 0, 1)}}));
    assert(exactSupportSize(program) == 1);
    assertMatchesCykBruteforce(grammar, 1);
}

void testWrongLengthGivesEmpty()
{
    const NonterminalId S = 10;
    const TerminalId a = 1;
    const CFG grammar = cfg(S, {S}, {a}, {{S, a}}, {});

    const PlusTimesProgram program = compileCFGToPlusTimes(grammar, 2);

    assert(program.root == EMPTY_NODE);
    assert(exactSupportSize(program) == 0);
    assertMatchesCykBruteforce(grammar, 2);
}

void testSimpleConcatenation()
{
    const NonterminalId S = 10;
    const NonterminalId A = 11;
    const NonterminalId B = 12;
    const TerminalId a = 1;
    const TerminalId b = 2;
    const CFG grammar = cfg(S, {S, A, B}, {a, b}, {{A, a}, {B, b}}, {{S, A, B}});

    const PlusTimesProgram program = compileCFGToPlusTimes(grammar, 2);

    assert((exactSupport(program) == Support{{varId(0, 0, 2), varId(1, 1, 2)}}));
    assert(exactSupportSize(program) == 1);
    assertMatchesCykBruteforce(grammar, 2);
}

void testTwoAlternatives()
{
    const NonterminalId S = 10;
    const NonterminalId A = 11;
    const NonterminalId B = 12;
    const NonterminalId C = 13;
    const TerminalId a = 1;
    const TerminalId b = 2;
    const TerminalId c = 3;
    const CFG grammar = cfg(
        S,
        {S, A, B, C},
        {a, b, c},
        {{A, a}, {B, b}, {C, c}},
        {{S, A, B}, {S, A, C}});

    const PlusTimesProgram program = compileCFGToPlusTimes(grammar, 2);
    const Support expected = {
        {varId(0, 0, 3), varId(1, 1, 3)},
        {varId(0, 0, 3), varId(2, 1, 3)},
    };

    assert((exactSupport(program) == expected));
    assert(exactSupportSize(program) == 2);
    assertMatchesCykBruteforce(grammar, 2);
}

void testAmbiguityDoesNotDuplicateWords()
{
    const NonterminalId S = 10;
    const NonterminalId A = 11;
    const NonterminalId B = 12;
    const NonterminalId C = 13;
    const NonterminalId D = 14;
    const TerminalId a = 1;
    const TerminalId b = 2;
    const CFG grammar = cfg(
        S,
        {S, A, B, C, D},
        {a, b},
        {{A, a}, {B, b}, {C, a}, {D, b}},
        {{S, A, B}, {S, C, D}});

    const PlusTimesProgram program = compileCFGToPlusTimes(grammar, 2);

    assert((exactSupport(program) == Support{{varId(0, 0, 2), varId(1, 1, 2)}}));
    assert(exactSupportSize(program) == 1);
    assertMatchesCykBruteforce(grammar, 2);
}

void testOrderMatters()
{
    const NonterminalId S = 10;
    const NonterminalId A = 11;
    const NonterminalId B = 12;
    const TerminalId a = 1;
    const TerminalId b = 2;
    const CFG grammar = cfg(
        S,
        {S, A, B},
        {a, b},
        {{A, a}, {B, b}},
        {{S, A, B}, {S, B, A}});

    const PlusTimesProgram program = compileCFGToPlusTimes(grammar, 2);
    const Support expected = {
        {varId(0, 0, 2), varId(1, 1, 2)},
        {varId(1, 0, 2), varId(0, 1, 2)},
    };

    assert((exactSupport(program) == expected));
    assert(exactSupportSize(program) == 2);
    assertMatchesCykBruteforce(grammar, 2);
}

void testLengthSplits()
{
    const NonterminalId S = 10;
    const NonterminalId A = 11;
    const NonterminalId B = 12;
    const TerminalId a = 1;
    const TerminalId b = 2;
    const CFG grammar = cfg(
        S,
        {S, A, B},
        {a, b},
        {{A, a}, {B, b}},
        {{S, A, B}, {A, A, A}});

    const PlusTimesProgram program = compileCFGToPlusTimes(grammar, 3);

    assert((exactSupport(program) == Support{{varId(0, 0, 2), varId(0, 1, 2), varId(1, 2, 2)}}));
    assert(exactSupportSize(program) == 1);
    assertMatchesCykBruteforce(grammar, 3);
}

} // namespace

int main()
{
    testOneTerminal();
    testWrongLengthGivesEmpty();
    testSimpleConcatenation();
    testTwoAlternatives();
    testAmbiguityDoesNotDuplicateWords();
    testOrderMatters();
    testLengthSplits();

    return 0;
}
