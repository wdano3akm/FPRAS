#include "../dnnf/DNNFParser.hpp"

#include <cassert>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

void expectParseError(const std::string& text)
{
    std::istringstream input(text);
    bool threw = false;
    try {
        (void)parseDNNF(input);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    assert(threw);
}

} // namespace

int main()
{
    std::istringstream input(
        "c comments and blank lines are ignored\n"
        "\n"
        "o 1 0\n"
        "o 2 0\n"
        "o 3 0\n"
        "t 4 0\n"
        "3 4 -2 3 0\n"
        "3 4 2 0\n"
        "2 3 -1 0\n"
        "2 4 1 0\n"
        "1 2 0\n");

    const DNNF dnnf = parseDNNF(input);

    assert(dnnf.root == 1);
    assert(dnnf.numVars == 3);
    assert(dnnf.nodes.size() == 4);

    assert(dnnf.nodes[0].id == 1);
    assert(dnnf.nodes[0].kind == Or);
    assert(dnnf.nodes[0].degree == 1);
    assert(dnnf.nodes[0].children[0].target == 2);
    assert(dnnf.nodes[0].children[0].literals.empty());

    assert(dnnf.nodes[2].id == 3);
    assert(dnnf.nodes[2].degree == 2);
    assert(dnnf.nodes[2].children[0].target == 4);
    assert(dnnf.nodes[2].children[0].literals.size() == 2);
    assert(dnnf.nodes[2].children[0].literals[0] == -2);
    assert(dnnf.nodes[2].children[0].literals[1] == 3);

    assert(dnnf.nodes[3].id == 4);
    assert(dnnf.nodes[3].kind == TrueVar);
    assert(dnnf.nodes[3].degree == 0);

    const DNNF fileDnnf = parseDNNF("dnnf/input.nnf");
    assert(fileDnnf.root == 1);
    assert(fileDnnf.nodes.size() == 4);

    expectParseError("o 1 0\no 1 0\n");
    expectParseError("o 1 0\n1 2 0\n");
    expectParseError("t 1 0\no 2 0\n1 2 0\n");
    expectParseError("o 1 0\no 2 0\n");

    return 0;
}
