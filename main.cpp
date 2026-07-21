#include "algorithm/algorithm.hpp"
#include "cfg/CFGParser.hpp"
#include "dnnf/DNNFParser.hpp"
#include "plustimes/plustimes.hpp"
#include "plustimes/ptfromdnnf.hpp"
#include "utils/dnnf_ops.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

constexpr double epsilon = 0.5;
constexpr double delta = 0.25;
constexpr int defaultCfgLength = 3;

double runCFG(const std::string& path, int length)
{
    const CFG cfg = parseCFG(path);
    const PlusTimesProgram program = compileCFGToPlusTimes(cfg, length);
    return counter(program, epsilon, delta);
}

double runDNNF(const std::string& path)
{
    const DNNF dnnf = parseDNNF(path);
    const DNNF smoothDnnf = smoothDNNF(dnnf);
    const PlusTimesProgram program = compileDNNFtoPlusTimes(smoothDnnf);
    return counter(program, epsilon, delta);
}

std::string lowerExtension(const std::string& path)
{
    std::string extension = std::filesystem::path(path).extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return extension;
}

int parseCfgLength(const char* argument)
{
    const std::string text(argument);
    std::size_t parsed = 0;
    unsigned long value = 0;
    try {
        value = std::stoul(text, &parsed, 10);
    } catch (const std::exception&) {
        throw std::invalid_argument("CFG length must be a positive integer");
    }

    if (parsed != text.size() || value == 0 ||
        value > static_cast<unsigned long>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("CFG length must be a positive integer");
    }
    return static_cast<int>(value);
}

void printUsage(const char* executable)
{
    std::cerr << "usage: " << executable << " <input.nnf|input.cfg> [cfg-length]\n"
              << "       cfg-length defaults to " << defaultCfgLength << '\n';
}

} // namespace

int main(int argc, char* argv[])
{
    if (argc < 2 || argc > 3) {
        printUsage(argv[0]);
        return 2;
    }

    try {
        const std::string path(argv[1]);
        const std::string extension = lowerExtension(path);
        double value = 0.0;

        if (extension == ".nnf") {
            if (argc != 2) {
                throw std::invalid_argument("a DNNF input does not take a CFG length");
            }
            value = runDNNF(path);
        } else if (extension == ".cfg") {
            const int length = argc == 3 ? parseCfgLength(argv[2]) : defaultCfgLength;
            value = runCFG(path, length);
        } else {
            throw std::invalid_argument("input file must have a .nnf or .cfg extension");
        }

        std::cout << "value: " << value << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
