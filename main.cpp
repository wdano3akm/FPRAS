#include "algorithm/algorithm.hpp"
#include "cfg/CFGParser.hpp"
#include "dnnf/DNNFParser.hpp"
#include "plustimes/plustimes.hpp"
#include "plustimes/ptfromdnnf.hpp"
#include "utils/dnnf_ops.hpp"

#include <algorithm>
#include<cmath>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <limits>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

constexpr double epsilon = 0.5;
constexpr double delta = 0.25;
constexpr int defaultCfgLength = 3;
constexpr std::size_t N = 4;

double runCFG(const std::string& path, int length, bool programSizeOnly)
{
    const CFG cfg = parseCFG(path);
    const PlusTimesProgram program = compileCFGToPlusTimes(cfg, length);
    std::cout << "|P|: " << program.nodes.size() << '\n';
    if (programSizeOnly) {
        return 0.0;
    }

    const std::size_t Psize = program.nodes.size();
    const std::size_t n = static_cast<std::size_t>(program.degree);
    const auto ns = 4;
    const auto nt = 8;
    const double kappa =  0.25f / (4.0 * static_cast<double>(n));

    const cpp_int theta = cpp_int(512) * ns * nt * n * Psize;
    return countCore(program, ns, nt, theta, kappa);
}

double runDNNF(const std::string& path, bool programSizeOnly)
{
    const DNNF dnnf = parseDNNF(path);
    const DNNF smoothDnnf = smoothDNNF(dnnf);
    const PlusTimesProgram program = compileDNNFtoPlusTimes(smoothDnnf);
    std::cout << "|P|: " << program.nodes.size() << '\n';
    if (programSizeOnly) {
        return 0.0;
    }
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
    std::cerr << "usage: " << executable << " [--program-size-only] <input.nnf|input.cfg> [cfg-length]\n"
              << "       --program-size-only prints |P| and exits before counting\n"
              << "       cfg-length defaults to " << defaultCfgLength << '\n';
}

} // namespace

int main(int argc, char* argv[])
{
    const bool programSizeOnly = argc >= 2 && std::string(argv[1]) == "--program-size-only";
    const int firstInputArgument = programSizeOnly ? 2 : 1;
    if (argc < firstInputArgument + 1 || argc > firstInputArgument + 2) {
        printUsage(argv[0]);
        return 2;
    }

    try {
        const std::string path(argv[firstInputArgument]);
        const std::string extension = lowerExtension(path);
        std::vector<double> values(N);
        std::exception_ptr workerError;
        std::mutex workerErrorMutex;

        const auto runTest = [&](std::size_t run) {
            try {
                if (extension == ".nnf") {
                    values[run] = runDNNF(path, programSizeOnly);
                } else {
                    const int length = argc == firstInputArgument + 2
                        ? parseCfgLength(argv[firstInputArgument + 1])
                        : defaultCfgLength;
                    values[run] = runCFG(path, length, programSizeOnly);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(workerErrorMutex);
                if (!workerError) {
                    workerError = std::current_exception();
                }
            }
        };

        if (extension == ".nnf") {
            if (argc != firstInputArgument + 1) {
                throw std::invalid_argument("a DNNF input does not take a CFG length");
            }
        } else if (extension == ".cfg") {
            if (argc == firstInputArgument + 2) {
                parseCfgLength(argv[firstInputArgument + 1]);
            }
        } else {
            throw std::invalid_argument("input file must have a .nnf or .cfg extension");
        }

        if (programSizeOnly) {
            runTest(0);
            if (workerError) {
                std::rethrow_exception(workerError);
            }
            return 0;
        }

        std::vector<std::thread> workers;
        workers.reserve(N);
        for (std::size_t run = 0; run < N; ++run) {
            workers.emplace_back(runTest, run);
        }
        for (std::thread& worker : workers) {
            worker.join();
        }
        if (workerError) {
            std::rethrow_exception(workerError);
        }

        for (std::size_t run = 0; run < N; ++run) {
            std::cout << "value [" << run << "]: " << values[run] << '\n';
            std::cout << "error [" << run << "]: "
                      << std::abs(values[run] - std::exp2(22)) / std::exp2(22)
                      << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
