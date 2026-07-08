#pragma once

#include "../plustimes/plustimes.hpp"
#include "../utils/alg_checks.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <cstddef>
#include <cstdint>
#include <optional>

using boost::multiprecision::cpp_int;

double countCore(
    const PlusTimesProgram& P,
    std::size_t ns,
    std::size_t nt,
    cpp_int theta,
    double kappa);

double countCoreWithSupportThreshold(
    const PlusTimesProgram& P,
    std::size_t ns,
    std::size_t nt,
    cpp_int theta,
    double kappa,
    std::size_t supportThreshold);

double counter(const PlusTimesProgram& P, double epsilon, double delta);

std::optional<std::size_t> supportSizeBounded(
    const PlusTimesProgram& P,
    NodeId o,
    std::size_t threshold);

void seedSamplingForTesting(std::uint32_t seed);
