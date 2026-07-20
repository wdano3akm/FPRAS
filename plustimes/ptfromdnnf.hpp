#pragma once

#include "../dnnf/DNNF.hpp"
#include "plustimes.hpp"

PlusTimesProgram compileDNNFtoPlusTimes(const DNNF& dnnf);
