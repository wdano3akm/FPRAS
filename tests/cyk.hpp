#ifndef TESTS_CYK_HPP
#define TESTS_CYK_HPP

#include "../cfg/CFG.hpp"

#include <vector>

bool cykAccepts(const CFG& cfg, const std::vector<TerminalId>& word);

#endif
