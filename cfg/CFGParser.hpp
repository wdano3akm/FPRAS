#ifndef CFG_PARSER_HPP
#define CFG_PARSER_HPP

#include "CFG.hpp"

#include <istream>
#include <string>

CFG parseCFG(const std::string& path);
CFG parseCFG(std::istream& input);

#endif
