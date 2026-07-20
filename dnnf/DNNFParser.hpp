#ifndef DNNF_PARSER_HPP
#define DNNF_PARSER_HPP

#include "DNNF.hpp"

#include <istream>
#include <string>

DNNF parseDNNF(const std::string& path);
DNNF parseDNNF(std::istream& input);

#endif
