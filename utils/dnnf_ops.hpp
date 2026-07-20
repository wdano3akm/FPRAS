#ifndef DNNF_OPS_HPP
#define DNNF_OPS_HPP

#include "../dnnf/DNNF.hpp"

// Tests both conditions of Darwiche's Definition 6.
bool isSmoothDNNF(const DNNF& dnnf);

// Returns an equivalent smooth DNNF without modifying the input.
DNNF smoothDNNF(const DNNF& dnnf);

// Backwards-compatible mutating spelling retained from the work in progress.
DNNF smooth_DNNF(DNNF& dnnf);

#endif
