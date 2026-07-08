#pragma once

#include "../plustimes/plustimes.hpp"

#include <unordered_set>
#include <vector>

using VarId = PolyVarId;
using VariableSets = std::vector<std::unordered_set<VarId>>;

void requireReady(bool condition, const char* message);

std::vector<int> computeDegreesOrThrow(const PlusTimesProgram& P);

bool isTopologicallyOrdered(const PlusTimesProgram& P);
bool noConstants(const PlusTimesProgram& P);
bool allProductsBinary(const PlusTimesProgram& P);
bool noSumHasSumChild(const PlusTimesProgram& P);
bool isHomogeneous(const PlusTimesProgram& P, const std::vector<int>& degrees);
bool disjointVarSets(
    const std::unordered_set<VarId>& left,
    const std::unordered_set<VarId>& right);
bool isStructurallyMultilinear(const PlusTimesProgram& P);

int rootHeightOrThrow(const PlusTimesProgram& P);

void assertReadyForFPRAS(const PlusTimesProgram& P);
