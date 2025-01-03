#include "MeshUtils.hpp"
#include <iostream>
#include <cassert>

void testSpaceTimeMesh() {
    SpaceTimeMesh stm(0.0, 1.0, 2.0, 11, 20);
    assert(stm.get_N() == 10 && "get_N failed");
    assert(stm.get_N_T() == 20 && "get_N_T failed");
    assert(stm.get_R() == 1.0 && "get_R failed");
    assert(stm.get_T() == 2.0 && "get_T failed");
    assert(std::abs(stm.get_dx() - 0.1) < 1e-6 && "get_dx failed");
    assert(std::abs(stm.get_dt() - 0.1) < 1e-6 && "get_dt failed");
    auto coords = stm.getCoords(5, 10);
    assert(std::abs(coords.first - 0.5) < 1e-6 && "getCoords x failed"); // precision 10^-6?
    assert(std::abs(coords.second - 1.0) < 1e-6 && "getCoords t failed");
}
