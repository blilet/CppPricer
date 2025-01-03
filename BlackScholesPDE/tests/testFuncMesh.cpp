#include "MeshUtils.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

void testFunctionMesh() {
    SpaceTimeMesh stm(0.0, 1.0, 2.0, 11, 20);
    FunctionMesh fm(stm);
    assert(fm.getNumRows() == stm.get_N() && "getNumRows failed");
    assert(fm.getNumCols() == stm.get_N_T() && "getNumCols failed");
    fm.setMeshData(5, 10, 3.14);
    assert(std::abs(fm.getMeshData(5, 10) - 3.14) < 1e-6 && "getMeshData failed");
    std::vector<std::vector<bool>> contour(10, std::vector<bool>(20, false));
    contour[5][10] = true;
    std::function<double( double, double)>  func = [](double x, double y) { return x * y; };
    BoundaryConditions bc(contour, func);
    fm.applyBoundaryConditions(bc);

    assert(std::abs(fm.getMeshData(5, 10) - 5.0) < 1e-6 && "Boundary condition application failed");
    fm.logMesh();
}
