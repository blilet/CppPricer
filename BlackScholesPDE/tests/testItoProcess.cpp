#include "ItoProcess.hpp"
#include "MeshUtils.hpp"
#include <cassert>
#include <iostream>

void testItoProcess() {
    SpaceTimeMesh stm(0.0, 1.0, 2.0, 11, 20);
    ItoProcess process(stm);
    std::vector<std::vector<bool>> contour(10, std::vector<bool>(20, false));
    contour[0][0] = true;
    std::function<double(double, double)>  boundaryFunc = [](double x, double t) { return x + t; };
    BoundaryConditions bc(contour, boundaryFunc);

    std::function<double(double, double, double)>  partial_t = [](double t, double x, double p) { return t + x; };
    std::function<double(double, double, double)>  partial_x = [](double t, double x, double p) { return p * x; };
    ItoDynamics dynamics(partial_t, partial_x);

    process.solve(bc, dynamics);
    double val = process.getVal(0, 0);
    assert(std::abs(val - boundaryFunc(0.0, 0.0)) < 1e-6 && "getVal failed");
    const FunctionMesh& mesh = process.getProcessMesh();
    assert(mesh.getNumRows() == stm.get_N() && "Process mesh row size mismatch");
    assert(mesh.getNumCols() == stm.get_N_T() && "Process mesh column size mismatch");
    process.logMesh();

}
