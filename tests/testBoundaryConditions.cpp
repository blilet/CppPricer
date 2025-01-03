#include "MeshUtils.hpp"
#include <iostream>
#include <cassert>
#include <cmath>

void testBoundaryConditions() {
    std::vector<std::vector<bool>> contour = {{true, false}, {false, true}};
    std::function<double(double, double)>  func = [](double x, double y) { return x + y; };
    BoundaryConditions bc(contour, func);
    assert(bc.apply(1.0, 2.0) == 3.0 && "apply function failed");
    assert(bc.check(0, 0) == true && "Initial check failed");
    bc.uncheck(0, 0);
    assert(bc.check(0, 0) == false && "Uncheck failed");
    bc.ToggleDir(true, true);
    bc.logMesh();
}
