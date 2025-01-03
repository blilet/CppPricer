#include "ItoProcess.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

void testItoDynamics() {
    std::function<double(double, double, double)>  partial_t = [](double t, double x, double p) { return t + x + p; };
    std::function<double(double, double, double)>  partial_x = [](double t, double x, double p) { return t * x * p; };
    ItoDynamics dynamics(partial_t, partial_x);
    
    double drift = dynamics.getDrift(1.0, 2.0, 3.0);
    double pseudoVol = dynamics.getPseudoVol(1.0, 2.0, 3.0);
    assert(std::abs(drift - 6.0) < 1e-6 && "getDrift failed");
    assert(std::abs(pseudoVol - 6.0) < 1e-6 && "getPseudoVol failed");
    
    auto [converted_t, converted_x] = ItoDynamics::convert_S_to_x(partial_t, partial_x);
    assert(std::abs(converted_t(1.0, 2.0, 3.0) - 6.0) < 1e-6 && "convert_S_to_x partial_t failed");
    assert(std::abs(converted_x(1.0, 2.0, 3.0) - 6.0) < 1e-6 && "convert_S_to_x partial_x failed");
    
}
