#pragma once

#include "MeshUtils.hpp"
#include <set>
#include <stdexcept>
#include <functional>
#include <utility>

class ItoDynamics {

private:
    const std::function<double(double, double, double)> partial_t;
    const std::function<double(double, double, double)> partial_x;

public:
    static std::pair<std::function<double(double, double, double)>, std::function<double(double, double, double)>>
    convert_S_to_x(const std::function<double(double, double, double)>& f,
                   const std::function<double(double, double, double)>& g);

    ItoDynamics(const std::function<double(double, double, double)>& partial_t,
                const std::function<double(double, double, double)>& partial_x);

    double getDrift(double t, double x, double p) const;
    double getPseudoVol(double t, double x, double p) const;
    

};

class ItoProcess {
private:
    FunctionMesh processMesh;
public:
    ItoProcess(const SpaceTimeMesh& stm);
    void solve(const BoundaryConditions& bc, const ItoDynamics& dynamics);
    double getVal(std::size_t i, std::size_t j);
    void logMesh() const;
    const FunctionMesh& getProcessMesh() const;

};

