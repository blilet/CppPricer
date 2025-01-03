#include "ItoProcess.hpp"
#include <cmath>

std::pair<std::function<double(double, double, double)>, std::function<double(double, double, double)>>
ItoDynamics::convert_S_to_x(const std::function<double(double, double, double)>& f,
                            const std::function<double(double, double, double)>& g) {

    auto func1 = [&f, &g](double t, double x, double v) -> double {
        double exp_x = std::exp(x);
        return f(t, exp_x, v) + g(t, exp_x, v) * exp_x / 2;
    };
    
    auto func2 = [&f, &g](double t, double x, double v) -> double {
        double exp_x = std::exp(x);
        return g(t, exp_x, v) * exp_x;
    };
    
    return {func1, func2};
}

ItoDynamics::ItoDynamics(const std::function<double(double, double, double)>& partial_t,
                         const std::function<double(double, double, double)>& partial_x)
    : partial_t(partial_t), partial_x(partial_x) {}

double ItoDynamics::getDrift(double t, double x, double p) const {
    return partial_t(t, x, p);
}

double ItoDynamics::getPseudoVol(double t, double x, double p) const {
    return partial_x(t, x, p);
}

ItoProcess::ItoProcess(const SpaceTimeMesh& stm) :processMesh(stm) {

}

void ItoProcess::solve(const BoundaryConditions& bc, const ItoDynamics& dynamics) {
    processMesh.applyBoundaryConditions(bc);
    std::set<std::pair<int, int>> visited;
    std::set<std::pair<int, int>> not_visited;

    for (int x = 0; x < processMesh.getNumRows(); x++) {
        for (int y = 0; y < processMesh.getNumCols(); y++) {
            if (bc.check(x, y)) {
                visited.insert({x, y});
            } else {
                not_visited.insert({x, y});
            }
        }
    }

    std::size_t overflow_counter = 0;
    std::set<std::pair<int, int>> dirs = {{0, 1}, {0, -1}, {-1, 0}, {1, 0}};
    double dx = processMesh.getSpaceTimeMesh().get_dx();
    double dt = processMesh.getSpaceTimeMesh().get_dt();
    std::size_t max_iter = (processMesh.getNumCols()+1)*(processMesh.getNumRows() +1);
    while (!not_visited.empty()) {
        if (overflow_counter > max_iter) {
            throw std::runtime_error("Convergence error: Boundaries aren't sufficient.");
            return;
        }
        for (auto& p : visited) {
            for (auto& dir : dirs) {
                if (not_visited.find({p.first + dir.first, p.second + dir.second}) != not_visited.end()) {
                    not_visited.erase({p.first + dir.first, p.second + dir.second});
                    visited.insert({p.first + dir.first, p.second + dir.second});
                    std::pair<double, double> spaceTimeCoords = processMesh.getSpaceTimeMesh().getCoords(p.first, p.second);
                    double new_val = processMesh.getMeshData(p.first,p.second) +
                    dir.first * dynamics.getDrift(spaceTimeCoords.first, spaceTimeCoords.second, processMesh.getMeshData(p.first,p.second))*dx +
                                                  dir.second * dynamics.getPseudoVol(spaceTimeCoords.first, spaceTimeCoords.second, processMesh.getMeshData(p.first,p.second)) *dt;
                    processMesh.setMeshData(p.first + dir.first,p.second + dir.second, new_val);
                }
            }
        }
        overflow_counter++;
    }
}
double ItoProcess::getVal(std::size_t i, std::size_t j){
    return processMesh.getMeshData(i, j);
}
void ItoProcess::logMesh() const{
    processMesh.logMesh();
}
const FunctionMesh& ItoProcess::getProcessMesh() const{
    return processMesh;
}

