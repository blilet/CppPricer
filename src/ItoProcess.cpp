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
    return partial_t(t, x, p); // attention t,x,p non pas x,t,p (alors que functionmesh etc est x,t)
}

double ItoDynamics::getPseudoVol(double t, double x, double p) const {
    return partial_x(t, x, p);
}

ItoProcess::ItoProcess(const SpaceTimeMesh& stm) :processMesh(stm) {

}

void ItoProcess::solve(const BoundaryConditions& bc, const ItoDynamics& dynamics) {
    processMesh.applyBoundaryConditions(bc);
    double dx = processMesh.getSpaceTimeMesh().get_dx();
    double dt = processMesh.getSpaceTimeMesh().get_dt();
    bool flag1 = true; // check for easy case, itoprocess given at t= 0
    for (int x = 0; x < processMesh.getNumRows(); x++) {
            if (!bc.check(x, 0)) {
                flag1 = false;
            }
    }
    bool flag2 = true; // check for easy case, itoprocess given at t= T
    for (int x = 0; x < processMesh.getNumRows(); x++) {
            if (!bc.check(x, processMesh.getNumCols()-1)) {
                flag2 = false;
            }
    }
    bool flag3 = true; // check for easy case, itoprocess given at x= 0 (i.e inf x)
    for (int y= 0; y < processMesh.getNumCols(); y++) {
            if (!bc.check(0, y)) {
                flag3 = false;
            }
    }
    bool flag4 = true; // check for easy case, itoprocess given at x= -1 (i.e sup x)
    for (int y= 0; y < processMesh.getNumCols(); y++) {
        if (!bc.check(processMesh.getNumRows()-1, y)) {
                flag4 = false;
            }
    }
    
    if (flag1){
        for (int y = 1; y<processMesh.getNumCols(); y++){
            for (int x = 0; x<processMesh.getNumRows();x++){
                std::pair<double, double> spaceTimeCoords = processMesh.getSpaceTimeMesh().getCoords(x, y);
                processMesh.setMeshData(x, y,  processMesh.getMeshData(x, y-1) + dt*dynamics.getDrift(spaceTimeCoords.second, spaceTimeCoords.first, processMesh.getMeshData(x, y-1)));
            }
        }
    }else if (flag2){
        for (int y = processMesh.getNumCols()-2; y>=0; y-- ){
            for (int x= 0; x<processMesh.getNumRows(); x++){
                std::pair<double, double> spaceTimeCoords = processMesh.getSpaceTimeMesh().getCoords(x, y);
                processMesh.setMeshData(x, y,  processMesh.getMeshData(x, y+1) - dt*dynamics.getDrift(spaceTimeCoords.second, spaceTimeCoords.first, processMesh.getMeshData(x, y+1)));
            }
        }
        
    }else if (flag3){
        for (int y =0; y<processMesh.getNumCols(); y++){
            for (int x= 1; x<processMesh.getNumRows(); x++){
                std::pair<double, double> spaceTimeCoords = processMesh.getSpaceTimeMesh().getCoords(x, y);
                processMesh.setMeshData(x, y,  processMesh.getMeshData(x-1, y) + dx*dynamics.getPseudoVol(spaceTimeCoords.second, spaceTimeCoords.first, processMesh.getMeshData(x-1, y)));
            }
        }
        
    }else if (flag4){
        for (int y =0; y<processMesh.getNumCols(); y++){
            for (int x= processMesh.getNumRows()-2; x>=0; x++){
                std::pair<double, double> spaceTimeCoords = processMesh.getSpaceTimeMesh().getCoords(x, y);
                processMesh.setMeshData(x, y,  processMesh.getMeshData(x+1, y) - dx*dynamics.getPseudoVol(spaceTimeCoords.second, spaceTimeCoords.first, processMesh.getMeshData(x+1, y)));
            }
        }
        
    }else{
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
                        dir.first * dynamics.getPseudoVol(spaceTimeCoords.second,spaceTimeCoords.first, processMesh.getMeshData(p.first,p.second))*dx +
                                                      dir.second * dynamics.getDrift(spaceTimeCoords.second,spaceTimeCoords.first, processMesh.getMeshData(p.first,p.second)) *dt;
                        processMesh.setMeshData(p.first + dir.first,p.second + dir.second, new_val);
                    }
                }
            }
            overflow_counter++;
        }
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

