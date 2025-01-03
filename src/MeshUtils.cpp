#include "MeshUtils.hpp"

template <typename T>
void logMatrix(const std::vector<std::vector<T>>& mat) {
    for (std::size_t i=0; i<mat[0].size();i++){
        std::cout<<"--";
    }
    std::cout<<"-> (axe des t)"<<std::endl;
    for (const auto& row : mat) {
        std::cout << "| ";
        for (const auto& e : row) {
            std::cout <<e << " ";
        }
        
        std::cout << '\n';
    }
    std::cout<<"v (axe des x)";
    std::cout<< "\n"<<std::endl;
}

BoundaryConditions::BoundaryConditions(const std::vector<std::vector<bool>>& contour, const std::function<double(double, double)>& function)
    : frontier(contour), frontier_function(function) {}

BoundaryConditions::BoundaryConditions(std::size_t X, std::size_t Y, const std::function<double(double, double)>& function)
    : frontier(X, std::vector<bool>(Y, false)), frontier_function(function) {}
BoundaryConditions::BoundaryConditions(const BoundaryConditions& other, const std::function<double(double, double)>& new_function)
        : frontier(other.frontier), frontier_function(new_function) {
    }
double BoundaryConditions::apply(double x, double y) const {
    return frontier_function(x, y);
}

bool BoundaryConditions::check(std::size_t x, std::size_t y) const {
    return frontier[x][y];
}
void BoundaryConditions::uncheck(std::size_t x, std::size_t y) {
    frontier[x][y] = false;
}
void BoundaryConditions::ToggleDir(bool dir, bool pos) {
    
    if (dir && pos) {for (std::size_t x = 0; x < frontier.size(); x++) {frontier[x][frontier[0].size() - 1] = true;}}
    else if (dir && !pos) {for (std::size_t x = 0; x < frontier.size(); x++) {frontier[x][0] = true;}}
    else if (!dir && pos) {for (std::size_t y = 0; y < frontier[0].size(); y++) {frontier[frontier.size() - 1][y] = true;}}
    else {for (std::size_t y = 0; y < frontier[0].size(); y++) {frontier[0][y] = true;}}
}

SpaceTimeMesh::SpaceTimeMesh(double x0, double R, double T, int N, int N_T): x0(x0), R(R), T(T), N(N), N_T(N_T) {
    assert( ((N & 1) == 1) && (N >= 2));
    assert(N_T >= 2);
}

std::size_t SpaceTimeMesh::get_N() const {
    return N;
}

std::size_t SpaceTimeMesh::get_N_T() const {
    return N_T;
}
double SpaceTimeMesh::get_T() const {
    return T;
}

double SpaceTimeMesh::get_R() const {
    return R;
}
double SpaceTimeMesh::get_dx() const{
    return R/(N-1);
}
double SpaceTimeMesh::get_dt() const{
    return T/(N_T-1);
}
std::pair<double, double> SpaceTimeMesh::getCoords(size_t i, size_t n) const {
    return {
        x0 + R * (static_cast<int>(i) - static_cast<int>(N) / 2) / static_cast<double>(N - 1),
        (static_cast<double>(n) / (N_T - 1)) * T
    };
}

FunctionMesh::FunctionMesh(const SpaceTimeMesh& stm): spaceTimeMesh(stm) {
    mesh_data = std::vector<std::vector<double>>(stm.get_N(), std::vector<double>(stm.get_N_T(), 0));
}

void FunctionMesh::applyBoundaryConditions(const BoundaryConditions& bc){
    for (int x = 0; x < mesh_data.size(); x++) {
        for (int y = 0; y < mesh_data[0].size(); y++) {
            if (bc.check(x, y)) {
                std::pair<double, double> coords = spaceTimeMesh.getCoords(x, y);
                mesh_data[x][y] = bc.apply(coords.first, coords.second);
            }
        }
    }
}
void FunctionMesh::logMesh() const{
    logMatrix(mesh_data);
}

void BoundaryConditions::logMesh() const{
    logMatrix(frontier);
}
void FunctionMesh::setMeshData(std::size_t i, std::size_t j, double val){
    mesh_data[i][j] = val;
}
double FunctionMesh::getMeshData(std::size_t i, std::size_t j) const{ return  mesh_data[i][j];}
std::size_t FunctionMesh::getNumRows() const{return mesh_data.size();}
std::size_t FunctionMesh::getNumCols() const{return mesh_data[0].size();}
const SpaceTimeMesh& FunctionMesh::getSpaceTimeMesh() const{ return spaceTimeMesh;}
