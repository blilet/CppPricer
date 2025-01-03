#pragma once
#include<iostream>
#include<vector>
#include<functional>
#include<cassert>

// for printing out the different meshes
template <typename T>
void logMatrix(const std::vector<std::vector<T>>& mat);

class BoundaryConditions {
private:
    // frontier_function is copied for stability
    std::vector<std::vector<bool>> frontier; // all vector matrices can be changed with one continuous vector (but for readability sake we opt for this approach)
    std::function<double(double, double)> frontier_function;
	
public:
    BoundaryConditions(const std::vector<std::vector<bool>>& contour, const std::function<double(double, double)>& function);
    BoundaryConditions(std::size_t X, std::size_t Y, const std::function<double(double, double)>& function);
    BoundaryConditions(const BoundaryConditions& other, const std::function<double(double, double)>& new_function); // necessary for vega calculus
    double apply(double x, double y) const;
    bool check(std::size_t x, std::size_t y) const;
    void uncheck(std::size_t x, std::size_t y);
    void ToggleDir(bool dir, bool pos);
    void logMesh() const;
};

class SpaceTimeMesh {
private:
    double x0;
    double R;
    double T;
    std::size_t N;
    std::size_t N_T;

public:
    SpaceTimeMesh(double x0, double R, double T, int N, int N_T);
    std::size_t get_N() const;
    std::size_t get_N_T() const;
    double get_T() const;
    double get_R() const;
    double get_dx() const;
    double get_dt() const;
    std::pair<double, double> getCoords(std::size_t i, std::size_t n) const;
};

class FunctionMesh {
    std::vector<std::vector<double>> mesh_data;
    const SpaceTimeMesh& spaceTimeMesh;
public:
    FunctionMesh(const SpaceTimeMesh& stm);
    void applyBoundaryConditions(const BoundaryConditions& bc);
    void logMesh() const;
    void setMeshData(std::size_t i, std::size_t j, double val);
    double getMeshData(std::size_t i, std::size_t j) const;
    std::size_t getNumRows() const;
    std::size_t getNumCols() const;
    const SpaceTimeMesh& getSpaceTimeMesh() const;
};
