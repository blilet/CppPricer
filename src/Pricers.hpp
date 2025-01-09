#pragma once
#include "Asset.hpp"
#include "MeshUtils.hpp"

double norm_cdf(double x);
std::pair<long double,long double> solve_Mx_b(long double& A, long double& B, long double& C, long double& D, long double& E, long double& F);
class DiscretePricer {
private:
    int N;
    int N_T;
    double sigma_0;
    const Contract& contract;
    const BoundaryConditions& volBC;
    const BoundaryConditions& rateBC;
    
    double current_theta;
    ItoProcess volApprox;
    ItoProcess rateApprox;
    const SpaceTimeMesh& stm;
    FunctionMesh contractPrices;

public:
    const BoundaryConditions& additionalBC;
    DiscretePricer(int N, int N_T, const Contract& contract, double sigma_0, const BoundaryConditions& volBC,
                   const BoundaryConditions& driftBC, const BoundaryConditions& additionalBC, const SpaceTimeMesh& stm);

    void price(double theta);
    const ItoProcess& getVolApprox() const;
    const ItoProcess& getRateApprox() const;
    double getPrice();
    double delta();
    double gamma();
    double theta();
    double vega(double d_sigma =10e-3);
    void logMesh();
};

class BlackScholesCallPricer {
private:
    double S0;
    double K;
    double T;
    double r;
    double sigma;
    double contractPrice;

public:
    BlackScholesCallPricer(double S0, double K, double T, double r, double sigma);

    void price();
    double getPrice();
    double delta();
    double gamma();
    double theta();
    double vega();
    
};

