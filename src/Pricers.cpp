#include "Pricers.hpp"

double norm_cdf(double x) {
    return 0.5 * erfc(-x / sqrt(2));
}
std::pair<long double, long double> solve_Mx_b(long double& A, long double& B, long double& C, long double& D, long double& E, long double& F) {
    long double det = A * E - B * D;
    if (det == 0) {
        throw std::invalid_argument("det nul, system non détérminé.");
    }
    long double x = (C * E - B * F) / det;
    long double y = (A * F - C * D) / det;

    return {x, y};
}

DiscretePricer::DiscretePricer(int N, int N_T, const Contract& contract, double sigma_0, const BoundaryConditions& volBC,
                               const BoundaryConditions& rateBC, const BoundaryConditions& additionalBC, const SpaceTimeMesh& stm)
    : N(N), N_T(N_T), contract(contract), sigma_0(sigma_0),
      stm(stm),
       contractPrices(stm), volBC(volBC), rateBC(rateBC), additionalBC(additionalBC),
        current_theta(0.5), volApprox(stm), rateApprox(stm) {
            assert(stm.get_N() == N);
            assert(stm.get_N_T() == N_T);
        volApprox.solve(volBC, contract.getUnderlying().getVolDynamics());
        rateApprox.solve(rateBC, contract.getUnderlying().getRateDynamics());
        // boundary conditions (not payoff), in our case we suppose that it is x = x_0 = inf_{x_r\in mesh} x_r
        contractPrices.applyBoundaryConditions(additionalBC);
        // other boundary condition which is the payoff hence f_0 and f^T are supposed available
        for (std::size_t i = 0; i < N; i++) {
            if (i == N-1 && std::abs(contractPrices.getMeshData(i, stm.get_N_T()-1) - contract.getPayoff()(std::exp(stm.getCoords(i, stm.get_N_T() - 1).first)))<10e-3){
                std::cerr << "Warning, boundary condition x = inf x (f_0) and t=T (payoff) don't coincide, we take the value of payoff..."<<std::endl;
            }
            contractPrices.setMeshData(i,stm.get_N_T() - 1, contract.getPayoff()(std::exp(stm.getCoords(i, stm.get_N_T() - 1).first)));
        }
            
}


void DiscretePricer::price(double theta) {
    current_theta = theta;
    assert(theta <= 1 && theta >= 0);
    double dx = stm.get_dx();
    double dt = stm.get_dt();
    
    
    // det of f_1^T-1 , f_2^T-1 and subsequent f_1, f_2
    
    for (int n = static_cast<int> (stm.get_N_T() - 2); n>=0; n--){ // int cause size_t -- >=0 gets stuck at 0
        long double a = -0.5*std::pow(volApprox.getVal(1, n),2)/(dx*dx) + 0.25*std::pow(volApprox.getVal(1, n),2)/(dx) - 0.5*rateApprox.getVal(1, n)/dx;
        long double b = rateApprox.getVal(1, n)+std::pow(volApprox.getVal(1, n),2)/(dx*dx);
        long double c = -0.5*std::pow(volApprox.getVal(1, n),2)/(dx*dx) - 0.25*std::pow(volApprox.getVal(1, n),2)/(dx) + 0.5*rateApprox.getVal(1, n)/dx;
        long double A= a*theta;
        long double B =(b*theta+1/dt);
        long double C = contractPrices.getMeshData(1, n+1)/dt - c*contractPrices.getMeshData(0, n) -
        (1-theta)*( a*contractPrices.getMeshData(2, n+1) + b*contractPrices.getMeshData(1, n+1) +c*contractPrices.getMeshData(0, n+1));
        long double D = 0.5*std::pow(volApprox.getVal(1, n),2);
        long double E = -(std::pow(volApprox.getVal(1, n),2)/(dx*dx) +(std::pow(volApprox.getVal(1, n),2) - rateApprox.getVal(1, n))/dx);
        long double F = - (contractPrices.getMeshData(0, n+1) - contractPrices.getMeshData(0, n))/dt +(
        rateApprox.getVal(1, n) - (0.5*std::pow(volApprox.getVal(1, n),2)/(dx*dx) + (0.5* std::pow(volApprox.getVal(1, n),2)-rateApprox.getVal(1, n))/dx))*
        rateApprox.getVal(1, n);
        auto p  = solve_Mx_b(A, B, C, D, E, F);
        // std::cout<< n <<"  :  "<<C<< " ,  " << contractPrices.getMeshData(1, n+1) <<"   , " << contractPrices.getMeshData(2, n+1)<< std::endl;
        contractPrices.setMeshData(2,n , p.first);
        contractPrices.setMeshData(1,n , p.second);
    }
    
    // det of rest
    for (size_t i = 3; i<stm.get_N()-1; i++){
        for (int n = static_cast<int> (stm.get_N_T() - 2); n>=0; n--){ // same
            long double a = -0.5*std::pow(volApprox.getVal(i, n),2)/(dx*dx) + 0.25*std::pow(volApprox.getVal(i, n),2)/(dx) - 0.5*rateApprox.getVal(i, n)/dx;
            long double b = rateApprox.getVal(i, n)+std::pow(volApprox.getVal(i, n),2)/(dx*dx);
            long double c = -0.5*std::pow(volApprox.getVal(i, n),2)/(dx*dx) - 0.25*std::pow(volApprox.getVal(i, n),2)/(dx) + 0.5*rateApprox.getVal(i, n)/dx;
            long double A= a*theta;
            long double B =(b*theta+1/dt);
            long double C = contractPrices.getMeshData(i, n+1)/dt - c*contractPrices.getMeshData(i, n) -
            (1-theta)*( a*contractPrices.getMeshData(i+1, n+1) + b*contractPrices.getMeshData(i, n+1) +c*contractPrices.getMeshData(i-1, n+1));
            long double D = 0.5*std::pow(volApprox.getVal(i, n),2);
            long double E = -(std::pow(volApprox.getVal(i, n),2)/(dx*dx) +(std::pow(volApprox.getVal(i, n),2) - rateApprox.getVal(i, n))/dx);
            long double F = - (contractPrices.getMeshData(i-1, n+1) - contractPrices.getMeshData(i-1, n))/dt +(
            rateApprox.getVal(1, n) - (0.5*std::pow(volApprox.getVal(i, n),2)/(dx*dx) + (0.5* std::pow(volApprox.getVal(i, n),2)-rateApprox.getVal(i, n))/dx))*
            rateApprox.getVal(i, n);
            auto p  = solve_Mx_b(A, B, C, D, E, F);
            // std::cout<< n <<"  :  "<<C<< " ,  " << contractPrices.getMeshData(1, n+1) <<"   , " << contractPrices.getMeshData(2, n+1)<< std::endl;
            contractPrices.setMeshData(i,n , p.second);
            if (i == stm.get_N() - 2){
                contractPrices.setMeshData(i+1,n ,p.first);
            }
        }
    }
}
const ItoProcess& DiscretePricer::getVolApprox() const{
    return volApprox;
}
const ItoProcess& DiscretePricer::getRateApprox() const{
    return rateApprox;
}
double DiscretePricer::getPrice(){
    return contractPrices.getMeshData(stm.get_N() / 2,0);
}
double DiscretePricer::delta() {
    return (contractPrices.getMeshData(stm.get_N() / 2 + 1,0) - contractPrices.getMeshData(stm.get_N() / 2 - 1,0)) / (2*stm.get_dx());
}

double DiscretePricer::gamma() {
    return (contractPrices.getMeshData(stm.get_N() / 2 + 1,0) + contractPrices.getMeshData(stm.get_N() / 2 - 1,0) - 2 * contractPrices.getMeshData(stm.get_N() / 2 ,0)) / (stm.get_dx() * stm.get_dx());
}

double DiscretePricer::theta() {
    return (contractPrices.getMeshData(stm.get_N() / 2 ,0) - contractPrices.getMeshData(stm.get_N() / 2 ,1)) / stm.get_dt();
}

double DiscretePricer::vega(double d_sigma) {
    d_sigma = sigma_0/100;
    std::function<double(double, double)> function_perturbed = [this, d_sigma](double x, double y) {
        return volBC.apply(x, y) + d_sigma;
    };
    BoundaryConditions volBC_perturbed(volBC, function_perturbed);
    DiscretePricer perturbed(N, N_T, contract, sigma_0, volBC_perturbed, rateBC, additionalBC,stm);
    perturbed.price(current_theta);
    return (perturbed.getPrice() - this->getPrice()) / d_sigma;
}

BlackScholesCallPricer::BlackScholesCallPricer(double S0, double K, double T, double r, double sigma)
    : S0(S0), K(K), T(T), r(r), sigma(sigma), contractPrice(0) {}

void BlackScholesCallPricer::price() {
    double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);

    contractPrice = S0 * norm_cdf(d1) - K * exp(-r * T) * norm_cdf(d2);
}
double BlackScholesCallPricer::getPrice(){
    return contractPrice;
}
double BlackScholesCallPricer::delta() {
    double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    return norm_cdf(d1);
}

double BlackScholesCallPricer::gamma() {
    double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    return exp(-0.5 * d1 * d1) / (S0 * sigma * sqrt(T) * sqrt(2 * M_PI));
}

double BlackScholesCallPricer::theta() {
    double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);
    double term1 = -S0 * exp(-0.5 * d1 * d1) * sigma / (2 * sqrt(T) * sqrt(2 * M_PI));
    double term2 = r * K * exp(-r * T) * norm_cdf(d2);
    return term1 - term2;
}

double BlackScholesCallPricer::vega() {
    double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrt(T));
    return S0 * sqrt(T) * exp(-0.5 * d1 * d1) / sqrt(2 * M_PI);
}
void DiscretePricer::logMesh(){
    contractPrices.logMesh();
}
