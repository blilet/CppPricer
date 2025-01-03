#include "Pricers.hpp"

double norm_cdf(double x) {
    return 0.5 * erfc(-x / sqrt(2));
}

DiscretePricer::DiscretePricer(int N, int N_T, const Contract& contract, double sigma_0, const BoundaryConditions& volBC,
                               const BoundaryConditions& rateBC, const BoundaryConditions& additionalBC, const SpaceTimeMesh& stm)
    : N(N), N_T(N_T), contract(contract), sigma_0(sigma_0),
      stm(stm),
       contractPrices(stm), volBC(volBC), rateBC(rateBC), additionalBC(additionalBC),
        current_theta(0.5), volApprox(stm), rateApprox(stm) {
        volApprox.solve(volBC, contract.getUnderlying().getVolDynamics());
        rateApprox.solve(volBC, contract.getUnderlying().getRateDynamics());
}


void DiscretePricer::price(double theta) {
    current_theta = theta;
    assert(theta <= 1 && theta >= 0);
    //double dx = stm.get_dx();
    //double dt = stm.get_dt();
    contractPrices.applyBoundaryConditions(additionalBC);
    for (std::size_t i = 0; i < stm.get_N(); i++) {
        contractPrices.setMeshData(i,stm.get_N_T() - 1, contract.getPayoff()(std::exp(stm.getCoords(i, stm.get_N_T() - 1).first)));
    }
    // TBD
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
