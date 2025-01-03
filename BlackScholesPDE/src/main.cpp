#include <iostream>
#include "ItoProcess.hpp"
#include "Asset.hpp"
#include "Pricers.hpp"
int main(int argc, const char * argv[]) {
    
    // Example file (replace main) for pricing european calls and comparing with BS closed form model
    // For simplicty's sake, and for a fair comparaison, we assume constant rate/vol
    
    // vol ito dynamics and boundary conditions
    std::function<double(double, double, double)> volDrift = [](double t, double x, double p) { return 0; }; // cste vol
    std::function<double(double, double, double)> volStochastic = [](double t, double x, double p) { return 0; }; //cste vol
    
    // rate ito dynamics
    std::function<double(double, double, double)> rateDrift = [](double t, double x, double p) { return 0; }; // cste rate
    std::function<double(double, double, double)> rateStochastic = [](double t, double x, double p) { return 0; }; //cste rate
    
    // init of constant ItoDynamics
    ItoDynamics volDynamics(volDrift, volStochastic);
    ItoDynamics  rateDynamics(rateDrift, rateStochastic);
    
    
    // init of underlying (we take AAPL as an example)
    double S0 = 250;
    Asset underlying(S0, volDynamics,rateDynamics); // cste dynamics

    
    // init of a contract (European Call for example)
    int K = 250;
    double T = 31; // month contract
    std::function<double(double)> payoff = [&K](double S) { return std::abs(static_cast<int>(S>=K)* (S-K)); }; // european call payoff
    Contract contract(underlying, payoff, T);
    
    // pricing and creating mesh
    
    int N = 101; //nbr of steps in x
    int N_T = 31; //nbr of steps in t (we take as delta T = 1 day)
    double sigma_0 = 0.19; // cste vol sigma_0
    double r_0 = 0.05; // cste risk-free rate
    
    // defining vol and rate boundary conditions (in this case we suppose sigma(0,.) as well as r(0,.) are given and are constant)
    std::function<double(double, double)> csteVol = [&sigma_0](double t, double x) { return sigma_0; };
    std::function<double(double, double)> csteRate = [&r_0](double t, double x) { return r_0; };
    
    BoundaryConditions volBoundaries(N, N_T, csteVol);
    BoundaryConditions rateBoundaries(N, N_T, csteRate);
    volBoundaries.ToggleDir(true, false); // this is the section t = 0 ...
    rateBoundaries.ToggleDir(true, false);

    // defining additional contract boundaries (not payoff)
    std::function<double(double, double)> zeroPayoff= [](double t, double x) { return 0; };
    BoundaryConditions contractAdditionalBoundaries(N,N_T, zeroPayoff);
    contractAdditionalBoundaries.ToggleDir(false, false);
    
    // init mesh
    SpaceTimeMesh stm(std::log(contract.getUnderlying().getS0()), 5*sigma_0 * std::sqrt(contract.getMaturity()), contract.getMaturity(), N, N_T);
    
    // init pricer
    DiscretePricer pricer(N, N_T, contract, sigma_0, volBoundaries, rateBoundaries, contractAdditionalBoundaries,stm);
    pricer.price(0.5);
    //pricer.logMesh();  //uncomment to visualize the pricer mesh
    
    std::cout << "Crank-Nicholson Scheme's price: "<< pricer.getPrice()<<std::endl;
    std::cout << "Crank-Nicholson Scheme's delta: "<< pricer.delta()<<std::endl;
    std::cout << "Crank-Nicholson Scheme's gamma: "<< pricer.gamma()<<std::endl;
    std::cout << "Crank-Nicholson Scheme's theta: "<< pricer.theta()<<std::endl;
    std::cout << "Crank-Nicholson Scheme's vega: "<< pricer.vega()<<std::endl;
    std::cout<<std::endl;
    
    // compare w/ BlackScholes
    BlackScholesCallPricer bsPricer(S0, K, T, r_0, sigma_0);
    bsPricer.price();
    std::cout<< "Black Scholes closed form's price: " << bsPricer.getPrice()<<std::endl;
    std::cout<< "Black Scholes closed form's delta: " << bsPricer.delta()<<std::endl;
    std::cout<< "Black Scholes closed form's gamma: " << bsPricer.gamma()<<std::endl;
    std::cout<< "Black Scholes closed form's theta: " << bsPricer.theta()<<std::endl;
    std::cout<< "Black Scholes closed form's vega: " << bsPricer.vega()<<std::endl;
    
    return 0;
}
    

