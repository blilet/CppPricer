#include <iostream>
#include <algorithm>
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
    
    
    // init of underlying
    double S0 = 200;
    Asset underlying(S0, volDynamics,rateDynamics); // cste dynamics

    
    // init of a contract (European Call for example)
    int K = 250;
    double T = 1; // month contract
    std::function<double(double)> payoff = [&K](double S) { return std::max(S-K,0.); }; // european call payoff
    Contract contract(underlying, payoff, T);
    
    // pricing and creating mesh
    
    int N = 1001; //nbr of steps in x
    int N_T = 1000; //nbr of steps in t (we take as delta T = 1 day)
    double sigma_0 = .05; // cste vol sigma_0
    double r_0 = 0.2; // cste risk-free rate
    
    // defining vol and rate boundary conditions (in this case we suppose sigma(0,.) as well as r(0,.) are given and are constant)
    std::function<double(double, double)> csteVol = [&sigma_0](double t, double x) { return sigma_0; };
    std::function<double(double, double)> csteRate = [&r_0](double t, double x) { return r_0; };
    
    BoundaryConditions volBoundaries(N, N_T, csteVol);
    BoundaryConditions rateBoundaries(N, N_T, csteRate);
    volBoundaries.ToggleDir(true, false); // this is the section t = 0 ...
    rateBoundaries.ToggleDir(true, false);
    // volBoundaries.logMesh();//uncomment to visualize

    // defining additional contract boundaries (not payoff)
    std::function<double(double, double)> zeroPayoff= [&T](double t, double x) { return 0; };
    std::function<double(double, double)> bsBoundaries= [&](double t, double x) {
        BlackScholesCallPricer bsP(std::exp(x), K, T-t, r_0, sigma_0);
        bsP.price();
        return bsP.getPrice();
    };
    BoundaryConditions contractAdditionalBoundaries(N,N_T, bsBoundaries);
    contractAdditionalBoundaries.ToggleDir(false, false);
    // contractAdditionalBoundaries.logMesh(); // f_0 uncomment to log
    
    // init mesh
    SpaceTimeMesh stm(std::log(contract.getUnderlying().getS0()), 5*sigma_0 * std::sqrt(contract.getMaturity()), contract.getMaturity(), N, N_T);
    // init pricer
    std::cout << "dx: " <<stm.get_dx()<<std::endl;
    std::cout << "dt: " <<stm.get_dt()<<std::endl;

    

    for (int i =0; i<10; i++){
        DiscretePricer pricer(N, N_T, contract, sigma_0, volBoundaries, rateBoundaries, contractAdditionalBoundaries,stm);
        pricer.price(static_cast<double>(i)/10);
        std::cout << "Crank-Nicholson Scheme's price for theta = " <<static_cast<double>(i)/10<<" : "<< pricer.getPrice()<<std::endl;
    }
    std::cout <<std::endl;
    DiscretePricer pricer(N, N_T, contract, sigma_0, volBoundaries, rateBoundaries, contractAdditionalBoundaries,stm);
    pricer.price(0.5);
    //pricer.logMesh();  //uncomment to visualize the pricer mesh
    std::cout << "Crank-Nicholson Scheme's price for theta = 0.5 : "<< pricer.getPrice()<<std::endl;
    std::cout << "Crank-Nicholson Scheme's delta for theta = 0.5 : "<< pricer.delta()<<std::endl;
    std::cout << "Crank-Nicholson Scheme's gamma  for theta = 0.5 : "<< pricer.gamma()<<std::endl;
    std::cout << "Crank-Nicholson Scheme's theta for theta (scheme) = 0.5 : "<< pricer.theta()<<std::endl;
    std::cout << "Crank-Nicholson Scheme's vega for theta = 0.5 : "<< pricer.vega()<<std::endl;
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
    

