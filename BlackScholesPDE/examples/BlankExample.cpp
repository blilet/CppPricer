#include <iostream>
#include "ItoProcess.hpp"
#include "Asset.hpp"
#include "Pricers.hpp"
int blankExample(int argc, const char * argv[]) {

    // v: volatily's ito dynamics and boundary conditions
    // we suppose dv = volDrift(t,x,v) dt+ volStochastic(t,x,v)dx (attention: x not S, you can use _to_x in ItoProcess)
    
    std::function<double(double, double, double)> volDrift; // ! INIT
    std::function<double(double, double, double)> volStochastic; // ! INIT
    
    // r: rate ito dynamics
    // we suppose dr = rateDrift(t,x,r) dt+ rateStochastic(t,x,r) dx
    
    std::function<double(double, double, double)> rateDrift; // ! INIT
    std::function<double(double, double, double)> rateStochastic; // ! INIT
    
    // init of constant ItoDynamics
    ItoDynamics volDynamics(volDrift, volStochastic);
    ItoDynamics  rateDynamics(rateDrift, rateStochastic);
    
    
    // init of underlying
    double S0; // ! INIT
    Asset underlying(S0, volDynamics,rateDynamics);

    
    // init of a contract (European Call for example)
    double T ; // ! INIT
    std::function<double(double)> payoff; // ! INIT
    Contract contract(underlying, payoff, T);
    
    // pricing and creating mesh
    
    int N ; // ! INIT : nbr of steps in x
    int N_T = 31; // ! INIT: nbr of steps in t (we take as delta T = 1 day)
    double sigma_0; // ! INIT : initial
    
    // defining vol and rate boundary conditions (in this case we suppose sigma(0,.) as well as r(0,.) are given and are constant)
    std::function<double(double, double)> boundaryVolFunction ; // ! INIT
    std::function<double(double, double)> boundaryRateFunction; // ! INIT
    
    BoundaryConditions volBoundaries(N, N_T, boundaryVolFunction);
    BoundaryConditions rateBoundaries(N, N_T, boundaryRateFunction);
    volBoundaries.ToggleDir(true, false); // this is the section t = 0 ...
    rateBoundaries.ToggleDir(true, false);

    // defining additional contract boundaries (not payoff)
    std::function<double(double, double)> additionalFunc; // ! INIT
     BoundaryConditions contractAdditionalBoundaries(N,N_T, additionalFunc);
    contractAdditionalBoundaries.ToggleDir(false, false); // this is x= x_min, ...
    
    // init mesh
    SpaceTimeMesh stm(std::log(contract.getUnderlying().getS0()), 5*sigma_0 * std::sqrt(contract.getMaturity()), contract.getMaturity(), N, N_T);
    
    // init pricer
    DiscretePricer pricer(N, N_T, contract, sigma_0, volBoundaries, rateBoundaries, contractAdditionalBoundaries,stm);
    double theta; // ! INIT
    pricer.price(theta);
    
    //pricer.logMesh();  //uncomment to visualize the pricer mesh
    
    std::cout << "Crank-Nicholson Scheme's price: "<< pricer.getPrice()<<std::endl;
    std::cout << "Crank-Nicholson Scheme's delta: "<< pricer.delta()<<std::endl;
    std::cout << "Crank-Nicholson Scheme's gamma: "<< pricer.gamma()<<std::endl;
    std::cout << "Crank-Nicholson Scheme's theta: "<< pricer.theta()<<std::endl;
    std::cout << "Crank-Nicholson Scheme's vega: "<< pricer.vega()<<std::endl;
    std::cout<<std::endl;
    
    return 0;
    
}
