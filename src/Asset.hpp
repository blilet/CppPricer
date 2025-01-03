#pragma once
#include "ItoProcess.hpp"
#include <functional>


class Asset {
private:
    double S0;
    ItoDynamics volDynamics;
    ItoDynamics rateDynamics;
    
public:
    Asset(double S0, const ItoDynamics& volDynamics, const ItoDynamics& rateDynamics);

    double getS0() const;
    const ItoDynamics& getVolDynamics() const;
    const ItoDynamics& getRateDynamics() const;

};

class Contract {
private:
    const Asset& underlying; // possibility of many contracts for same underlying
    std::function<double(double)> payoff;
    double T;
    
public:
    Contract(const Asset& underlying, const std::function<double(double)>& payoff, double maturity);
    const Asset& getUnderlying() const;
    const std::function<double(double)>& getPayoff() const;
    double getMaturity() const;


};
