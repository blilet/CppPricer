#include "Asset.hpp"

Asset::Asset(double S0, const ItoDynamics& volDynamics, const ItoDynamics& rateDynamics)
    : S0(S0), volDynamics(volDynamics), rateDynamics(rateDynamics) {}

double Asset::getS0() const {
    return S0;
}

const ItoDynamics& Asset::getVolDynamics() const {
    return volDynamics;
}

const ItoDynamics& Asset::getRateDynamics() const {
    return rateDynamics;
}

Contract::Contract(const Asset& underlying, const std::function<double(double)>& payoff, double maturity)
    : underlying(underlying), payoff(std::move(payoff)), T(maturity) {}

const Asset& Contract::getUnderlying() const {
    return underlying;
}

const std::function<double(double)>& Contract::getPayoff() const {
    return payoff;
}

double Contract::getMaturity() const {
    return T;
}
