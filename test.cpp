#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "FunctionalUtilities.h"
#include <iostream>
#include "CFDistUtilities.h"
#include "FangOost.h"
#include <complex>

TEST_CASE("Test computeVaR", "[CFDistUtilities]"){
    const double mu=2;
    const double sigma=5;
    const int numU=64;
    const double xMin=-20;
    const double xMax=25;
    const double alpha=.05;
    auto normCF=[&](const auto& u){ //normal distribution's CF
        return exp(u*mu+.5*u*u*sigma*sigma);
    };      
    const auto qnormReference=6.224268;
    double prec=.0000001;
    auto myqNorm=cfdistutilities::computeVaR(alpha, prec, xMin, xMax, numU, normCF);
    REQUIRE(myqNorm==Approx(qnormReference));
} 
/*
TEST_CASE("Test computeES", "[CFDistUtilities]"){
    const double mu=2;
    const double sigma=1;
    const int numX=5;
    const int numU=64;
    const double xMin=-3;
    const double xMax=7;
    const double alpha=.05;
    auto normCF=[&](const auto& u){ //normal distribution's CF
        return exp(u*mu+.5*u*u*sigma*sigma);
    };      
    const auto reference=-0.06271281;
    double prec=.00001;
    auto myqNorm=cfdistutilities::computeES(.3551464, xMin, xMax, numU, normCF);
    std::cout<<myqNorm<<std::endl;
    for(int i=0; i<numX; ++i){
        REQUIRE(myqNorm==Approx(reference));
    }  
} */
TEST_CASE("Test computeESandVaR", "[CFDistUtilities]"){
    const double mu=2;
    const double sigma=5;
    const int numU=64;
    const double xMin=-20;
    const double xMax=25;
    const double alpha=.05;
    auto normCF=[&](const auto& u){ //normal distribution's CF
        return exp(u*mu+.5*u*u*sigma*sigma);
    };      
    //const auto reference=-0.06271281;
    const auto reference=8.313564;
    double prec=.0000000001;
    auto myqNorm=cfdistutilities::computeES(alpha, prec, xMin, xMax, numU, normCF);
    REQUIRE(myqNorm==Approx(reference).epsilon(.0001)); 
} 
TEST_CASE("Test computeEL", "[CFDistUtilities]"){
    const double mu=2;
    const double sigma=5;
    const int numU=64;
    const double xMin=-20;
    const double xMax=25;
    const double alpha=.05;
    auto normCF=[&](const auto& u){ //normal distribution's CF
        return exp(u*mu+.5*u*u*sigma*sigma);
    };      
    //const auto reference=-0.06271281;
    double prec=.0000000001;
    auto myqNorm=cfdistutilities::computeEL(alpha, prec, xMin, xMax, numU, normCF);
    REQUIRE(myqNorm==Approx(mu)); 
} 