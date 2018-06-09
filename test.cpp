#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "FunctionalUtilities.h"
#include <iostream>
#include "CFDistUtilities.h"
#include "CharacteristicFunctions.h"
#include "FangOost.h"
#include <complex>

/**NOTE That this is to test distributions for the option dashboard.  See issue
 * https://github.com/phillyfan1138/levy-functions/issues/27
 * 
 * */
auto cfLogBase(
    double T
){
    return [=](
        const auto& u,
        double lambda, 
        double muJ, double sigJ,
        double sigma, double v0, 
        double speed,double adaV, 
        double rho
    ){
        
        return chfunctions::cirLogMGF(
            -chfunctions::mertonLogRNCF(u, lambda, muJ, sigJ, 0.0, sigma),
            speed, 
            speed-adaV*rho*u*sigma,
            adaV,
            T,
            v0
        );
        
    };
}

auto cf(
    double r,
    double T
){
    //trivially copyable...this is the SV3 of the following paper:
    //https://pdfs.semanticscholar.org/67cd/b553e2624c79a960ff79d0dfe6e6833690a7.pdf 
    return [=](
        double lambda,
        double muJ, 
        double sigJ,
        double sigma,
        double v0,
        double speed,
        double adaV,
        double rho
    ){
        auto cfLogTmp=cfLogBase(T);
        return [=, cfLog=std::move(cfLogTmp)](const auto& u){
            return exp(r*T*u+
                cfLog(u, lambda, muJ, sigJ, sigma, v0, speed, adaV, rho)
            );
        };
    };
}
auto get_jump_diffusion_vol(double sigma, double lambda, double muJ, double sigJ, double T){
    return sqrt((sigma*sigma+lambda*(muJ*muJ+sigJ*sigJ))*T);
}
/**
Compile your application with -g, then you'll have debug symbols in the binary file.
Use gdb to open the gdb console.
Use file and pass it your application's binary file in the console.
Use run and pass in any arguments your application needs to start.
Do something to cause a Segmentation Fault.
Type bt in the gdb console to get a stack trace of the Segmentation Fault.*/
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
TEST_CASE("Test computeVaR Difficult Distribution", "[CFDistUtilities]"){
    const int numU=256;
    const double alpha=.01;
    const double r=.004;
    const double sigma=.3183;
    const double S0=191.96;
    const double sigJ=.220094;
    const double muJ=-.302967;
    const double lambda=.204516;
    const double speed=2.6726;
    const double v0=.237187;
    const double rho=-.182754;
    const double T=.187689;
    const double adaV=0;
    const double xMax=get_jump_diffusion_vol(sigma, lambda, muJ, sigJ, T)*5.0;
    const double xMin=-xMax;
    auto cfInst=cf(r, T)(lambda, muJ, sigJ, sigma, v0, speed, adaV, rho);
    //const auto qnormReference=6.224268;
    double prec=.0000001;
    auto myq=cfdistutilities::computeVaR(alpha, prec, xMin, xMax, numU, cfInst);
    auto myES=cfdistutilities::computeES(alpha, prec, xMin, xMax, numU, cfInst);
    std::cout<<"myq: "<<myq<<std::endl;
    //std::cout<<"myES: "<<std::get<cfdistutilities::ES>(myES)<<std::endl;
    REQUIRE(std::get<cfdistutilities::VAR>(myES)==Approx(myq));
} 
TEST_CASE("Test computeVaRNewton", "[CFDistUtilities]"){
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
    auto myqNorm=cfdistutilities::computeVaRNewton(alpha, prec, prec, xMin, xMax, mu, numU, normCF);
    REQUIRE(myqNorm==Approx(qnormReference));
} 
TEST_CASE("Test cdf", "[CFDistUtilities]"){
    const double mu=2;
    const double sigma=5;
    const int numU=64;
    const double xMin=-20;
    const double xMax=25;
    const double alpha=.05;
    auto normCF=[&](const auto& u){ //normal distribution's CF
        return exp(u*mu+.5*u*u*sigma*sigma);
    };      
    const auto pnormReference=.6554217; 
    const int xDiscrete=46;
    auto myPnorm=cfdistutilities::computeCDF(xDiscrete, numU, xMin, xMax, normCF);
    REQUIRE(myPnorm[24]==Approx(pnormReference));
} 
TEST_CASE("Test cdf discrete", "[CFDistUtilities]"){
    const double mu=2;
    const double sigma=5;
    const int numU=64;
    const double xMin=-20;
    const double xMax=25;
    const double alpha=.05;
    auto normCF=[&](const auto& u){ //normal distribution's CF
        return exp(u*mu+.5*u*u*sigma*sigma);
    };      
    const auto pnormReference=.6554217; 
    //corresponds with idndex 25
    const int xDiscrete=46;
    auto myPnorm=cfdistutilities::computeCDF(xDiscrete, xMin, xMax, fangoost::computeDiscreteCFReal(xMin, xMax, numU, normCF));
    REQUIRE(myPnorm[24]==Approx(pnormReference));
} 
TEST_CASE("Test cdf point", "[CFDistUtilities]"){
    const double mu=2;
    const double sigma=5;
    const int numU=64;
    const double xMin=-20;
    const double xMax=25;
    const double alpha=.05;
    auto normCF=[&](const auto& u){ //normal distribution's CF
        return exp(u*mu+.5*u*u*sigma*sigma);
    };      
    const auto pnormReference=.6554217; 
    //corresponds with idndex 25
   
    const auto actualResult=cfdistutilities::computeCDFAtPoint(4.0, numU, xMin, xMax, normCF);
    REQUIRE(actualResult==Approx(pnormReference));
} 



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
    auto myqNorm=std::get<cfdistutilities::ES>(cfdistutilities::computeES(alpha, prec, xMin, xMax, numU, normCF));
    REQUIRE(myqNorm==Approx(reference).epsilon(.0001)); 
} 
/*
TEST_CASE("Test computeESandVaRSeperate", "[CFDistUtilities]"){
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
    auto computeVar=cfdistutilities::computeVaR(alpha, prec, xMin, xMax, numU, normCF);
    auto myqNorm=std::get<0>(cfdistutilities::computeES(alpha,  xMin, xMax,computeVar, numU, normCF, true));
    //auto computeVar
    REQUIRE(myqNorm==Approx(reference).epsilon(.0001)); 
} */

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
    auto myqNorm=cfdistutilities::computeEL( xMin, xMax, numU, normCF);
    REQUIRE(myqNorm==Approx(mu)); 
} 