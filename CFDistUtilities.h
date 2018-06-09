#ifndef __CFDISTUTILITIES_H_INCLUDED__
#define __CFDISTUTILITIES_H_INCLUDED__
//https://github.com/phillyfan1138/cfdistutilities.git
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#include <cmath>
#include "FunctionalUtilities.h"
#include "FangOost.h"
#include "Newton.h"
#include <tuple>

namespace cfdistutilities {
    /**
        Function to compute the CDF of a distribution; see 
        http://danielhstahl.com/static/media/CreditRiskExtensions.c31991d2.pdf
        
    */
    template<typename X, typename Number, typename Index>
    auto VkCDF(const Number& u, const X& x, const Number& a, const Number& b, const Index& k){
        return k==0?x-a:sin((x-a)*u)/u;
    }
    
    /**helper functions*/
    template<typename Number>
    auto powTwo(const Number& x){
        return x*x;
    }
    template<typename Number, typename X>
    auto diffPow(const X& x, const Number& a){
        return .5*(powTwo(x)-powTwo(a));
    }
 
    /**
        Function to compute the partial expectation of a distribution; see
        http://danielhstahl.com/static/media/CreditRiskExtensions.c31991d2.pdf.
    */
    template<typename Number, typename Index>
    auto VkE(const Number& u, const Number& x, const Number& a, const Number& b, const Index& k){
        auto arg=(x-a)*u;
        auto uDen=1.0/u;
        return k==0?diffPow(x, a):x*sin(arg)*uDen+powTwo(uDen)*(cos(arg)-1.0);
    }
    //this is a helper function.  It bisects until it finds the point such that the CDF is equal to alpha 
    template<typename Number, typename CFDiscrete>
    auto computeVaRHelper(const Number& alpha, const Number& xMin, const Number& xMax, CFDiscrete&& discreteCF, const Number& prec1, const Number& prec2){

        return -newton::bisect([&](const auto& pointInX){

            return fangoost::computeExpectationPointDiscrete(pointInX, xMin, xMax, std::move(discreteCF), [&](const auto& u, const auto& x, const auto& index){
                return VkCDF(u, x, xMin, xMax, index);
            })-alpha;
        }, xMin, xMax, prec1, prec2);
    }

    template<typename Number, typename CFDiscrete, typename Index>
    auto computeCDF(const Index& numXDiscrete, const Number& xMin, const Number& xMax, CFDiscrete&& discreteCF){
        return fangoost::computeExpectationDiscrete(numXDiscrete, xMin, xMax, std::move(discreteCF), [&](const auto& u, const auto& x, const auto& index){
            return VkCDF(u, x, xMin, xMax, index);
        });
    }
    template<typename Number, typename CF, typename Index>
    auto computeCDF(const Index& numXDiscrete, const Index& numU, const Number& xMin, const Number& xMax, CF&& cf){
        return fangoost::computeExpectation(numXDiscrete, numU, xMin, xMax, std::move(cf), [&](const auto& u, const auto& x, const auto& index){
            return VkCDF(u, x, xMin, xMax, index);
        });
    }
    template<typename Number, typename CFDiscrete>
    auto computeCDFAtPoint(const Number& xValue, const Number& xMin, const Number&xMax, CFDiscrete&& cfDiscrete){
        return fangoost::computeExpectationPointDiscrete(xValue, xMin, xMax, std::move(cfDiscrete), [&](const auto& u, const auto& x, const auto& index){
            return VkCDF(u, x, xMin, xMax, index);
        });
    }
    template<typename Number, typename CF, typename Index>
    auto computeCDFAtPoint(const Number& xValue,const Index& numU, const Number& xMin, const Number&xMax, CF&& cf){
        return fangoost::computeExpectationPoint(xValue, xMin, xMax, numU, std::move(cf), [&](const auto& u, const auto& x, const auto& index){
            return VkCDF(u, x, xMin, xMax, index);
        });
    }


    template<typename Number, typename CFDiscrete>
    auto computeVaRNewtonHelper(const Number& alpha, const Number& xMin, const Number& xMax, const Number& guess, CFDiscrete&& discreteCF, const Number& prec1, const Number& prec2){
        return -newton::zeros([&](const auto& pointInX){
            return fangoost::computeExpectationPointDiscrete(pointInX, xMin, xMax, std::move(discreteCF), [&](const auto& u, const auto& x, const auto& index){
                return VkCDF(u, x, xMin, xMax, index);
            })-alpha;
        }, guess, prec1, prec2, 50);
    }

    template<typename Number, typename CF, typename Index>
    auto computeVaR(const Number& alpha, const Number& prec, const Number& xMin, const Number& xMax, const Index& numU, CF&& cf){
        return computeVaRHelper(alpha, xMin, xMax, fangoost::computeDiscreteCFReal(xMin, xMax, numU, std::move(cf)), prec, prec);
    }
    /**Newton is faster but not as stable*/
    template<typename Number, typename CF, typename Index>
    auto computeVaRNewton(const Number& alpha, const Number& prec1, const Number& prec2, const Number& xMin, const Number& xMax, const Number& guess, const Index& numU, CF&& cf){
        return computeVaRNewtonHelper(alpha, xMin, xMax, guess, fangoost::computeDiscreteCFReal(xMin, xMax, numU, std::move(cf)), prec1, prec2);
    }

    template<typename Number, typename CFDiscrete>
    auto computeVaRNewtonDiscrete(const Number& alpha, const Number& prec1, const Number& prec2, const Number& xMin, const Number& xMax, const Number& guess, CFDiscrete&& cf){
        return computeVaRNewtonHelper(alpha, xMin, xMax, guess, std::move(cf), prec1, prec2);
    }

    template<typename Number, typename CFDiscrete>
    auto computeVaRDiscrete(const Number& alpha, const Number& prec, const Number& xMin, const Number& xMax, CFDiscrete&& cf){
        return computeVaRHelper(alpha, xMin, xMax, std::move(cf), prec, prec);
    }

    constexpr int ES=0;
    constexpr int VAR=1;
    /**
     * returns tuple of ES and VaR
     */
    template<typename Number,typename CFDiscrete>
    auto computeESDiscrete(const Number& alpha, const Number& prec, const Number& xMin, const Number& xMax, CFDiscrete&& cf){
        const auto VaR=computeVaRDiscrete(alpha, prec,
            xMin, xMax, 
            cf
        );
        return std::make_tuple(
            -fangoost::computeExpectationPointDiscrete(
            -VaR,
            xMin, 
            xMax,
            std::move(cf), 
            [&](const auto& u, const auto& x, const auto& index){
                return VkE(u, x, xMin, xMax, index);
            }
        )/alpha, VaR);
    }
    /**
     * returns tuple of ES and VaR
     */
    template<typename Number, typename CF, typename Index>
    auto computeES(const Number& alpha, const Number& prec, const Number& xMin, const Number& xMax, const Index& numU, CF&& cf){
        return computeESDiscrete(alpha, prec, xMin, xMax, 
            fangoost::computeDiscreteCFReal(xMin, xMax, numU, std::move(cf))
        );
    }
    
    template<typename Number, typename CF, typename Index>
    auto computeEL(const Number& xMin, const Number& xMax, const Index& numU, CF&& cf){
        return fangoost::computeExpectationPoint(xMax, xMin, xMax, numU, std::move(cf), [&](const auto& u, const auto& x, const auto& index){
            return VkE(u, x, xMin, xMax, index);
        });
    }
    template<typename Number, typename CFDiscrete>
    auto computeELDiscrete(const Number& xMin, const Number& xMax, CFDiscrete&& cf){
        return fangoost::computeExpectationPointDiscrete(xMax, xMin, xMax, std::move(cf), [&](const auto& u, const auto& x, const auto& index){
            return VkE(u, x, xMin, xMax, index);
        });
    }


}


#endif