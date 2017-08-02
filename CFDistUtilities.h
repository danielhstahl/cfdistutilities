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

namespace cfdistutilities {
    /**
        Function to compute the CDF of a distribution; see 
        http://danielhstahl.com/static/media/CreditRiskExtensions.c31991d2.pdf
        
    */
    template<typename Number, typename Index>
    auto VkCDF(const Number& u, const Number& x, const Number& a, const Number& b, const Index& k){
        return k==0?x-a:sin((x-a)*u)/u;
    }
    
    /**helper functions*/
    template<typename Number>
    auto powTwo(const Number& x){
        return x*x;
    }
    template<typename Number>
    auto diffPow(const Number& x, const Number& a){
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
    auto computeVaRDiscrete(const Number& alpha, const Number& xMin, const Number& xMax, CFDiscrete&& discreteCF, const Number& prec1, const Number& prec2){
        return -newton::bisect([&](const auto& pointInX){
            return fangoost::computeConvolutionAtPoint(pointInX, xMin, xMax, discreteCF, [&](const auto& u, const auto& x, const auto& index){
                return VkCDF(u, x, xMin, xMax, index);
            })-alpha;
        }, xMin, xMax, prec1, prec2);
    }

    template<typename Number, typename CF, typename Index>
    auto computeVaR(const Number& alpha, const Number& prec, const Number& xMin, const Number& xMax, const Index& numU, CF&& cf){
        return computeVaRDiscrete(alpha, xMin, xMax, fangoost::halfFirstIndex(fangoost::computeDiscreteCFReal(xMin, xMax, numU, cf)), prec, prec);
    }

    template<typename Number, typename CFDiscrete>
    auto computeVaRDiscrete(const Number& alpha, const Number& prec, const Number& xMin, const Number& xMax, CFDiscrete&& cf){
        return computeVaRDiscrete(alpha, xMin, xMax, fangoost::halfFirstIndex(cf), prec, prec);
    }

    /**note that in actual implementation we probably want 
    to return both VaR and ES from one function since they
    both get computed in this function*/
    template<typename Number, typename CF, typename Index>
    auto computeES(const Number& alpha, const Number& prec, const Number& xMin, const Number& xMax, const Index& numU, CF&& cf){
        auto computeVaRAndES=[&](auto&& discreteCF){
            return fangoost::computeConvolutionAtPoint(
                -computeVaRDiscrete(alpha, 
                    xMin, xMax, 
                    discreteCF, prec, prec
                ),
                xMin, 
                xMax,
                discreteCF, 
                [&](const auto& u, const auto& x, const auto& index){
                    return VkE(u, x, xMin, xMax, index);
                }
            );
        };
        return -computeVaRAndES(
            fangoost::halfFirstIndex(
                fangoost::computeDiscreteCFReal(xMin, xMax, numU, cf)
            )
        )/alpha;
    }
    template<typename Number, typename CFDiscrete>
    auto computeESDiscrete(const Number& alpha, const Number& prec, const Number& xMin, const Number& xMax, CFDiscrete&& cf){
        return fangoost::computeConvolutionAtPoint(
            -computeVaRDiscrete(alpha, prec,
                xMin, xMax, 
                cf
            ),
            xMin, 
            xMax,
            fangoost::halfFirstIndex(cf), 
            [&](const auto& u, const auto& x, const auto& index){
                return VkE(u, x, xMin, xMax, index);
            }
        )/alpha;
    }
    template<typename Number, typename CF, typename Index>
    auto computeEL(const Number& alpha, const Number& prec, const Number& xMin, const Number& xMax, const Index& numU, CF&& cf){
        return fangoost::computeExpectationPoint(xMax, xMin, xMax, numU, cf, [&](const auto& u, const auto& x, const auto& index){
            return VkE(u, x, xMin, xMax, index);
        });
    }


}


#endif