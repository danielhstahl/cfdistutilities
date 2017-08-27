| [Linux][lin-link] | [Windows][win-link] | [Codecov][cov-link] |
| :---------------: | :-----------------: | :-------------------: |
| ![lin-badge]      | ![win-badge]        | ![cov-badge]          |

[lin-badge]: https://travis-ci.org/phillyfan1138/cfdistutilities.svg?branch=master "Travis build status"
[lin-link]:  https://travis-ci.org/phillyfan1138/cfdistutilities "Travis build status"
[win-badge]: https://ci.appveyor.com/api/projects/status/y36u1hdyxjj9r2a0?svg=true "AppVeyor build status"
[win-link]:  https://ci.appveyor.com/project/phillyfan1138/cfdistutilities "AppVeyor build status"
[cov-badge]: https://codecov.io/gh/phillyfan1138/cfdistutilities/branch/master/graph/badge.svg
[cov-link]:  https://codecov.io/gh/phillyfan1138/cfdistutilities

# CFDistUtilities
This repository holds utilities for computing risk metrics like VaR and expected shortfall for any distribution with an analytic Characteristic Function.  See the [tests](./test.cpp) for examples.

In numerical tests it is extremely quick: for a portfolio of 1,000,000 loans the computation of VaR was faster than the computation of the graph of the density.  See my [FaaSDemo](https://github.com/phillyfan1138/ModelFaaSDemo).

## Potential limitations

* For densities without derivatives of all orders, the convergence may be slow.  For example, Beta distributions may not converge at all when the mode of the distribution is near zero or one.  
* The VaR technique works when the cumulative density is monotonic.  While this is manifestly the case for any density, it is not the case that the cosine approximation is monotonic.  Hence the VaR may not converge to the actual VaR.  However, in tests it appears that it does; at least for "nice" distributions.  
* The Newton VaR technique is faster than the bisection technique; however for large domains the derivative is tiny and Newton's algorithm diverges.  The bisect method is safer.