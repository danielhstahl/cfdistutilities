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
This repository holds utilities for computing risk metrics like VaR and expected shortfall for any distribution with an analytic Characteristic Function.

In numerical tests it is extremely quick: for a portfolio of 1,000,000 loans the computation of VaR was faster than the computation of the graph of the density.  See my [FaaSDemo](https://github.com/phillyfan1138/ModelFaaSDemo).