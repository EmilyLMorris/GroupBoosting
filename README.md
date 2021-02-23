
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GroupBoosting

<!-- badges: start -->

<!-- badges: end -->

Group boosting is an algorithm developed to perform variable selection
in a scalar-on-network regression setting, where the network has
underlying group structure.

## Installation

You can install the released version of GroupBoosting from
[CRAN](https://CRAN.R-project.org) with:

``` r
devtools::install_github("EmilyLMorris/GroupBoosting")
```

## Example

Here is a very simple example to show how to use the method:

``` r
library(GroupBoosting)

data1 = simul_dat_group(nodes = 6, n = 100, num.groups = 2, q.groups = 1, sparse_g = 0, dense_g = 1, effect_size = 5)
test1 = group_boosting_fit(data1$X, data1$Y, data1$group, total_steps=5000, step_size=1e-4, adj_var = 999, stop_tol=-1e-7, gamma = 0, lasso_lambda = 0.0314, weighted = FALSE)

cbind(test1$beta, data1$beta)
#>           [,1] [,2]
#>  [1,] 4.864109    5
#>  [2,] 5.330865    5
#>  [3,] 0.000000    0
#>  [4,] 0.000000    0
#>  [5,] 0.000000    0
#>  [6,] 5.203695    5
#>  [7,] 0.000000    0
#>  [8,] 0.000000    0
#>  [9,] 0.000000    0
#> [10,] 0.000000    0
#> [11,] 0.000000    0
#> [12,] 0.000000    0
#> [13,] 0.000000    0
#> [14,] 0.000000    0
#> [15,] 0.000000    0
```

If you have other adjustment variables to account for in addition to the
network this method, you can use the following example.

``` r

data2 = simul_dat_group(nodes = 6, n = 100, num.groups = 2, q.groups = 1, sparse_g = 0, dense_g = 1, effect_size = 5, adj_var = 3)
test2 = group_boosting_fit(data2$X, data2$Y, data2$group, total_steps=5000, step_size=1e-4, adj_var = c(1:3), stop_tol=-1e-7, gamma = 0, lasso_lambda = 0.0314, weighted = FALSE)

cbind(test2$beta, data2$beta[4:18])
#>           [,1] [,2]
#>  [1,] 5.850681    5
#>  [2,] 4.146541    5
#>  [3,] 0.000000    0
#>  [4,] 0.000000    0
#>  [5,] 0.000000    0
#>  [6,] 5.584251    5
#>  [7,] 0.000000    0
#>  [8,] 0.000000    0
#>  [9,] 0.000000    0
#> [10,] 0.000000    0
#> [11,] 0.000000    0
#> [12,] 0.000000    0
#> [13,] 0.000000    0
#> [14,] 0.000000    0
#> [15,] 0.000000    0
```
