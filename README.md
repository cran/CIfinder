
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CIfinder

<!-- badges: start -->
<!-- badges: end -->

Provide functionality for finding confidence intervals for parameters of
interests based on different methods. Those functions are for interval
use only.

## Installation

You can install the development version of `CIfinder` like so:

1.  Download the `CIfinder_0.0.0.9003.tgz` from the releases tags

2.  Install the package in command line:

``` r
R CMD INSTALL CIfinder_0.0.0.9003.tgz
```

or install it in R console

``` r
install.packages("CIfinder_0.0.0.9003.tgz", repos = NULL, type="source")
```

## Examples

### Generate the confidence interval for ![PPV](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;PPV "PPV") or ![NPV](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;NPV "NPV")

- ![sensitivity = x_1/n_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;sensitivity%20%3D%20x_1%2Fn_1 "sensitivity = x_1/n_1")

- ![specificity = x_0/n_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;specificity%20%3D%20x_0%2Fn_0 "specificity = x_0/n_0")

- ![ppv=\frac{\rho}{\rho+(1-\rho)\phi\_{ppv}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;ppv%3D%5Cfrac%7B%5Crho%7D%7B%5Crho%2B%281-%5Crho%29%5Cphi_%7Bppv%7D%7D "ppv=\frac{\rho}{\rho+(1-\rho)\phi_{ppv}}")

  where
  ![\rho](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Crho "\rho")
  is the prevalence and
  ![\phi\_{ppv}=\frac{1-specificity}{sensitivity}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cphi_%7Bppv%7D%3D%5Cfrac%7B1-specificity%7D%7Bsensitivity%7D "\phi_{ppv}=\frac{1-specificity}{sensitivity}")

- ![npv=\frac{1-\rho}{(1-\rho)+\rho\phi\_{npv}}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;npv%3D%5Cfrac%7B1-%5Crho%7D%7B%281-%5Crho%29%2B%5Crho%5Cphi_%7Bnpv%7D%7D "npv=\frac{1-\rho}{(1-\rho)+\rho\phi_{npv}}")

  where
  ![\rho](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Crho "\rho")
  is the prevalence and
  ![\phi\_{npv}=\frac{1-sensitivity}{specificity}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cphi_%7Bnpv%7D%3D%5Cfrac%7B1-sensitivity%7D%7Bspecificity%7D "\phi_{npv}=\frac{1-sensitivity}{specificity}")

``` r
library(CIfinder)
ppv_npv_ci(x1 = 53, n1 = 57, x0 = 113, n0 = 113, prevalence = 0.017, method = "gart and nam")
#> $method
#> [1] "gart and nam"
#> 
#> $sensitivity
#> [1] 0.9298246
#> 
#> $specificity
#> [1] 1
#> 
#> $phi_ppv
#> phi_ppv_est   phi_ppv_l   phi_ppv_u phi_ppv_mle 
#>    0.000000    0.000000    0.026435    0.001581 
#> 
#> $ppv
#>   ppv_est     ppv_l     ppv_u   ppv_mle 
#> 1.0000000 0.3954812 1.0000000 0.9162384 
#> 
#> $phi_npv
#> phi_npv_est   phi_npv_l   phi_npv_u phi_npv_mle 
#>  0.07017544  0.02318500  0.15959300  0.07267400 
#> 
#> $npv
#>   npv_est     npv_l     npv_u   npv_mle 
#> 0.9987879 0.9972476 0.9995992 0.9987448
ppv_npv_ci(x1 = 53, n1 = 57, x0 = 113, n0 = 113, prevalence = 0.017, method = "walter")
#> $method
#> [1] "walter"
#> 
#> $sensitivity
#> [1] 0.9298246
#> 
#> $specificity
#> [1] 1
#> 
#> $phi_ppv
#>          phi        lower        upper 
#> 0.0000000000 0.0002976938 0.0753020260 
#> 
#> $ppv
#>     ppv_est ppv_l.upper ppv_u.lower 
#>   1.0000000   0.1867683   0.9830776 
#> 
#> $phi_npv
#>        phi      lower      upper 
#> 0.07017544 0.03223337 0.19001312 
#> 
#> $npv
#>     npv_est npv_l.upper npv_u.lower 
#>   0.9987879   0.9967247   0.9994429
```

### Confidence intervals for single proportion

``` r
single_prop_ci(53, 57)
#> $clopper_ci
#>  cp_lower  cp_upper 
#> 0.8299602 0.9805496 
#> 
#> $wilson
#> wilson_lower wilson_upper 
#>    0.8329983    0.9723736 
#> 
#> $wilson.correct
#> wilson.correct_lower wilson.correct_upper 
#>            0.8216814            0.9772799 
#> 
#> $wald
#> wald_lower wald_upper 
#>  0.8635108  0.9961383 
#> 
#> $agresti
#> agresti_lower agresti_upper 
#>     0.8282120     0.9771599 
#> 
#> $beta_binomial
#> beta_lower beta_upper 
#>  0.8327319  0.9714140
```
