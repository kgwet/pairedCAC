Testing the difference of 2 agreement coefficients for statistical
significance
================
Kilem L. Gwet, Ph.D.
2021-09-16

<!-- README.md is generated from README.Rmd. Please edit that file -->

# pairedCAC

<!-- badges: start -->
<!-- badges: end -->

The goal of **pairedCAC** is to allow you to test the difference between
2 chance-corrected agreement coefficients for statistical significance.

``` r
library(pairedCAC)
```

## Installation

You can install the most current version of **pairedCAC** from the
GitHub repositure as follows:

``` r
devtools::install_github("kgwet/pairedCAC")
```

## Abstract

The **pairedCAC** is an R package that provides a series of functions
for testing the difference between 2 **correlated** chance-corrected
agreement coefficients for statistical significance. This package
closely follows the general framework discussed by Gwet (2016) and
expanded more recently by Gwet (2021).

Typically, 2 agreement coefficients are correlated when they are based
on 2 overlapping samples of subjects or 2 overlapping rosters of raters.
For uncorrelated coefficients, the testing is straightforward and
discussed in Gwet(2021, section 9.3).

## Example

Suppose that Fleiss’ generalized coefficient is calculated on 2
occasions using the 2 datasets of ratings *ratings1* and *ratings2*
included in this package. How do you determine whether the difference
between these 2 Fleiss’ kappa coefficients is statistically significant?
Proceed as follows:

``` r
fleiss <- ttest.fleiss(ratings1,ratings2)
fleiss$test
#>   fleiss.coeff1 fleiss.coeff2 coeff.diff   std.err   t.stat    p.value n.obs
#> 1     0.3856655     0.7248495   0.339184 0.1574322 2.154477 0.04910569    15
#>   n.raters1 n.raters2
#> 1         4         4
```

The function `ttest.fleiss()` returns a 2-element list, the first of
which is a data frame named `test` that you can display with the
expression `fleiss$test`. This data frame contains the test statistics
shown above.

Fleiss’ generalized kappa coefficients associated with the datasets 1
and 2 are respectively given by **fleiss.coeff1 =0.386** and
**fleiss.coeff2 = 0.725**. The difference between the second and first
coefficient that is to be tested for statistical significance is given
by **coeff.diff = 0.339** and its standard error is **std.err = 0.157**.

Now, a key element for testing the statistical significance of the
difference is the *Test Statistic* *t.stat = 2.154* that you must
compare to the 97.5-th of the Standard Normal distribution, which is
1.96. Since the test statistic exceeds this threshold, you may conclude
that the difference is statistically significant.

The `ttest.fleiss()` function also outputs the p-value associated with
the difference between the 2 Fleiss’ kappa coefficients. In our example,
it is given by **p-value = 0.0491**, which is smaller than the standard
threshold of 0.05. This is an indication that the difference is
statistically significant.

Here are the 2 data frames `ratings1` and `ratings2` and the 6 functions
available for testing the difference between agreement coefficients for
statistical significance.

``` r
data.frame(ratings1,ratings2)
#>    RaterA RaterB RaterC RaterD RaterA.1 RaterB.1 RaterC.1 RaterD.1
#> 1       1      1      1      1        1        1        1        1
#> 2       1      1      2      3        3        3        3        2
#> 3       1      1      2      3        3        3        3        3
#> 4       2      2      2      3        1        1        2        1
#> 5       1      2      1      1        1        1        1        1
#> 6       2      1      1      1        1        1        1        1
#> 7       2      2      2      2        2        2        2        2
#> 8       2      3      2      2        2        2        1        1
#> 9       3      3      3      3        3        3        3        3
#> 10      1      1      1      1        1        1        1        3
#> 11      1      1      2      3        2        2        2        2
#> 12      3      3      3      3        3        3        3        3
#> 13      3      2      3      3        3        3        3        3
#> 14      1      2      1      1        1        2        1        1
#> 15      2      1      1      1        3        3        3        3
fleiss <- ttest.fleiss(ratings1,ratings2[,1:3])
ac2 <- ttest.ac2(ratings1,ratings2)
conger <- ttest.conger(ratings1,ratings2)
alpha <- ttest.alpha(ratings1,ratings2)
bp <- ttest.bp(ratings1,ratings2)
pa <- ttest.pa(ratings1,ratings2)
fleiss
#> $test
#>   fleiss.coeff1 fleiss.coeff2 coeff.diff   std.err  t.stat    p.value n.obs
#> 1     0.3856655     0.7942073  0.4085418 0.1612858 2.53303 0.02389035    15
#>   n.raters1 n.raters2
#> 1         4         3
#> 
#> $weights
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
ac2
#> $test
#>   ac2.coeff1 ac2.coeff2 coeff.diff   std.err   t.stat    p.value n.obs
#> 1  0.4069193   0.737382  0.3304628 0.1452689 2.274835 0.03917311    15
#>   n.raters1 n.raters2
#> 1         4         4
#> 
#> $weights
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
conger
#> $test
#>   conger.coeff1 conger.coeff2 coeff.diff   std.err   t.stat    p.value n.obs
#> 1     0.3932584     0.7250859  0.3318275 0.1527202 2.172781 0.04745695    15
#>   n.raters1 n.raters2
#> 1         4         4
#> 
#> $weights
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
alpha
#> $test
#>   alpha.coeff1 alpha.coeff2 coeff.diff   std.err   t.stat    p.value n.obs
#> 1    0.3959044    0.7294354  0.3335309 0.1548083 2.154477 0.04910569    15
#>   n.raters1 n.raters2
#> 1         4         4
#> 
#> $weights
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
bp
#> $test
#>   bp.coeff1 bp.coeff2 coeff.diff   std.err   t.stat    p.value n.obs n.raters1
#> 1       0.4 0.7333333  0.3333333 0.1477342 2.256304 0.04056859    15         4
#>   n.raters2
#> 1         4
#> 
#> $weights
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
pa
#> $test
#>   pa.coeff1 pa.coeff2 coeff.diff    std.err   t.stat    p.value n.obs n.raters1
#> 1       0.6 0.8222222  0.2222222 0.09848947 2.256304 0.04056859    15         4
#>   n.raters2
#> 1         4
#> 
#> $weights
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
```

# References:

1.  Gwet, K.L. (2016), [Testing the Difference of Correlated Agreement
    Coefficients for Statistical
    Significance](https://agreestat.com/papers/correlated_agreement_coefficients_educational_and_psychological_measurement_2016.pdf),
    *Educational and Psychological Measurement*, Vol. 76(4) 609–637.

2.  Gwet, K.L. (2021,
    [ISBN:978-1792354632](https://www.amazon.com/Handbook-Inter-Rater-Reliability-Categorical-Definitive/dp/1792354630/)).
    “*Handbook of Inter-Rater Reliability, Volume 1: Analysis of
    Categorical Ratings*,” 5th Edition. AgreeStat Analytics
