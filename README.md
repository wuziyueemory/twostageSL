
<!-- README.md is generated from README.Rmd. Please edit that file -->

# twostageSL: Two-stage Super Learner for predicting healthcare cost with zero-inflation

<!-- badges: start -->

<!-- badges: end -->

> Tools for analyzing healthcare cost data with zero inflation and
> skewed distribution via two stage super learner

**Author:** Ziyue Wu

## Motivation

The study of healthcare expenditures has become an important area in
epidemiological and public health researches as healthcare utilization
and associated costs have increased rapidly in recent years. Estimating
healthcare expenditures is challenging due to heavily skewed
distributions and zero-inflation. Myriad methods have been developed for
analyzing cost data, however, a priori determination of an appropriate
method is difficult. Super-learning, a technique that considers an
ensemble of methods for cost estimation, provides an interesting
alternative for modeling healthcare expenditures as it addresses the
difficulties in choosing the correct model form. The super learner has
shown benefits over a single method in recent studies. Nevertheless, it
is notable that none of the healthcare cost studies have explored the
use of super learner in the case of zero-inflated cost data. Here we
proposed a new method called the two-stage super learner that
implemented the super-learning in a two-part model framework. The
two-stage super learner had natural support for cases of zero healthcare
utilization and extended the standard super learner by providing an
additional layer of the ensemble, thus allowing for greater flexibility
and improving the cost prediction.

## Two-stage Super Learner

To deal with zero-inflation in the positively skewed data, we designed a
method that implemented the super learner under the Two-part model
framework (Two-stage Super Learner). Using
\(E[Y|X] = Pr(Y>0|X)E[Y|Y>0,X]\), the objective of estimating \(E[Y|X]\)
could be broken into two separate pieces: estimating \(Pr(Y>0|X)\) and
\(E[Y|Y>0,X]\). This is done by adding an additional layer of the
ensemble to tackle the problem of a point mass at zero. Rather than
specifying a single library of prediction algorithms as in standard
super learner, we specify two separate libraries of prediction
algorithms, one for each of the two stages. Specifically, we specify the
first stage library by positing a collection of different estimators
predicting the probability of the outcome being positive and specify the
second stage library by positing another collection of estimators
predicting the mean of positive outcomes. Assuming the first stage
library includes \(K1\) methods and the second stage library includes
\(K2\) methods, then the two-stage super learner’s “whole library” would
contain \(K1*K2\) candidate estimators with each one representing a
specific combination of methods from the first stage and second stage.
The two-stage super learner is also an ensembling method built through a
combination of algorithms that provide the best fit to the data.

## Features

`twostageSL` is an R package that provides tools for analyzing health
care cost data useing the super learning technique in the two stage
framework.

  - Automatic optimal predictor ensembling via cross-validation with one
    line of code.
  - Dozens of algorithms: XGBoost, Random Forest, GBM, Lasso, SVM, BART,
    KNN, Decision Trees, Neural Networks, and additional algorithms
    specifically for heavily skewed zero-inflation data: Tweedie, Tobit,
    ZIP, Cox Hazard, Adaptive Hazard, Quantile Regression, Acclerated
    Failture Time.
  - Two options for conducting either a two stage super learner or a
    standard superlearner. Include two seperate libraries for prediction
    algorithms at two stages, accompanied with another library with
    prediction algorithms for standard superlearner.
  - Screen variables (feature selection) based on univariate
    association, Random Forest, Elastic Net, et al. or custom screening
    algorithms.
  - Includes framework to quickly add custom algorithms to the ensemble.

## Installation

twostageSL can be installed from Github using the following command:

``` r
# install.packages("devtools")

library(devtools)
devtools::install_github("wuziyueemory/twostageSL")

library(twostageSL)
```

## Example

To show an example on how to use this package, we simulated an example
with training data containing 1000 observations, 5 predictors and
testing data containing 100 observations, 5 predictions, all following a
normal distribution. The outcome is continuous and non-negative with
relatively high proportion of zeros.

Generate training and testing set.

``` r
library(twostageSL)
set.seed(121)
## training set
n <- 1000
p <- 5
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X) <- paste("X", 1:p, sep="")
X <- data.frame(X)
Y <- rep(NA,n)
## probability of outcome being zero
prob <- plogis(1 + X[,1] + X[,2] + X[,1]*X[,2])
g <- rbinom(n,1,prob)
## assign zero outcome
ind <- g==0
Y[ind] <- 0
## assign non-zero outcome
ind <- g==1
Y[ind] <- 10 + X[ind, 1] + sqrt(abs(X[ind, 2] * X[ind, 3])) + X[ind, 2] - X[ind, 3] + rnorm(sum(ind))

## test set
m <- 100
newX <- matrix(rnorm(m*p), nrow = m, ncol = p)
colnames(newX) <- paste("X", 1:p, sep="")
newX <- data.frame(newX)
newY <- rep(NA,m)
## probability of outcome being zero
newprob <- plogis(1 + newX[,1] + newX[,2] + newX[,1]*newX[,2])
newg <- rbinom(n,1,newprob)
## assign zero outcome
newind <- newg==0
newY[newind] <- 0
## assign non-zero outcome
newind <- g==1
newY[newind] <- 10 + newX[newind, 1] + sqrt(abs(newX[newind, 2] * newX[newind, 3])) + newX[newind, 2] - X[newind, 3] + rnorm(sum(newind))
```

The training data looks like

    #>            X1         X2          X3         X4         X5
    #> 1 -0.25536047  0.5210512  0.07889016  1.5794606 -0.2476091
    #> 2  0.10837472  0.8525208  0.26945885 -0.4203506  0.5246956
    #> 3  0.12778056 -1.9302377 -0.33684853 -0.9691443 -0.6481140
    #> 4 -0.08107345 -1.4007554  0.23814700 -0.2929321 -1.9521166
    #> 5 -0.71368372  0.4402426 -1.15347200 -0.5136347 -1.3213276
    #> 6  1.61533208  0.2069617 -0.04666425  0.2745282 -2.6145711

The proportion of zeros in outcome Y.

    #> [1] 0.368

The distribution of outcome Y, which is zero-inflated and heavily
skewed. <img src="man/figures/README-example 4-1.png" width="100%" />

`twostageSL` is the core function to fit the two stage super learner. At
a minimum for implementation, you need to specify the predictor matrix
X, outcome variable Y, library of prediction algorithms at two stages,
boolean for two stage superlearner or standard superlearner,
distribution at two stages and number of folds for cross-validation. The
package incldues all the prediction algorithms from package
`SuperLearner` and also includes additional prediction algorithms from
package `slcost` that specifically used at stage 2 to deal with heavily
skewed data.

Generate the library and run the two stage super learner

``` r
## generate the Library
twostage.library <- list(stage1=c("SL.mean","SL.glm","SL.earth"),
                       stage2=c("SL.mean","SL.glm","SL.earth","SL.coxph"))

onestage.library <- c("SL.mean","SL.glm","SL.earth")


## run the twostage super learner
two <- twostageSL(Y=Y,
                   X=X,
                   newX = newX,
                   library.2stage <- twostage.library,
                   library.1stage <- onestage.library,
                   twostage = TRUE,
                   cvControl = list(V = 5))
two
#> 
#> Call:  
#> twostageSL(Y = Y, X = X, newX = newX, library.2stage = library.2stage <- twostage.library,  
#>     library.1stage = library.1stage <- onestage.library, twostage = TRUE,  
#>     cvControl = list(V = 5)) 
#> 
#> 
#>                                         Risk      Coef
#> S1: SL.mean_All + S2: SL.mean_All   30.99158 0.0000000
#> S1: SL.mean_All + S2: SL.glm_All    26.02499 0.0000000
#> S1: SL.mean_All + S2: SL.earth_All  26.34461 0.0000000
#> S1: SL.mean_All + S2: SL.coxph_All  26.25908 0.0000000
#> S1: SL.glm_All + S2: SL.mean_All    22.91554 0.0000000
#> S1: SL.glm_All + S2: SL.glm_All     20.80390 0.0000000
#> S1: SL.glm_All + S2: SL.earth_All   20.97040 0.0000000
#> S1: SL.glm_All + S2: SL.coxph_All   20.80807 0.2111024
#> S1: SL.earth_All + S2: SL.mean_All  21.59218 0.0000000
#> S1: SL.earth_All + S2: SL.glm_All   19.64819 0.4474539
#> S1: SL.earth_All + S2: SL.earth_All 19.64437 0.1748323
#> S1: SL.earth_All + S2: SL.coxph_All 19.65313 0.0000000
#> Single: SL.mean_All                 30.99155 0.0000000
#> Single: SL.glm_All                  21.34758 0.0000000
#> Single: SL.earth_All                20.13396 0.1666114
```

When setting `twostage` to FALSE, we get the result from standard super
learner

``` r
## run the standard super learner
one <- twostageSL(Y=Y,
                  X=X,
                  newX = newX,
                  library.2stage <- twostage.library,
                  library.1stage <- onestage.library,
                  twostage = FALSE,
                  cvControl = list(V = 5))
one
#> 
#> Call:  
#> SuperLearner(Y = Y, X = X, family = family.single, SL.library = library.1stage,  
#>     method = method.CC_LS.scale, verbose = verbose, control = list(saveCVFitLibrary = T),  
#>     cvControl = cvControl) 
#> 
#> 
#>                  Risk       Coef
#> SL.mean_All  30.99155 0.01539633
#> SL.glm_All   21.34758 0.27146007
#> SL.earth_All 20.13396 0.71314359
```

We could also use `CV.twostageSL` to get v-fold cross validated risk
estimate for two stage super learner, standard super learner, and all
prediction algorithms in the library. In this case we set numbert of
folds to be 2.

``` r
two_cv <- CV.twostageSL(
  Y = Y, X = X,
  library.2stage = twostage.library,
  library.1stage = onestage.library,
  cvControl = list(V = 2),
  innerCvControl = list(list(V = 5),
                        list(V = 5))
)
two_cv
#> 
#> Call:  
#> CV.twostageSL(Y = Y, X = X, library.2stage = twostage.library, library.1stage = onestage.library,  
#>     cvControl = list(V = 2), innerCvControl = list(list(V = 5), list(V = 5))) 
#> 
#> 
#> 
#> Cross-validated predictions from the two-stage Super Learner:  SL.predict 
#> 
#> Cross-validated predictions from the discrete super learner (cross-validation selector):  discreteSL.predict 
#> 
#> Which library algorithm was the discrete super learner:  whichDiscreteSL 
#> 
#> Cross-validated prediction for all algorithms in the library and standard super learner:  library.predict
```

The summary and plots are shown below, with the two stage super learner
outperforming all other prediction algorithms.

    #> 
    #> Call:  
    #> CV.twostageSL(Y = Y, X = X, library.2stage = twostage.library, library.1stage = onestage.library,  
    #>     cvControl = list(V = 2), innerCvControl = list(list(V = 5), list(V = 5))) 
    #> 
    #> Risk is based on: Mean Squared Error
    #> 
    #> All risk estimates are based on V =  2 
    #> 
    #>                            Algorithm    Ave      se    Min    Max
    #>              Two-stage Super Learner 19.514 0.72115 19.013 20.015
    #>                          Discrete SL 19.895 0.84434 18.929 20.861
    #>    S1: SL.mean_All + S2: SL.mean_All 30.991 0.64242 30.953 31.030
    #>     S1: SL.mean_All + S2: SL.glm_All 26.036 0.49851 25.657 26.415
    #>   S1: SL.mean_All + S2: SL.earth_All 26.168 0.50502 25.790 26.545
    #>   S1: SL.mean_All + S2: SL.coxph_All 26.206 0.48504 25.714 26.697
    #>     S1: SL.glm_All + S2: SL.mean_All 23.099 0.67940 22.807 23.391
    #>      S1: SL.glm_All + S2: SL.glm_All 20.889 0.71529 20.883 20.894
    #>    S1: SL.glm_All + S2: SL.earth_All 20.952 0.71781 20.923 20.981
    #>    S1: SL.glm_All + S2: SL.coxph_All 20.919 0.70313 20.908 20.929
    #>   S1: SL.earth_All + S2: SL.mean_All 21.466 0.79182 19.988 22.945
    #>    S1: SL.earth_All + S2: SL.glm_All 19.921 0.84741 18.837 21.005
    #>  S1: SL.earth_All + S2: SL.earth_All 19.895 0.84434 18.929 20.861
    #>  S1: SL.earth_All + S2: SL.coxph_All 19.888 0.83303 18.854 20.923
    #>                  Single: SL.mean_All 30.991 0.64242 30.953 31.030
    #>                   Single: SL.glm_All 21.451 0.71995 20.942 21.961
    #>                 Single: SL.earth_All 21.004 0.78972 19.124 22.884
    #>                Standard SuperLearner 20.184 0.69532 19.268 21.101

<img src="man/figures/README-example 8-1.png" width="100%" />

We could use `predict` to obtain predictions on a new dataset from the
twostageSL fit. The `onlySL` argument controls whether to compute
predictions for all algorithms or for algorithms with non-zero
coefficients in the two stage super learner fit.

``` r
newdata <- newX[c(1:50),]
prediction <- predict(two,newdata = newX, onlySL = TRUE)
#> Predictions from the two-stage Super Learner:  pred 
#> 
#> Prediction for all algorithms in the library:  library.predict
```

This returns a list with predictions from two stage super learner
(`pred`) and predictions for each algorithm in library
(`library.predict`).
