## DependentBootstrap

A module for the Julia language that implements several varieties of the dependent statistical bootstrap as well as the corresponding block-length selection procedures.


## Main Features

This module allows Julia users to estimate the distribution of a test statistic using any of several varieties of the dependent bootstrap. It also includes the corresponding block-length selection procedure for each of these varieties of dependent bootstrap.

The module is implemented entirely in Julia. Each bootstrap procedure and each block-length selection procedure is implemented as its own type. Various functions then operate on these types via multiple dispatch in order to return bootstrap indices, data, statistics, or distribution parameter estimates.

## What this package does not (yet) include

FILL IN HERE


## How to use DependentBootstrap

#### Installation

This package is not yet registered, so you will not be able to get it using `Pkg.add("DependentBootstrap")`. Further, it also depends on another package of mine that is not yet registered called [KernelStat](https://github.com/colintbowers/KernelStat.jl). This means that if you want to try out DependentBootstrap, you'll need to first clone KernelStat using `Pkg.clone("https://github.com/colintbowers/KernelStat.jl")` at the julia REPL, and then clone DependentBootstrap using `Pkg.clone("https://github.com/colintbowers/DependentBootstrap.jl")` at the julia REPL. This should install both packages in your .julia folder. You should then be able to load either of them into an active julia session with `using(KernelStat)` or `using(DependentBootstrap)` respectively. Note, `using(DependentBootstrap)` will automatically run `using(KernelStat)`.


#### Quick Start

For those who simply want to bootstrap some parameter of the distribution of the statistic of interest, e.g. the variance of the test statistic, or confidence intervals for the test statistic, and don't want to do anything fancier, you can simply input a vector of observable data `x` into the `dBootstrap` function and then use the various keyword arguments to specify the details. A short-list (a full list is provided later in this document) of the most useful key-word arguments and their default values follows:
* numResample::Int=600, the number of bootstrap resamples to construct)
* bootstrapMethod::ASCIIString="stationary", the type of dependent bootstrap to use. A full list of options is provided later in this document, but the most popular choices are `"stationary"` for the stationary bootstrap, `"circularBlock"` for the circular-block bootstrap, or `"movingBlock"` for the moving-block bootstrap.
*  blockLength::Int-1, the length of block (or expected length for stationary bootstrap) to use with the dependent bootstrap. A negative value for blockLength (default) means the function will attempt to estimate the optimal block-length using a data-driven method appropriate for the specified bootstrap method.
* statistic::ASCIIString="mean", the test statistic of interest. When using keyword arguments, the user is limited to the following options: `"mean"`, `"median"`, `"variance"`, `"std"` (standard deviation), `"sum"`, `"quantileXX"` (where XX should be set to the percentile of interest, e.g. `"quantile95"` for the 95% quantile). Users can input their own functions in this field, but this option is not compatible with key-word arguments yet (ie you'll have to keep reading this document if you want to do this).
* distributionParam::ASCIIString="variance", the distribution parameter of the test statistic that the user wants to estimate. It is an estimate of the distribution parameter that the `dBootstrap` function returns. When using keyword arguments, the user is limited to the following options: `"mean"`, `"median"`, `"variance"`, `"std"` (standard deviation), `"quantileXX"` (where XX should be set to the percentile of interest, e.g. `"quantile95"` for the 95% quantile), or `"conf"` for a 95% confidence interval (for now, it has to be 95% - I'll make this more flexible soon). The vast majority of users will only be interested in `"variance"`, `"std"`, and `"conf"`. Users can input their own functions in this field, but this option is not compatible with key-word arguments yet (ie you'll have to keep reading this document if you want to do this).

For example, a user wishing to estimate 95% confidence intervals for a sample mean test statistic using a dependent bootstrap would, given observable data `x::Vector{T}`, where `T<:Number`, use:

(LB, UB) = dBootstrap(x, distributionParam="conf")

where `LB` and `UB` are the returned lower and upper bounds of the confidence interval respectively.

A second example: a user wishing to estimate the variance of a sample variance test statistic using a circular block bootstrap, and where some higher being has informed them that the best block length to use is 5, would use:

testStatVar = dBootstrap(x, bootstrapMethod="circularBlock", blockLength=5, statistic="variance")

where `testStatVar` is the returned esimtate of the variance of the sample variance test statistic.

Note that in both these examples I did not need to specify the keyword arguments where I was happy with the default value.

This concludes the quick-start. Many users will not read beyond this point. However, for those who wish to understand the full flexibility offered by this package, read on!


of Politis, Romano (1994) "The Stationary Bootstrap",
of Politis, Romano (1992) "A circular block resampling procedure for stationary data", 
 Kunsch (1989) "The jackknife and the bootstrap for general stationary observations" and (independently) Liu, Singh (1992) "Moving blocks jackknife and bootstrap capture weak dependence".

