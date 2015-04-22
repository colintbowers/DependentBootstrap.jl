## DependentBootstrap

A module for the Julia language that implements several varieties of the dependent statistical bootstrap as well as the corresponding block-length selection procedures.


## Main Features

This module allows Julia users to estimate the distribution of a test statistic using any of several varieties of the dependent bootstrap. It also includes the corresponding block-length selection procedure for each of these varieties of dependent bootstrap.

The module is implemented entirely in Julia. Each bootstrap procedure and each block-length selection procedure is implemented as its own type. Various functions then operate on these types via multiple dispatch in order to return bootstrap indices, data, statistics, or distribution parameter estimates.

AFAIK, this module is the most comprehensive library of bootstrapping techniques for *univariate* time-series available. It includes moving-block, circular-block, stationary, and tapered-block bootstraps, as well as block-length selection procedures of Politis, White (2004), Patton, Politis, and White (2009), and Paparoditis, Politis (2002).

## What this package does not (yet) include

This package currently only supports bootstrapping of *univariate* time-series. I would be happy to provide whatever support I can to any users who are interested in extending this package to support multivariate time-series. Otherwise, hopefully I will get to this sometime this year (2015).

## How to use DependentBootstrap

#### Installation

This package is not yet registered, so you will not be able to get it using `Pkg.add("DependentBootstrap")`. Further, it also depends on another package of mine that is not yet registered called [KernelStat](https://github.com/colintbowers/KernelStat.jl). This means that if you want to try out DependentBootstrap, you'll need to first clone KernelStat using `Pkg.clone("https://github.com/colintbowers/KernelStat.jl")` at the julia REPL, and then clone DependentBootstrap using `Pkg.clone("https://github.com/colintbowers/DependentBootstrap.jl")` at the julia REPL. This should install both packages in your .julia folder. You should then be able to load either of them into an active julia session with `using(KernelStat)` or `using(DependentBootstrap)` respectively. Note, `using(DependentBootstrap)` will automatically run `using(KernelStat)`.


#### Quick Start

For those who simply want to bootstrap some parameter of the distribution of the statistic of interest, e.g. the variance of the test statistic, or confidence intervals for the test statistic, and don't want to do anything fancier, you can simply input a vector of observable data `x` into the `dBootstrap` function and then use the various keyword arguments to specify the details. A short-list (a full list is provided later in this document) of the most useful key-word arguments and their default values follows:
* numResample::Int=600, the number of bootstrap resamples to construct
* bootstrapMethod::ASCIIString="stationary", the type of dependent bootstrap to use. A full list of options is provided later in this document, but the most popular choices are `"stationary"` for the stationary bootstrap, `"circularBlock"` for the circular-block bootstrap, or `"movingBlock"` for the moving-block bootstrap.
*  blockLength::Int-1, the length of block (or expected length for stationary bootstrap) to use with the dependent bootstrap. A negative value for blockLength (default) means the function will attempt to estimate the optimal block-length using a data-driven method appropriate for the specified bootstrap method.
* statistic::ASCIIString="mean", the test statistic of interest. When using keyword arguments, the user is limited to the following options: `"mean"`, `"median"`, `"variance"`, `"std"` (standard deviation), `"sum"`, `"quantileXX"` (where XX should be set to the percentile of interest, e.g. `"quantile95"` for the 95% quantile). Users can input their own functions in this field, but this option is not compatible with key-word arguments yet (ie you'll have to keep reading this document if you want to do this).
* distributionParam::ASCIIString="variance", the distribution parameter of the test statistic that the user wants to estimate. It is an estimate of the distribution parameter that the `dBootstrap` function returns. When using keyword arguments, the user is limited to the following options: `"mean"`, `"median"`, `"variance"`, `"std"` (standard deviation), `"quantileXX"` (where XX should be set to the percentile of interest, e.g. `"quantile95"` for the 95% quantile), or `"conf"` for a 95% confidence interval (for now, it has to be 95% - I'll make this more flexible soon). The vast majority of users will only be interested in `"variance"`, `"std"`, and `"conf"`. Users can input their own functions in this field, but this option is not compatible with key-word arguments yet (ie you'll have to keep reading this document if you want to do this).

For example, a user wishing to estimate 95% confidence intervals for a sample mean test statistic using a dependent bootstrap would, given observable data `x::Vector{T}`, where `T<:Number`, use:

    (LB, UB) = dBootstrap(x, distributionParam="conf")

where `LB` and `UB` are the returned lower and upper bounds of the confidence interval respectively.

A second example: a user wishing to estimate the standard deviation of a sample variance test statistic using a circular block bootstrap, and where some higher being has informed them that the best block length to use is 5, would use:

    testStatVar = dBootstrap(x, bootstrapMethod="circularBlock", blockLength=5, statistic="variance", distributionParam="std")

where `testStatVar` is the returned esimtate of the variance of the sample variance test statistic.

Note that in both these examples I did not need to specify the keyword arguments where I was happy with the default value.

This concludes the quick-start. Many users will not read beyond this point. However, for those who wish to understand the full flexibility offered by this package, read on!


#### A quick note on block lengths

In the next section, we discuss the types used by this module, including a type for each bootstrap method. The types that represent bootstrap methods invariably (except for one special case) include a field that describes the block length to use with the bootstrapping procedure. 

IMPORTANT: If the user specifies a positive value for this field, then that value will be used as the block length in any calling function. If the user specifies a non-positive value for this field, or does not specify a value for this field (in this case it defaults to a non-positive value), then all calling functions will take that to imply that the user wants the function to estimate an appropriate block length from the supplied data, using either the specified block length estimation procedure, or else defaulting to the procedure most appropriate for the supplied bootstrap method.

In summary, if some higher being tells you what block length to use, then specify it explicitly, otherwise the software will try and automatically detect the block length for you.


#### Types

###### BootstrapParam

The most important type in this package is named `BootstrapParam`. This type stores all the information needed to run a specific bootstrap procedure from start to finish, excluding the actual data that is to be bootstrapped, and every function will accept a `BootstrapParam`, along with the underlying data, as inputs. The fields of `BootstrapParam` are listed now:
* numObsData::Int, the number of observations in your dataset.
* numObsResample::Int, the number of observations that would want *per resample*. The vast majority of users will want numObsData = numObsResample. However, there are legitimate cases where the two should differ.
* numResample::Int, the number of resamples. The module also includes the constant defaultNumResample which is set equal to 600.
* bootstrapMethod::BootstrapMethod, the bootstrap method the user wants to use. The type `BootstrapMethod` is an abstract type that stores as each of its sub-types a specific bootstrap method type. They are discussed in more detail later in this document.
* blockLengthMethod::BlockLengthMethod, the method for automatically choosing an appropriate block length the user wants to use. As with `BootstrapMethod`, `BlockLengthMethod` is an abstract type that stores as each of its sub-types a specific block length method type. They are discussed in more detail later in this document.
* statistic::Union(Function, ASCIIString), the test statistic of interest, ie the user is interested in the distribution of the test statistic stored in this field, which is why they are bootstrapping. If this field is set to a function, it must be a function that converts a single input of type `Vector{T}`, where `T<:Number`, and return a `Float64`. If this field is set to `ASCIIString`, it must be equal to one of the following hard-coded values: `"mean"`, `"median"`, `"variance"`, `"std"` (standard deviation), `"sum"`, `"quantileXX"` (where XX should be set to the percentile of interest, e.g. `"quantile95"` for the 95% quantile)
* distributionParam::Union(Function, ASCIIString), the distribution parameter of the test statistic that the user is interested in. Most user will be interested in either the variance/standard deviation in this field, or else confidence intervals. If this field is set to a function, it must be a function that accepts a single input of type `Vector{Float64}` and it can return whatever the user wants. If this field is set to `ASCIIString`, it must be equal to one of the following hard-coded values: `"mean"`, `"median"`, `"variance"`, `"std"` (standard deviation), `"quantileXX"` (where XX should be set to the percentile of interest, e.g. `"quantile95"` for the 95% quantile), or `"conf"` for a 95% confidence interval (for now, it has to be 95% - I'll make this more flexible soon).

The fields of `BootstrapParam` contain two new types, `BootstrapMethod` and `BlockLengthMethod`. They are discussed now:


###### BootstrapMethod

`BootstrapMethod` is just an abstract type that nests a number of sub-types. Each of these sub-types corresponds to a different bootstrap method, and store as their fields the parameters needed to implement that method. The sub-types and their fields are listed now:
* `BootstrapIID`, the iid bootstrap of Efron (1979) "Bootstrap Methods Another Look at the Jackknife". This type does not need any fields and can therefore be constructed using `IIDBootstrap()`
* `BootstrapStationary`, the stationary bootstrap of Politis, Romano (1994) "The Stationary Bootstrap". This type has the following fields:
  * expectedBlockLength::Float64, the default value constructor of `BootstrapStationary()` intialises this field to -1.0, ie automatic block length selection using the procedure recommended in Politis, White (2004) and Patton, Politis, and White (2009).
* `BootstrapMovingBlock`, the moving block bootstrap propsed in  Kunsch (1989) "The jackknife and the bootstrap for general stationary observations" and (independently) in Liu, Singh (1992) "Moving blocks jackknife and bootstrap capture weak dependence". This type has the following fields:
  * blockLength::Int, the default constructor of `BootstrapMovingBlock()` intialises this field to -1,  ie automatic block length selection using the procedure recommended in Politis, White (2004) and Patton, Politis, and White (2009).
* `BootstrapCircularBlock`, the circular block bootstrap of Politis, Romano (1992) "A circular block resampling procedure for stationary data". This type has the following fields:
  * blockLength::Int, the default constructor of `BootstrapMovingBlock()` intialises this field to -1,  ie automatic block length selection using the procedure recommended in Politis, White (2004) and Patton, Politis, and White (2009).
* `BootstrapNonoverlappingBlock`, the nonoverlapping block bootstrap. A summary can be found in the textbook Lahiri (1999) *Resampling Method for Dependent Data*. This method is not usually used and is included in this module mostly as a curiosity. This type has the following fields:
  * blockLength::Int, the default constructor of `BootstrapMovingBlock()` intialises this field to -1. HOWEVER, I am not aware of an automatic block length selection procedure specifically designed for this method, so users of this method must explcitly set a strictly positive block length or else the software in this module will throw an error. In other words, use the full constructor, eg `BootstrapMovingBlock(5)` to use a block length of 5.
* `BootstrapTaperedBlock`, the tapered block bootstrap of Paparoditis, Politis (2002) "The tapered block bootstrap for general statistics from stationary sequences". Theoretically, this bootstrap procedure should provide superior estimates over all the other procedures in this module. HOWEVER, the data input to this procedure needs to be standardised using the appropriate influence function. I haven't implemented a general procedure for calculating influence functions yet, so currently this procedure is restricted to the case where statistic field of a `BootstrapParam` is equal to the functions `mean` or `sum` or the strings `"mean"` or `"sum"`. In both these cases, the underlying data MUST be de-meaned first. Anything other than what I've described above will throw an error. This type has the following fields:
  * blockLength::Int, as above, the default constructor `BootstrapTaperedBlock()` initialises this field to -1,  ie automatic block length selection using the procedure recommended in Paparoditis, Politis (2002).
  * kernelFunction::KernelFunction, This field requires a kernel function from the package `KernelStat`. Currently, this field must be either the type `KernelPP2002Trap` or `KernelPP2002Smooth`, which corresponds to the trapezoidal and smooth kernel functions discussed in Paparoditis, Politis (2002). The default constructors for these two kernel functions, ie `KernelPP2002Trap()` or `KernelPP2002Smooth()` will initialise the parameters of these kernel functions to the optimal values discussed in Paparoditis, Politis (2002). It is strongly recommended that any users use these default constructors. However, if you want to try something else, see the [docs for KernelStat](https://github.com/colintbowers/KernelStat.jl) for more detail. The default constructor `BootstrapTaperedBlock()` initialises this field to `KernelPP2002Trap`.

This concludes a description of the bootstrap methods. Next we discuss the types that represent different block length selection procedures:


#### BlockLengthMethod

`BlockLengthMethod` is an abstract type that nests a number of sub-types.

