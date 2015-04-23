## DependentBootstrap

A module for the Julia language that implements several varieties of the dependent statistical bootstrap as well as the corresponding block-length selection procedures.


## Main Features

This module allows Julia users to estimate the distribution of a test statistic using any of several varieties of the dependent bootstrap. It also includes the corresponding block-length selection procedure for each of these varieties of dependent bootstrap.

The module is implemented entirely in Julia. Each bootstrap procedure and each block-length selection procedure is implemented as its own type. Various functions then operate on these types via multiple dispatch in order to return bootstrap indices, data, statistics, or distribution parameter estimates.

AFAIK, this module is the most comprehensive library of bootstrapping techniques for *univariate* time-series available. It includes moving-block, circular-block, stationary, and tapered-block bootstraps, as well as block-length selection procedures of Politis, White (2004), Patton, Politis, and White (2009), and Paparoditis, Politis (2002).

## What this package does not (yet) include

This package currently only supports bootstrapping of *univariate* time-series. I would be happy to provide whatever support I can to any users who are interested in extending this package to support multivariate time-series or bootstrapping in a linear regression framework. Otherwise, hopefully I will get to this sometime this year (2015).

## A Quick Introduction to DependentBootstrap

#### Installation

This package is not yet registered, so you will not be able to get it using `Pkg.add("DependentBootstrap")`. Further, it also depends on another package of mine that is not yet registered called [KernelStat](https://github.com/colintbowers/KernelStat.jl). This means that if you want to try out DependentBootstrap, you'll need to first clone KernelStat using `Pkg.clone("https://github.com/colintbowers/KernelStat.jl")` at the julia REPL, and then clone DependentBootstrap using `Pkg.clone("https://github.com/colintbowers/DependentBootstrap.jl")` at the julia REPL. This should install both packages in your .julia folder. You should then be able to load either of them into an active julia session with `using(KernelStat)` or `using(DependentBootstrap)` respectively. Note, `using(DependentBootstrap)` will automatically run `using(KernelStat)`.


#### Quick Start

For those who simply want to bootstrap some parameter of the distribution of the statistic of interest, e.g. the variance of the test statistic, or confidence intervals for the test statistic, and don't want to do anything fancier, you can simply input a vector of observable data `x` into the `dBootstrap` function and then use the various keyword arguments to specify the details. A short-list (a full list is provided later in this document) of the most useful key-word arguments and their default values follows:
* numResample::Int=600, the number of bootstrap resamples to construct
* bootstrapMethod::ASCIIString="stationary", the type of dependent bootstrap to use. A full list of options is provided later in this document, but the most popular choices are `"stationary"` for the stationary bootstrap, `"circularBlock"` for the circular-block bootstrap, or `"movingBlock"` for the moving-block bootstrap.
*  blockLength::Int-1, the length of block (or expected length for stationary bootstrap) to use with the dependent bootstrap. A non-positive value for blockLength (default) means the function will attempt to estimate the optimal block-length using a data-driven method appropriate for the specified bootstrap method.
* statistic::ASCIIString="mean", the test statistic of interest. When using keyword arguments, the user is limited to the following options: `"mean"`, `"median"`, `"variance"`, `"std"` (standard deviation), `"sum"`, `"quantileXX"` (where XX should be set to the percentile of interest, e.g. `"quantile95"` for the 95% quantile). Users can input their own functions in this field, but this option is not compatible with key-word arguments yet (ie you'll have to keep reading this document if you want to do this).
* distributionParam::ASCIIString="variance", the distribution parameter of the test statistic that the user wants to estimate. It is an estimate of the distribution parameter that the `dBootstrap` function returns. When using keyword arguments, the user is limited to the following options: `"mean"`, `"median"`, `"variance"`, `"std"` (standard deviation), `"quantileXX"` (where XX should be set to the percentile of interest, e.g. `"quantile95"` for the 95% quantile), or `"conf"` for a 95% confidence interval (for now, it has to be 95% - I'll make this more flexible soon). The vast majority of users will only be interested in `"variance"`, `"std"`, and `"conf"`. Users can input their own functions in this field, but this option is not compatible with key-word arguments yet (ie you'll have to keep reading this document if you want to do this).

For example, a user wishing to estimate 95% confidence intervals for a sample mean test statistic using a dependent bootstrap would, given observable data `x::Vector{T}`, where `T<:Number`, use:

    (LB, UB) = dBootstrap(x, distributionParam="conf")

where `LB` and `UB` are the returned lower and upper bounds of the confidence interval respectively.

A second example: a user wishing to estimate the standard deviation of a sample variance test statistic using a circular block bootstrap, and where some higher being has informed them that the best block length to use is 5, would use:

    distParam = dBootstrap(x, bootstrapMethod="circularBlock", blockLength=5, statistic="variance", distributionParam="std")

where `distParam` is the returned esimtate of the standard deviation of the sample variance test statistic.

Note that in both these examples I did not need to specify the keyword arguments if I was happy with the default value.

Some users may wish instead to plot an empirical distribution of their bootstrapped test statistics, rather than obtain a specific distribution parameter. I have not yet implemented this capability, however, these users can use the function `dBootstrapStatistic` in order to return the bootstrapped test statistics (as type `Vector{Float64}`), eg

    bootTestStatVec = dBootstrapStatistic(x, ...)
    
and then build the histogram themselves. Note, the valid keyword arguments to `dBootstrapStatistic` are identical to those described above, with the exception of distributionParam, which is obviously no longer relevant since we're not computing any distribution parameters. In fact, the `dBoostrap` function is actually just a wrapper on `dBootstrapStatistic` that performs the additional step of estimating a distribution parameter from the bootstrapped test statistics. In a similar fashion, `dBootstrapStatistic` is a wrapper on `dBootstrapData`, and `dBootstrapData` is a wrapper on `dBootstrapIndex`. Interested users can explore some of these other functions, or else keep reading to find out more.

This concludes the quick-start. Many users will not read beyond this point. However, for those who wish to understand the full flexibility offered by this package, read on!


#### A quick note on block lengths

In the next section, we discuss the types used by this module, including a type for each bootstrap method. The types that represent bootstrap methods invariably (except for one special case) include a field that describes the block length to use with the bootstrapping procedure. 

IMPORTANT: If the user specifies a positive value for this field, then that value will be used as the block length in any calling function. This trumps any other consideration, eg a supplied block length method etc. If the user specifies a non-positive value for this field, or does not specify a value for this field (in this case it defaults to a non-positive value), then all calling functions will take that to imply that the user wants the function to estimate an appropriate block length from the supplied data, using either the specified block length estimation procedure, or else defaulting to the procedure most appropriate for the supplied bootstrap method.

In summary, if some higher being tells you what block length to use, then specify it explicitly, otherwise the software will try and automatically detect the block length for you.


## A complete description of DependentBootstrap

#### Type Definitions

###### BootstrapParam

The most important type in this package is `BootstrapParam`. This type stores all the information needed to run a specific bootstrap procedure from start to finish, excluding the actual data that is to be bootstrapped. Almost every function will accept a `BootstrapParam`, along with the underlying data (in either order) as inputs. The fields of `BootstrapParam` are listed now, along with their default values (supplied after the equals sign):
* numObsData::Int, the number of observations in your dataset. One way or another, this field must be supplied by the user, so it has no default value.
* numObsResample::Int=numObsData, the number of observations that the user wants *per resample*. The vast majority of users will want numObsResample=numObsData, hence this is the default value for this field. However, there are legitimate cases where the two should differ.
* numResample::Int=defaultNumResample, the number of resamples. The module includes a constant `defaultNumResample` which is set equal to 600.
* bootstrapMethod::BootstrapMethod=BootstrapStationary(defaultBlockLength)). This field describes the bootstrap method the user wants to use. The type `BootstrapMethod` is an abstract type that stores as each of its sub-types a specific bootstrap method type. They are discussed in more detail later in this section. For now, it suffices to say that the default for this field is the stationary bootstrap using the constant (for the module) `defaultBlockLength` which is set to -1. A non-positive block length is interpreted by the functions in this module as implying that you want to use an appropriate data-driven method to estimate the block length.
* blockLengthMethod::BlockLengthMethod=BlockLengthPPW2009(BandwidthP2003(), "stationary"). This field describes the method for automatically choosing an appropriate block length. As with `BootstrapMethod`, `BlockLengthMethod` is an abstract type that stores as each of its sub-types a specific block length method type. They are discussed in more detail later in this document. For now, it suffices to say that the default block length detection method is the on described in Politis, White (2004) "Automatic Block Length Selection For The Dependent Bootstrap", incorporating the correction in Patton, Politis, White (2009) "Correction To Automatic Block Length Selection For The Dependent Bootstrap". Specifically, it is the method that applies the stationary bootstrap, and uses the bandwidth estimator of Politis (2003) "Adaptive Bandwidth Choice".
* statistic::Union(Function, ASCIIString). This field stores the test statistic of interest, ie the test statistic that the user is interested bootstrapping, in order to estimate a parameter of the distribution of the test statistic. If this field is set to a function, then the function must accept a single input of type `Vector{T}`, where `T<:Number`, and return a `Float64` (for example, `mean`). If this field is set to `ASCIIString`, it must be equal to one of the following hard-coded values: `"mean"`, `"median"`, `"variance"`, `"std"` (standard deviation), `"sum"`, `"quantileXX"` (where XX should be set to the percentile of interest, e.g. `"quantile95"` for the 95% quantile)
* distributionParam::Union(Function, ASCIIString). This field stores the distribution parameter of the test statistic that the user is interested in. Most user will be interested in either the variance/standard deviation or else confidence intervals. If this field is set to a function, it must be a function that accepts a single input of type `Vector{Float64}` and it can return whatever the user wants. If this field is set to `ASCIIString`, it must be equal to one of the following hard-coded values: `"mean"`, `"median"`, `"variance"`, `"std"` (standard deviation), `"quantileXX"` (where XX should be set to the percentile of interest, e.g. `"quantile95"` for the 95% quantile), or `"conf"` for a 95% confidence interval (for now, it has to be 95% - I'll make this more flexible soon).

The fields of `BootstrapParam` contain two new abstract types, `BootstrapMethod` and `BlockLengthMethod`. They are discussed now:


###### BootstrapMethod

`BootstrapMethod` is an abstract type that nests a number of sub-types. Each of these sub-types corresponds to a different bootstrap method, and store as their fields the parameters needed to implement that method. The sub-types and their fields are listed now:
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


###### BlockLengthMethod

`BlockLengthMethod` is an abstract type that nests a number of sub-types. Each of these sub-types corresponds to a different block length estimation method, and store as their fields the parameters needed to implement that method. The sub-types and their fields are listed now:
* `BlockLengthPPW2009`, the block length selection procedure of Politis, White (2004) "Automatic Block Length Selection For The Dependent Bootstrap", incorporating the correction in Patton, Politis, White (2009) "Correction To Automatic Block Length Selection For The Dependent Bootstrap". This method is capable of estimating an appropriate block length for use with the moving block, circular block, and stationary bootstraps. It contains the fields:
  * bandwidthMethod::BandwidthMethod. This field specifies the method for estimating the bandwidth as part of the block-length selection procedure. The default constructor sets this field to `BandwidthP2003` (see the docs for the KernelStat package for more detail on this and other bandwidth selection procedures).
  * bootstrapMethod::ASCIIString. This field specifies the bootstrap method that the block length procedure is being used for. It can be set to "movingBlock", "circularBlock", or "stationary". There are small differences in the PPW2009 block-length selection procedure depending on which bootstrap method it is being used with (see Politis, White (2004) for more detail).
* `BlockLengthPP2002`, the block length selection procedure of Paparoditis, Politis (2002) "The tapered block bootstrap for general statistics from stationary sequences". Note, this procedure is designed for use with the tapered block bootstrap. `BlockLengthPP2002` contains the fields:
  * bandwidthMethod::BandwidthMethod. See description of bandwidthMethod field for `BlockLengthPPW2009` above as usage is identical.
  * kernelFunction::KernelFunction. This field specifies the kernel function being used with the tapered block bootstrap. There are slight variations in how the optimal block length is selected depending on the kernel function in use. Currently, the only valid values for this field are `KernelPP2002Trap` and `KernelPP2002Smooth` (the trapezoidal and smooth kernel function described in Paparoditis, Politis (2002)). See the KernelStat package for more detail on these kernel functions.

Note that there is no block length selection procedure designed for use with the nonoverlapping block bootstrap. I am not aware of one in the literature, so users who wish to use this bootstrap procedure are expected to either supply their own block length, or else use one of the procedures listed above.


#### Type Constructors

Sub-types of `BlockLengthMethod` and `BootstrapMethod` typically only have one or two fields, and so generally have two constructors, the classic inner constructor where all fields are specified, and a default constructor with no inputs, eg `BootstrapStationary()`, that initialises the fields to sensible default values. An exception to this rule is the bandwidthMethod field of the subtypes of `BlockLengthMethod` which can typically be constructed using an `ASCIIString` input, eg

    bLM = BlockLengthPPW2009("P2003", "stationary")

to construct a type representing the block length selection procedure of Politis, White (2004) for the stationary bootstrap using the bandwidth selection procedure of Politis (2003), or

    blM  = BlockLengthPP2002("bartlett", "PP2002Trap")

to construct a type representing the block length selection procedure of Paparoditis, Politis (2002), assuming the use of the trapezoidal kernel function with the tapered block bootstrap, and using the Bartlett method for estimating bandwidth.

In the above examples, the string is converted into the appropriate type using the default constructor for that type. To find the appropriate string to associate with any given type, use `string(SomeType)`, eg `string(BandwidthP2003())` will evaluate to `"P2003"`, and `string(KernelPP2002Trap())` will evaluate to `"PP2002Trap"`.

The type `BootstrapParam` is central to the entire module, since, as mentioned, nearly every function can accept a vector of data and a `BootstrapParam` as inputs (in either order). We list all constructors for this type now. A list of non-keyword constructors for `BootstrapParam` follows:

* `BootstrapParam(numObsData::Int)`. Number of observations supplied, and use default values (listed above) for all other fields.
* `BootstrapParam(numObsData::Int, statistic::Union(Function, ASCIIString), distributionParam::Union(Function, ASCIIString))`. Number of observations, test statistic of interest, and distribution parameter of test statistic, all supplied, and use default values for all other fields.
* `BootstrapParam(numObsData::Int, numResample::Int)`. Number of observations and number of resamples supplied, and use default values for all other fields.
* `BootstrapParam(numObsData::Int, numResample::Int, statistic::Union(Function, ASCIIString), distributionParam::Union(Function, ASCIIString))`. Number of observations, number of resamples, test statistic, and distribution parameter all supplied, and use default values for all other fields.
* `BootstrapParam{T<:Number}(x::Vector{T})`. Number of observations equal to `length(x)`, and use default values for all other fields, except block-length which the constructor attempts to estimate using the data in `x` and the default block-length method.
* `BootstrapParam{T<:Number}(x::Vector{T}, statistic::Union(Function, ASCIIString), distributionParam::Union(Function, ASCIIString))`. As above, but test statistic of interest, and distribution parameter of test statistic of interest, also supplied.
* `BootstrapParam{T<:Number}(x::Vector{T}, numResample::Int)`. Number of observations equal to `length(x)`, number of resamples supplied, and use default values for all other fields, except block-length which the constructor attempts to estimate using the data in `x` and the default block-length method.
* `BootstrapParam{T<:Number}(x::Vector{T}, numResample::Int, statistic::Union(Function, ASCIIString), distributionParam::Union(Function, ASCIIString))`. As above, except the test statistic and distribtion parameter are also supplied.

This type also allows key-word constructors. These are important, since all of the functions in this module that admit keyword arguments do this by inputting those keywords into the `BootstrapParam` keyword constructor, and then feeding that `BootstrapParam` into the appropriate function. The two keyword constructors are provided now:

* `BootstrapParam(numObsData::Int; numObsResample::Int=numObsData, numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003", statistic::ASCIIString="mean", distributionParam::ASCIIString="variance")`. The number of observations is supplied as the first input, and then the user can specify all other inputs using keyword arguments. Currently, this constructor limits the keyword arguments to numbers or strings. I plan on making this more flexible in the future.
* `BootstrapParam{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003", statistic::ASCIIString="mean", distributionParam::ASCIIString="variance")`. The observed data is supplied as the first input, and then the user can specify all other inputs using keyword arguments. If the supplied keyword argument blockLength is non-positive (or is not specified), then the constructor will attempt to estimate the block-length using the specified (or default) block length method and the input observed data.


#### Functions

DependentBootstrap exports several functions, most of which nest each other, adding one additional piece of functionality as we move up the chain. This can be conceptualised as follows:

`dBootstrapIndex` -> `dBootstrapData` -> `dBootstrapStatistic` -> `dBootstrap`

where additional functionality is provided as we move across the page. There is also the stand-alone function `dBootstrapBlockLength` which estimates the block length, and `dBootstrapVar`, `dBootstrapStd`, `dBootstrapConf`, all of which wrap `dBootstrap` but with hard-coded values for the distributionParam field, ie distributionParam field is set to variance, standard deviation, or confidence intervals respectively. There are a few other stand-alone functions that we will cover at the very end.

###### `dBootstrapIndex`

This function is arguably the most important function in the module. It utilises multiple dispatch over the inputs (::BootstrapMethod, numObsData::Int, numObsResample::Int, numResample::Int) in order to return a set of bootstrap indices. The bootstrap indices are used to index into the original vector of data in order to create all of the bootstrap re-samples.

This function has many methods. Interested users should refer to the source code. However, the vast majority of users who are interested in using this function on its own can gain all possible flexibility using just one method:

    dBootstrapIndex(bp::BootstrapParam)

which returns a `Matrix{Int}` of bootstrap indices that can be used to index into the original data vector to build bootstrap resamples (one resampled column of data per column of the output of `dBootstrapIndex`).

The fields of a `BootstrapParam` are sufficient to tell the function exactly which type of bootstrap indices to generate. Two other methods that might be of interest are based on keywords:
* `dBootstrapIndex{T<:Number}(numObsData::Int, blockLength::T; numObsResample::Int=numObsData, numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary")`
* `dBootstrapIndex{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003")`


###### `dBootstrapData`

This function wraps `dBootstrapIndex`, and in the vast majority of cases, the only functionality it adds is to use the output of `dBootstrapIndex` to index into the observable data vector. The one exception to this is for the tapered block bootstrap, where this function will also take care of applying the appropriate weighting to the returned bootstrapped data.

The method that provides maximum flexibility to this function is:

    dBootstrapData{T<:Number}(x::Vector{T}, bp::BootstrapParam)
    
which returns a `Matrix{T}` of bootstrapped data, where each column corresponds to resample of the underlying data vector `x`. Some users may also wish to use keyword functionality via the following method:

    dBootstrapData{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003")
    

###### `dBootstrapStatistic`

This function wraps `dBootstrapData`. The only functionality added by `dBootstrapStatistic` is to loop over the columns of the output of `dBootstrapData`, and to compute the test statistic of interest from each column.  Maximum flexibility can be obtained using the following method:

    dBootstrapStatistic{T<:Number}(x::Vector{T}, bp::BootstrapParam)
    
which returns a `Vector{Float64}` with length equal to numResample (a field in `bp`). That is, the output vector contains the re-sampled test statistics of interest. The user could, for example, use this output vector to construct a histogram or kernel density estimate of the distribution of the test statistic of interest.

Some users may also wish to use keyword functionality via the following method:

    dBootstrapStatistic{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003", statistic::ASCIIString="mean")
    
    
###### `dBootstrap`

This function wraps `dBootstrapStatistic`. The only functionality added by `dBootstrap` is to convert the resampled test statistics output by `dBootstrapStatistic` into an estimate of a distribution parameter of the test statistic, e.g. variance, or confidence intervals. Maximum flexibility can be obtained using the following method:

    dBootstrap{T<:Number}(x::Vector{T}, bp::BootstrapParam)
    
The output of this function can vary, although in most cases it will be a `Float64`. However, for example, in the case of confidence intervals, it will be the tuple `(Float64, Float64)`, providing lower and upper confidence bounds. Some users may also wish to use keyword functionality via the following method:

    dBootstrap{T<:Number}(x::Vector{T}; numObsResample::Int=length(x), numResample::Int=defaultNumResample, bootstrapMethod::ASCIIString="stationary", blockLength::Int=-1, blockLengthMethod::ASCIIString="PPW2009", bandwidthMethod::ASCIIString="P2003", statistic::ASCIIString="mean", distributionParam::ASCIIString="variance")
    

###### `dBootstrapBlockLength`

This function is a stand-alone function that can get called at various points by the other functions in the module if the block length is set to a non-positive number. There are multiple methods for this function with inputs `({T<:Number}(x::Vector{T}, bm::BlockLengthMethod)`. Multiple dispatch is used to ensure that the appropriate block length selection procedure is called. As with the other functions in this module, the easiest way to interact with this function is via the method:

    dBootstrapBlockLength{T<:Number}(x::Vector{T}, bp::BootstrapParam)

where `x` is the observable data vector. This will return an estimate of the block length expressed as a `Float64`. Some users may also wish to use keyword functionality via the following method:

    dBootstrapBlockLength{T<:Number}(x::Vector{T}; blockLengthMethod::ASCIIString="PPW2009", bootstrapMethod::ASCIIString="stationary", bandwidthMethod::ASCIIString="P2003")
    
###### In-place functions

All of the above functions, ie `dBootstrapIndex`, `dBootstrapData`, `dBootstrapStatistic`, `dBootstrap`, and `dBootstrapBlockLength`, all contain an "in-place" version, denoted `dBootstrapIndex!`, `dBootstrapData!`, `dBootstrapStatistic!`, `dBootstrap!`, and `dBootstrapBlockLength!`. 

These in-place versions always exhibit the method `{T<:Number}(x::Vector{T}, bp::BootstrapParam)` (or with the arguments supplied in the opposite order). These methods call the regular methods, but the only difference is that *if, and only if* the block length of the input `BootstrapParam` is non-positive, these methods will update in-place the block length in the input `BootstrapParam` with a new estimate of the block length based on the data-driven methods discussed earlier in this document.


###### Other functions

As discussed, `dBootstrapVar`, `dBootstrapStd`, and `dBootstrapConf`, are all just keyword wrappers on `dBootstrap` that automatically initialise the distributionParam field to variance, standard deviation, or confidence intervals respectively.

The only other exported functions are listed now:

* `getBlockLength(bp::BootstrapParam)`. Returns block-length in input `bp` as a `Float64`.
* `getBlockLength(bm::BootstrapMethod)`. Returns block-length in input bootstrap method subtype as a `Float64`.
* `replaceBlockLength!{T<:Number}(bp::BootstrapParam, newBlockLength{T})`. Replace the block-length in `bp` with `newBlockLength`.
* `replaceBlockLength!{T<:Number}(bm::BootstrapMethod, newBlockLength{T})`. Replace the block-length in `bm` with `newBlockLength`.
* `replaceNumObsData!(bp::BootstrapParam, newNumObsData::Int)`. Replace the field numObsData in `bp` with the new value `newNumObsData`.
* `replaceNumObsResample!(bp::BootstrapParam, newNumObsResample::Int)`. As above, but for the field numObsResample.
* `replaceNumResample!(bp::BootstrapParam, newNumResample::Int)`. As above, but for the field numResample.


This concludes the documentation for DependentBootstrap. If you have any other questions, don't hesitate to ask. If you find any bugs, please submit an issue. Any constructive recommendations would be most welcome, *especially* any places in which routines can be optimized for faster run-time.
