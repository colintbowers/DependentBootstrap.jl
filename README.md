DependentBootstrap.jl
=====================

[![Build Status](https://travis-ci.org/colintbowers/DependentBootstrap.jl.svg?branch=master)](https://travis-ci.org/colintbowers/DependentBootstrap.jl)

A module for the Julia language that implements several varieties of the dependent statistical bootstrap as well as the corresponding block-length selection procedures.

## WARNING

I have recently upgraded and simplified this package for use with Julia v0.5. Unfortunately, this update includes an enormous number of breaking changes. On the plus side, the package is currently in the process of becoming officially registered, and after this, any breaking changes will come with deprecation warnings.

## Main features

This module allows Julia users to estimate the distribution of a test statistic using any of several varieties of the dependent bootstrap.

The following bootstrap methods are implemented:
* the *iid* bootstrap proposed in Efron (1979) "Bootstrap Methods: Another Look at the Jackknife",
* the stationary bootstrap proposed in Politis, Romano (1994) "The Stationary Bootstrap"
* the moving block bootstrap proposed in Kunsch (1989) "The jackknife and the bootstrap for general stationary observations" and (independently) Liu, Singh (1992) "Moving blocks jackknife and bootstrap capture weak dependence",
* the circular block bootstrap proposed in Politis, Romano (1992) "A circular block resampling procedure for stationary data", and
* the non-overlapping block bootstrap described in Lahiri (1999) *Resampling Methods for Dependent Data* (this method is not usually used and is included mainly as a curiosity).

Some work has been done to implement the tapered block bootstrap of Paparoditis, Politis (2002) "The tapered block bootstrap for general statistics from stationary sequences", but, unfortunately it is not yet complete.

The module also implements the following block length selection procedures:
* the block length selection procedure proposed in Politis, White (2004) "Automatic Block Length Selection For The Dependent Bootstrap", including the correction provided in Patton, Politis, and White (2009), and
* the block length selection procedure proposed in Paparoditis, Politis (2002) "The tapered block bootstrap for general statistics from stationary sequences".

Bandwidth selection for these block length procedures is implemented using the method proposed in Politis (2003) "Adaptive Bandwidth Choice".

The module is implemented entirely in Julia. Each bootstrap procedure and block-length selection procedure is implemented as its own type. Various functions then operate on these types via multiple dispatch in order to return bootstrap indices, data, statistics, or distribution parameter estimates, as well as optimal block length estimates. These functions operate on both univariate and multivariate datasets.

## What this package does not include

I have not included any procedures for bootstrapping confidence intervals in a linear regression framework. I believe that this functionality is better provided by a separate package (possibly GLM), that can use this package for the bootstrapping step.

I also have not included support for the jackknife, wild bootstrap, or subsampling procedures. I would be quite open to pull requests that add these methods to the present package, but have not had time to implement them myself. Work is ongoing on including the tapered block bootstrap, and ideally, the package will also eventually include the extended tapered block bootstrap.

## How to use this package

#### Installation

This package should be added using `Pkg.add("DependentBootstrap")`, and can then be called with `using DependentBootstrap`. The package only has two dependencies (currently): StatsBase and Distributions. Support for DataFrames may be added in the future, depending on demand (the package is structured so that this will be a relatively trivial addition)

#### Terminology

In what follows, I use the terminology from Lahiri (1999) *Resampling Methods for Dependent Data* and refer to the underlying test statistic of interest as a *level 1 statistic*, and the distribution parameter of the test statistic that is of interest as a *level 2 parameter*. For example, the user might be interested in the variance of the sample mean estimator, in which case the level 1 statistic is the sample mean, and the level 2 parameter is the variance.

#### Core Type

The package implements the core type `BootInput`. This type includes the following fields (in this order):
* `numobs::Int` -> Number of observations (per column) in the dataset that is to be bootstrapped
* `blocklength::Float64` -> Block length for the  bootstrap procedure
* `numresample::Int` -> Number of bootstrap resamples to use
* `bootmethod::BootMethod` -> Indicates which bootstrap methodology to use
* `blmethod::BLMethod` -> Indicates which block length selection procedure to use (if blocklength <= 0.0)
* `flevel1::Function` -> The level 1 statistic of interest
* `flevel2::Function` -> The level 2 parameter of interest
* `numobsperresample::Int` -> Number of observations per resample (almost always equal to numobs)

For convenience, a keyword constructor for `BootInput` is provided and the function signature follows:

`BootInput(data ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary, blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=data_length(data))`

Some remarks are useful here:

* `data` is the users dataset. Currently accepted types for `data` are `Vector{T<:Number}`, `Vector{Vector{T<:Number}}`, and `Matrix{T<:Number}`. Additional types for `data` are easy to add and so interested users should open an issue on github. In what follows, the type of `data` will be written as `T_data`
* The default value for `blocklength` is 0.0. Any input `blocklength` <= 0.0 tells the `BootInput` constructor to estimate the optimal block length using the input `data` and the supplied block length method, or, if no block length method is supplied, to choose the best available block length procedure for the given bootstrap methodology.
* The default value for numresample is the module constant `NUM_RESAMPLE` which is currently `1000`.
* `bootmethod` can be specified using a symbol, which is converted to the appropriate type (which will be a subtype of the abstract `BootMethod`) by the constructor. Valid values are `:iid`, `:stationary`, `:moving`, `:circular`, and `:nooverlap`. Hopefully `:tapered` will be available in future iterations of this package.
* `blmethod` can be specified using a symbol, which is converted to the appropriate type (which will be a subtype of the abstract `BLMethod`) by the constructor. The default value of `:dummy` tells the constructor to choose the most appropriate block length selection procedure for the given `bootmethod`. Note, if the input `blocklength` is > 0.0, then `blmethod` is irrelevant. Valid symbols for `blmethod` are `:ppw2009` and `:pp2002`. It is anticipated that the vast majority of users will want to use `:ppw2009`, and this will be auto-selected for all bootstrap methods, except `:tapered`, which, when operational, will use `:pp2002`.
* `flevel1` should be a function that accepts an input of type `T_data`, and can return whatever type the user wants (note, it does not need to be scalar). In what follows, the output type of `flevel1` will be written as `T_level1`. An informative error is thrown if `flevel1` cannot accept `T_data` as input.
* `flevel2` should be a function that accepts an input of type `T_level1`, and can return whatever type the user wants (again, it need not be scalar). In what follows, the output type of `flevel2` will be written as `T_level2`. An informative error is thrown if `flevel2` cannot accept `T_level1` as input.

Several fields of a `BootInput` can be modified using the following exported functions:
* `setblocklength!(b::BootInput, blocklength::Number)`
* `setnumresample!(b::BootInput, numresample::Number)`
* `setflevel1!(b::BootInput, flevel1::Function)`
* `setflevel2!(b::BootInput, flevel2::Function)`

If you wish to modify any other fields, it is better just to construct a new `BootInput` (a near instantaneous procedure).

#### Functions

All core functions exported by this package include (at least) the following core function signature:

`f(data, bi::BootInput)`

where `data` is the users underlying dataset. They also incldue the following keyword function signature:

`f(data ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary, blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=data_length(data))`

where this method simply wraps the keyword constructor for `BootInput` in the first step and then calls the core function signature above.
Several exported functions also have other methods, which users are welcome to explore using Julia's `methods` function, or by using docstrings (i.e. type `?f` at the REPL, where `f` is the function name). However, it is envisaged that most users will use one of the two signatures described above.

The package exports the following core functions:
* `dbootinds(...)::Vector{Vector{Int}}` -> Returns indices that can be used to index into the underlying data to obtain bootstrapped data. Note, each inner vector of the output corresponds to a single re-sample for the underlying data.
* `dbootdata(...)::Vector{T_data}` -> Returns the bootstrapped data. Each element of the output vector corresponds to one re-sampled dataset, and the output vector will have length equal to `numresample`.
* `dbootlevel1(...)::Vector{T_level1}` -> Returns a vector of bootstrapped level 1 statistics, where the output vector will have length equal to `numresample`
* `dbootlevel2(...)::T_level2` -> Returns the bootstrapped distribution parameter of the level 1 statistic.
* `dboot(...)::T_level2` -> Identical to dbootlevel2
* `dbootvar(...)::Float64` -> Identical to `dboot` but automatically sets `flevel2` to `var` (the sample variance function)
* `dbootconf(...)::Vector{Float64}` -> Identical to `dboot` but automatically sets `flevel2` to the anonymous function `x -> quantile(x, [0.025, 0.975])`, so the level 2 distribution parameter is a 95% confidence interval. In addition to the usual keywords, the keyword version of this function also accepts the keyword `alpha::Float64=0.05`, which controls the width of the confidence interval. Note, `0.05` corresponds to a 95% confidence interval, `0.1` to a 90% interval, and `0.01` to a 99% interval (and so on).
* `optblocklength(...)::Float64` -> Returns the optimal block length. If `data` is a multivariate dataset, then the core function signatures include an optional additional input `[, f::Function=median]`. In this situation, the optimal block length is calculated for each column of the input `data` and then `f` is applied to reduce this `Vector{Float64}` to a single `Float64` optimal block length estimate. As can be seen, the default transformation is `median`.

The function `bandwidth_politis_2003{T<:Number}(x::AbstractVector{T})::Tuple{Int, Float64, Vector{Float64}}` is not exported, but the docstrings can be accessed using `?DependentBootstrap.bandwidth_politis_2003` at the REPL. This function implements the bandwidth selection procedure from Politis (2003) discussed above, and may be of independent interest to some users.

#### Examples

Let `data::Vector{Float64}`.

The variance of the sample mean of `data` can be bootstrapped using a stationary bootstrap with optimally estimated block length using `dboot(data)` or `dbootvar(data)`.

A 90% confidence interval for the sample median using a circular block bootsrap with block length of 5 can be estimated using `dboot(data, blocklength=5, bootmethod=:circular, flevel1=median, flevel2=(x -> quantile(x, [0.05, 0.95])))` or `dbootconf(data, blocklength=5, bootmethod=:circular, flevel1=median, alpha=0.1)`.

Moving block bootstrap indices for generating bootstrapped data with optimally estimated block length can be obtained using `dbootinds(data, bootmethod=:moving)`, or if the user wants the bootstrapped data not the indices, `dbootdata(data, bootmethod=:moving)`. If the user wants bootstrapped sample medians of `data`, then use `dbootlevel1(data, bootmethod=:moving, flevel1=median)`.

If the user wants the optimal block length using the method proposed in Patton, Politis, and White (2009), use `optblocklength(data, blmethod=:ppw2009)`.

Now let `data::Matrix{Float64}`.

If the user wants the average optimal block length from each column of `data`, use `optblocklength(x, mean, blmethod=:ppw2009, )`.

If the user wants the median of the test statistic that is the maximum of the sample mean of each column, using a stationary bootstrap with optimal block length, then use `dboot(data, flevel1=(x -> maximum(mean(x, 1))), flevel2=median)`. If `data::Vector{Vector{Float64}}` instead, and the user wanted the 95% confidence interval, use `dbootconf(data, flevel1=(x -> maximum([ mean(x[k]) for k = 1:length(x) ])))`.
