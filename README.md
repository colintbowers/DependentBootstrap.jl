DependentBootstrap.jl
=====================

[![Build Status](https://travis-ci.org/colintbowers/DependentBootstrap.jl.svg?branch=master)](https://travis-ci.org/colintbowers/DependentBootstrap.jl)

A module for the Julia language that implements several varieties of the dependent statistical bootstrap as well as the corresponding block-length selection procedures.

## News

This package is compatible with julia v1.0+. If you are running v0.6, you will need to use `Pkg.pin("DependentBootstrap", v"0.1.1")` at the REPL, and if you are running v0.5, use `Pkg.pin("DependentBootstrap", v"0.0.1")`. Compability with versions before v0.5 is not available.  

## Main features

This module allows Julia users to estimate the distribution of a test statistic using any of several varieties of the dependent bootstrap.

The following bootstrap methods are implemented:
* the *iid* bootstrap proposed in Efron (1979) "Bootstrap Methods: Another Look at the Jackknife",
* the stationary bootstrap proposed in Politis, Romano (1994) "The Stationary Bootstrap"
* the moving block bootstrap proposed in Kunsch (1989) "The jackknife and the bootstrap for general stationary observations" and (independently) Liu, Singh (1992) "Moving blocks jackknife and bootstrap capture weak dependence",
* the circular block bootstrap proposed in Politis, Romano (1992) "A circular block resampling procedure for stationary data", and
* the non-overlapping block bootstrap described in Lahiri (1999) *Resampling Methods for Dependent Data* (this method is not usually used and is included mainly as a curiosity).

The module also implements the following block length selection procedures:
* the block length selection procedure proposed in Politis, White (2004) "Automatic Block Length Selection For The Dependent Bootstrap", including the correction provided in Patton, Politis, and White (2009)

Bandwidth selection for the block length procedures is implemented using the method proposed in Politis (2003) "Adaptive Bandwidth Choice".

Some work has been done to implement the tapered block bootstrap of Paparoditis, Politis (2002) "The tapered block bootstrap for general statistics from stationary sequences", along with corresponding block-length selection procedures, but it is not yet complete.

The module is implemented entirely in Julia.

## What this package does not include

I have not included any procedures for bootstrapping confidence intervals in a linear regression framework, or other parametric models. This functionality is provided by [Bootstrap.jl](https://github.com/juliangehring/Bootstrap.jl), and work is currently under way to add bootstrap methods from this package to the `Bootstrap` API.

I also have not included support for the jackknife, wild bootstrap, or subsampling procedures. I would be quite open to pull requests that add these methods to the present package, but have not had time to implement them myself. Work is ongoing to include the tapered block bootstrap, and ideally, the package will also eventually include the extended tapered block bootstrap. If you are interested in working on any of these projects, please feel free to contact me.

## How to use this package

#### Installation

This package should be added using `using Pkg ; Pkg.add("DependentBootstrap")`, and can then be called with `using DependentBootstrap`. The package depends on `StatsBase` and `Distributions` for some functionality, and on `DataFrames` and `TimeSeries` so that `DataFrame` and `TimeArray` datasets can be supported by this packages methods.

#### Terminology

In what follows, I use the terminology from Lahiri (1999) *Resampling Methods for Dependent Data* and refer to the underlying test statistic of interest as a *level 1 statistic*, and the distribution parameter of the test statistic that is of interest as a *level 2 parameter*. For example, the user might have some dataset `x` of type `T_data`, and be interested in the variance of the sample mean of `x`. In this case, the level 1 statistic is the sample mean function `mean`, and the level 2 parameter is the sample variance function `var`.

I use `T_data` to refer to the type of the users dataset, `T_level1` to refer to the output type obtained by applying the level 1 statistic function to the dataset, and `T_level2` to refer to the output type obtained by applying the level 2 statistic to a `Vector{T_level1}` (i.e. a vector of resampled level 1 statistics).

#### Exported functions

The package exports the following functions, all of which have docstrings that can be called interactively at the REPL:
* `dbootinds(...)::Vector{Vector{Int}}` -> Returns indices that can be used to index into the underlying data to obtain bootstrapped data. Note, each inner vector of the output corresponds to a single re-sample for the underlying data.
* `dbootdata(...)::Vector{T_data}` -> Returns the bootstrapped data. Each element of the output vector corresponds to one re-sampled dataset, and the output vector will have length equal to `numresample` (a parameter discussed later).
* `dbootlevel1(...)::Vector{T_level1}` -> Returns a vector of bootstrapped level 1 statistics, where the output vector will have length equal to `numresample`.
* `dbootlevel2(...)::T_level2` -> Returns the bootstrapped distribution parameter of the level 1 statistic.
* `dboot(...)::T_level2` -> Identical to dbootlevel2. Most users will want to use this function.
* `dbootvar(...)::Float64` -> Identical to `dboot` but automatically sets `flevel2` to `var` (the sample variance function)
* `dbootconf(...)::Vector{Float64}` -> Identical to `dboot` but automatically sets `flevel2` to the anonymous function `x -> quantile(x, [0.025, 0.975])`, so the level 2 distribution parameter is a 95% confidence interval. In addition to the usual keywords, the keyword version of this function also accepts the keyword `alpha::Float64=0.05`, which controls the width of the confidence interval. Note, `0.05` corresponds to a 95% confidence interval, `0.1` to a 90% interval, and `0.01` to a 99% interval (and so on).
* `optblocklength(...)::Float64` -> Returns the optimal block length.

The function `bandwidth_politis_2003{T<:Number}(x::AbstractVector{T})::Tuple{Int, Float64, Vector{Float64}}` is not exported, but the docstrings can be accessed using `?DependentBootstrap.bandwidth_politis_2003` at the REPL. This function implements the bandwidth selection procedure from Politis (2003) discussed above, and may be of independent interest to some users.

All of the above functions exhibit the following two core methods:

* `f(data ; kwargs...)`
* `f(data, bi::BootInput)`

where `data` is the users underlying dataset, `kwargs` is a collection of keyword arguments, and `bi::BootInput` is a core type exported by the module that will be discussed later (but can be safely ignored by most users). The following types for `data` are currently accepted:
* `Vector{<:Number}`,
* `Matrix{<:Number}` where rows are observations and columns are variables,
* `Vector{Vector{<:Number}}` where each inner vector is a variable,
* `DataFrame`
* `TimeArray`

Of the two core methods, most users will want the `kwargs` method. A list of valid keyword arguments and their default values follows:
* `blocklength`         <- Block length for bootstrapping procedure. The default value is `0`. Set to <= 0 to auto-estimate the optimal block length from the dataset. `Float64` inputs are allowed.
* `numresample`         <- Number of times to resample the input dataset. The default value is the module constant `NUM_RESAMPLE`, currently set to `1000`.
* `bootmethod`          <- Bootstrapping methodology to use. The default value is `:stationary` (for the stationary bootstrap).
* `blocklengthmethod`   <- Block length selection procedure to use if user wishes to auto-estimate the block length. Default value is `:ppw2009` (use the method described in Patton, Politis, and White (2009)).
* `flevel1`             <- A function that converts the input dataset to the estimator that the user wishes to bootstrap. The default value is `mean`.
* `flevel2`             <- A function that converts a vector of estimators constructed by `flevel1` into a distributional parameter. The default value is `var`.
* `numobsperresample`   <- Number of observations to be drawn (with replacement) per resample. The default value is the number of observations in the dataset (the vast majority of users will want this default value).
* `fblocklengthcombine` <- A function for converting a `Vector{Float64}` of estimated blocklengths to a single `Float64` blocklength estimate, which is necessary when the input dataset is a multivariate type. The default value is `median`.

A list of acceptable keyword arguments for `bootmethod` and `blocklengthmethod` follows. Note you can use either `String` or `Symbol` when specifying these arguments. For `bootmethod` we have:
* `:iid` or `:efron`                     <- IID bootstrap
* `:stationary`                          <- Stationary bootstrap
* `:movingblock` or `:moving`            <- Moving block bootstrap
* `:nonoverlappingblock` or `:nooverlap` <- Nonoverlapping block bootstrap
* `:circularblock` or `circular`         <- Circular block bootstrap

For `blocklengthmethod` we have:
* `:ppw2009` <- Block length selection method of Patton, Politis, and White (2009)

Acceptable arguments can also be examined interactively by examining the keys of the module dictionaries `BOOT_METHOD_DICT` and `BLOCKLENGTH_METHOD_DICT`.

In practice, the keyword argument method `f(data ; kwargs...)` actually just wraps a call to `f(data, BootInput(kwargs...))` under the hood. However, most users will not need to concern themselves with this level of detail.

For those who wish more fine-grained control, please use `?BootInput` at the REPL to get more information on this core module type.

#### Examples

Let `data::Vector{Float64}`.

The variance of the sample mean of `data` can be bootstrapped using a stationary bootstrap with optimally estimated block length using `dboot(data)` or `dbootvar(data)`.

A 90% confidence interval for the sample median using a circular block bootsrap with block length of 5 can be estimated using `dboot(data, blocklength=5, bootmethod=:circular, flevel1=median, flevel2=(x -> quantile(x, [0.05, 0.95])))` or `dbootconf(data, blocklength=5, bootmethod=:circular, flevel1=median, alpha=0.1)`.

Moving block bootstrap indices for generating bootstrapped data with optimally estimated block length can be obtained using `dbootinds(data, bootmethod=:moving)`, or if the user wants the bootstrapped data not the indices, `dbootdata(data, bootmethod=:moving)`. If the user wants bootstrapped sample medians of `data`, then use `dbootlevel1(data, bootmethod=:moving, flevel1=median)`.

If the user wants the optimal block length using the method proposed in Patton, Politis, and White (2009), use `optblocklength(data, blmethod=:ppw2009)`.

Now let `data::Matrix{Float64}`.

If the user wants the median optimal block length from each column of `data`, use `optblocklength(data, blmethod=:ppw2009)`. If the user wants the average optimal block length use `optblocklength(data, blmethod=:ppw2009, fblocklengthcombine=mean)`.

If the user wants the median of the test statistic that is the maximum of the sample mean of each column, using a stationary bootstrap with optimal block length, then use `dboot(data, flevel1=(x -> maximum(mean(x, dims=1))), flevel2=median)`. If `data::Vector{Vector{Float64}}` instead, and the user wanted the 95% confidence interval, use `dbootconf(data, flevel1=(x -> maximum([ mean(x[k]) for k = 1:length(x) ])))`.
