
#Local function for checking a BootInput is valid for a call to dbootvecinds
function check_bi_for_dbootvecinds(bi::BootInput)::Bool
    bi.numobs < 1 && error("Number of observations must be strictly positive")
    bi.numresample < 1 && error("Number of resamples must be strictly positive")
    bi.numobsperresample < 1 && error("Number of observations per resample must be strictly positive")
    bi.blocklength <= 0.0 && error("Block length must be strictly positive")
    isnan(bi.blocklength) && error("Block length is set to NaN")
    isinf(bi.blocklength) && error("Block length is infinite")
    return true
end

"""
dbootinds_one(bi::BootInput)::Vector{Int}
dbootinds_one(data::T; kwargs...)::Vector{Int}

Returns a single resampling index that, when used to index the original dataset,
will provide a single resampled dataset.

A keyword method that calls the keyword constructor for BootInput is also provided. Please
use ?BootInput at the REPL for more detail on feasible keywords.
"""
(dbootinds_one(bi::BootInput{BootIID})::Vector{Int}) = rand(1:bi.numobs, bi.numobsperresample)
function dbootinds_one(bi::BootInput{BootStationary})::Vector{Int}
    bi.blocklength <= 1.0 && return rand(1:bi.numobs, bi.numobsperresample)
    inds = zeros(Int, bi.numobsperresample)
    geo1 = Geometric(1 / bi.blocklength)
    (c, geodraw) = (1, 1)
    for n = 1:bi.numobsperresample
        if c == geodraw #Start a new block
            inds[n] = rand(1:bi.numobs)
            geodraw = rand(geo1) + 1
            c = 1
        else #Next obs in existing block
            inds[n-1] == bi.numobs ? (inds[n] = 1) : (inds[n] = inds[n-1] + 1)
            c += 1
        end
    end
    return inds
end
function dbootinds_one(bi::BootInput{BootMoving})::Vector{Int}
    bl = ceil(Int, bi.blocklength)
    bl == 1 && return rand(1:bi.numobs, bi.numobsperresample)
    inds = zeros(Int, bi.numobsperresample)
    blockstart_ub = max(1, bi.numobs-bl+1)
    for n = 1:bl:bi.numobsperresample
        inds[n] = rand(1:blockstart_ub) #Start of block
        for s = n+1:min(n+bl-1, bi.numobsperresample) #Iterate through block (use of min avoids bounds error)
            inds[s] = inds[s-1] + 1
        end
    end
	return inds
end
function dbootinds_one(bi::BootInput{BootCircular})::Vector{Int}
    bl = ceil(Int, bi.blocklength)
    bl == 1 && return rand(1:bi.numobs, bi.numobsperresample)
    inds = zeros(Int, bi.numobsperresample)
    for n = 1:bl:bi.numobsperresample
        inds[n] = rand(1:bi.numobs) #Start of block
        for s = n+1:min(n+bl-1, bi.numobsperresample) #Iterate through block (use of min avoids bounds error)
            inds[s-1] == bi.numobs ? (inds[s] = 1) : (inds[s] = inds[s-1] + 1)
        end
    end
	return inds
end
function dbootinds_one(bi::BootInput{BootNoOverlap})::Vector{Int}
    bl = ceil(Int, bi.blocklength)
    bl == 1 && return rand(1:bi.numobs, bi.numobsperresample)
    inds = zeros(Int, bi.numobsperresample)
	blockstart_ub = max(1, bi.numobs-bl+1)
	blockstartvalues = collect(1:bl:blockstart_ub) #Build valid set of start indices for any block
    length(blockstartvalues) == 1 && error("Not enough observations to perform non-overlapping block bootstrap given block length: num obs = $(bi.numobs), block length = $(bl)")
    (blockstartvalues[end] + bl - 1 > bi.numobs) && pop!(blockstartvalues)
    for n = 1:bl:bi.numobsperresample
        inds[n] = blockstartvalues[rand(1:length(blockstartvalues))] #Start of block
        for s = n+1:min(n+bl-1, bi.numobsperresample) #Iterate through block (use of min avoids bounds error)
            inds[s] = inds[s-1] + 1
        end
    end
	return inds
end
#(dbootinds_one(bi::BootInput{BootTapered})::Vector{Int}) = error("Routines for the tapered block bootstrap are currently not completed. Users interested in contributing should check the package github page.")
(dbootinds_one(bi::BootInput{BootDummy})::Vector{Int}) = error("Logic fail. It should not have been possible to call this method. Please file an issue with full stacktrace.")
(dbootinds_one(data, bi::BootInput)::Vector{Int}) = dbootinds_one(bi)
(dbootinds_one(data ; kwargs...)::Vector{Int}) = dbootinds_one(data, BootInput(data ; kwargs...))

"""
    dbootinds(data::T ; bi::BootInput)::Vector{Vector{Int}}
    dbootinds(data::T ; kwargs...)::Vector{Vector{Int}}

Each inner vector of the returned Vector{Vector{Int}} provides indices that, when used to
index the original dataset, will provide a single resampled dataset.

A keyword method that calls the keyword constructor for BootInput is also provided. Please
use ?BootInput at the REPL for more detail on feasible keywords.

Please use dbootinds_one if you only want to obtain a single Vector{Int} resampling index.
"""
(dbootinds(bi::BootInput)::Vector{Vector{Int}}) = check_bi_for_dbootvecinds(bi) ? [ dbootinds_one(bi) for n = 1:bi.numresample ] : error("Logic fail in check_bi_for_dbootvecinds. Please file an issue.")
(dbootinds(data, bi::BootInput)::Vector{Vector{Int}}) = dbootinds(bi)
(dbootinds(data ; kwargs...)::Vector{Vector{Int}}) = dbootinds(data, BootInput(data ; kwargs...))
