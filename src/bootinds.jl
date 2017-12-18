
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
dbootinds_one(bi::BootInput)::Vector{Int} \n
dbootinds_one(data::T; kwargs...)::Vector{Int} \n

Returns a single resampling index that, when used to index the original dataset,
will provide a single resampled dataset. \n
A keyword method that calls the keyword constructor for BootInput is also provided.
"""
dbootinds_one(bi::BootInput, bm::BootIID)::Vector{Int} = rand(1:bi.numobs, bi.numobsperresample)
function dbootinds_one(bi::BootInput, bm::BootStationary)::Vector{Int}
    bi.blocklength <= 1.0 && return dbootinds_one(bi, BootIID())
    inds = Array{Int}(bi.numobsperresample)
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
function dbootinds_one(bi::BootInput, bm::BootMoving)::Vector{Int}
    bl = Int(ceil(bi.blocklength))
    bl == 1 && return dbootinds_one(bi, BootIID())
    inds = Array{Int}(bi.numobsperresample)
    blockstart_ub = max(1, bi.numobs-bl+1)
    for n = 1:bl:bi.numobsperresample
        inds[n] = rand(1:blockstart_ub) #Start of block
        for s = n+1:min(n+bl-1, bi.numobsperresample) #Iterate through block (use of min avoids bounds error)
            inds[s] = inds[s-1] + 1
        end
    end
	return inds
end
function dbootinds_one(bi::BootInput, bm::BootCircular)::Vector{Int}
    bl = Int(ceil(bi.blocklength))
    bl == 1 && return(dbootinds_one(bi, BootIID()))
    inds = Array{Int}(bi.numobsperresample)
    for n = 1:bl:bi.numobsperresample
        inds[n] = rand(1:bi.numobs) #Start of block
        for s = n+1:min(n+bl-1, bi.numobsperresample) #Iterate through block (use of min avoids bounds error)
            inds[s-1] == bi.numobs ? (inds[s] = 1) : (inds[s] = inds[s-1] + 1)
        end
    end
	return inds
end
function dbootinds_one(bi::BootInput, bm::BootNoOverlap)::Vector{Int}
    bl = Int(ceil(bi.blocklength))
    bl == 1 && return(dbootinds(bi, BootIID()))
    inds = Array{Int}(bi.numobsperresample)
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
#dbootinds_one(bi::BootInput, bm::BootTapered)::Vector{Int} = dbootinds_one(bi, BootMoving())
dbootinds_one(bi::BootInput, bm::BootTapered)::Vector{Int} = error("The tapered block bootstrap is currently not available in this package. Users interested in contributing should check the package github page.")
dbootinds_one(bi::BootInput)::Vector{Int} = dbootinds_one(bi, bi.bootmethod)
(dbootinds_one(data::T, bi::BootInput)::Vector{Int}) where {T} = dbootinds_one(bi) #included for consistency
function dbootinds_one(data::T ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                    blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=numobs)::Vector{Vector{Int}} where {T}
    return dbootinds_one(BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=flevel2, numobsperresample=numobsperresample))
end
"""
dbootinds(bi::BootInput)::Vector{Vector{Int}} \n
dbootinds(data::T; kwargs...)::Vector{Vector{Int}} \n

Each inner vector of the returned Vector{Vector{Int}} provides indices that, when used to
index the original dataset, will provide a single resampled dataset. \n
A keyword method that calls the keyword constructor for BootInput is also provided. \n
"""
function dbootinds(bi::BootInput, bm::T)::Vector{Vector{Int}} where {T<:BootMethod}
    !check_bi_for_dbootvecinds(bi) && error("Logic fail in check_bi_for_dbootvecinds. Please file an issue.")
    return [ dbootinds_one(bi, bm)  for n = 1:bi.numresample ]
end
dbootinds(bi::BootInput)::Vector{Vector{Int}} = dbootinds(bi, bi.bootmethod)
(dbootinds(data::T, bi::BootInput)::Vector{Vector{Int}}) where {T} = dbootinds(bi) #This method is included for consistency
function dbootinds(data::T ; blocklength::Number=0.0, numresample::Number=NUM_RESAMPLE, bootmethod::Symbol=:stationary,
                    blmethod::Symbol=:dummy, flevel1::Function=mean, flevel2::Function=var, numobsperresample::Number=numobs)::Vector{Vector{Int}} where {T}
    return dbootinds(BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blmethod=blmethod, flevel1=flevel1, flevel2=flevel2, numobsperresample=numobsperresample))
end
