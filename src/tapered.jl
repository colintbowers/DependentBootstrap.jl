
#Local function used to transform data via appropriate influence functions for the case where bootstrap method is tapered block
function apply_influence_function(x::Vector{T}, bi::BootInput)::Vector{Float64} where {T<:Number}
	if bi.flevel1 == mean
		x_if = x - mean(x, 1)
	elseif bi.flevel1 == sum
		x_if = length(x) * (x - mean(x))
	else
		error("Tapered block bootstrap only implemented for a limited number of cases with known influence functions")
	end
	return(x_if)
end
(apply_influence_function(x::Vector{Vector{T}}, bi::BootInput)::Vector{Vector{Float64}}) where {T<:Number} = Vector{Float64}[ apply_influence_function(x[k], bi) for k = 1:length(x) ]

#Local function used to weight data for the case where bootstrap method is tapered block
function dboot_kernel_weights(bi::BootInput)::Vector{Float64}
    bL = Int(ceil(bi.blocklength))
	kernelInput = Float64[ (1 / bL) * (n - 0.5) for n = 1:bL ]
    kernelWeight = kernel_func_pp2002(kernelInput, bi.bootmethod)
	normTerm = sqrt(bL) / norm(kernelWeight, 2)
    kernelWeight .*= normTerm
	return(kernelWeight)
end
function dboot_weight!(x::Vector{Vector{T}}, bi::BootInput)::Vector{Vector{Float64}} where {T<:Number}
    bL = Int(ceil(bi.blocklength))
    bL <= 1 && return(x)
    w = dboot_kernel_weights(bm)
	length(w) != bL && error("Logic fail. Incorrect length output from dboot_kernel_weights function. Please file an issue.")
    (num_repeat, num_remain) = divrem(size(x, 1), bL)
    wLong = vcat(repeat(w, outer=num_repeat), w[1:num_remain])
    for n = 1:length(x)
        x[n] = wLong .* x[n]
    end
	return(x)
end
#Local function for the two kernel functions proposed in Paparoditis and Politis (2002) "The tapered block bootstrap for general statistics from stationary sequences"
function kernel_func_pp2002_trap(x::Float64)::Float64
    p = 0.43 #Optimal value from PP (2002)
    x < 0 && return(0.0)
    x < p && return(x / p)
    x < 1 - p && return(1.0)
    x < 1 && return((1 - x) / p)
    return(0.0)
end
kernel_func_pp2002_trap(x::Vector{Float64})::Vector{Float64} = Float64[ kernel_func_pp2002_trap(x[k]) for k = 1:length(x) ]
function kernel_func_pp2002_smooth(x::Float64)::Float64
    p = 1.3 #Optimal value from PP (2002)
    (0.0 <= x <= 1.0) && return(1 - abs(2*x - 1)^p)
    return(0.0)
end
kernel_func_pp2002_smooth(x::Vector{Float64})::Vector{Float64} = Float64[ kernel_func_pp2002_smooth(x[k]) for k = 1:length(x) ]
function kernel_func_pp2002(x::Vector{Float64}, bm::BootTapered)::Vector{Float64}
    bm.kernelfunction == :trap && return(kernel_func_pp2002_trap(x))
    bm.kernelfunction == :smooth && return(kernel_func_pp2002_smooth(x))
    error("Invalid kernel function of $(bm.kernelfunction)")
end
