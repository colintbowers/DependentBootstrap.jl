
#Abstract supertypes
"BootMethod <- Abstract supertype for all dependent bootstrap methods"
abstract type BootMethod end
"BlockLengthMethod <- Abstract supertype for all block length selection methods"
abstract type BlockLengthMethod end
"BandwidthMethod <- Abstract supertype for all bandwidth selection methods"
abstract type BandwidthMethod end
"KernelFunctionMethod <- Abstract supertype for all kernel function methods"
abstract type KernelFunctionMethod end

#Kernel function types
"KernelDummy <- Dummy type for kernel functions"
struct KernelDummy <: KernelFunctionMethod ; end
struct KernelTrap <: KernelFunctionMethod ; end
struct KernelSmooth <: KernelFunctionMethod ; end
Base.show(io::IO, x::KernelDummy) = print(io, "Dummy kernel function method")
Base.show(io::IO, x::KernelTrap) = print(io, "Trapezoidal kernel function method from Paparoditis, Politis (2002)")
Base.show(io::IO, x::KernelSmooth) = print(io, "Smooth kernel function method from Paparoditis, Politis (2002)")
kernel_function_dict_input = Pair{Symbol,KernelFunctionMethod}[
	:dummy => KernelDummy(),
	:trap => KernelTrap(),
	:smooth => KernelSmooth()
]::Vector{Pair{Symbol,KernelFunctionMethod}}
"KERNEL_FUNCTION_DICT <- Dictionary for converting string or symbol inputs into kernel function types"
const KERNEL_FUNCTION_DICT = Dict{Union{Symbol,String},KernelFunctionMethod}()::Dict{Union{Symbol,String},KernelFunctionMethod}
sizehint!(KERNEL_FUNCTION_DICT, 2*length(kernel_function_dict_input))
for kf in kernel_function_dict_input ; KERNEL_FUNCTION_DICT[kf[1]] = kf[2] ; end
for kf in kernel_function_dict_input ; KERNEL_FUNCTION_DICT[string(kf[1])] = kf[2] ; end

#Bootstrap method types
"BootDummy <- Dummy type used within the module. Should never be seen by the end user"
struct BootDummy <: BootMethod ; end
"BootIID <- Type for using multiple dispatch to get the IID boostrap"
struct BootIID <: BootMethod ; end
"BootStationary <- Type for using multiple dispatch to get the stationary bootstrap"
struct BootStationary <: BootMethod ; end
"BootMoving <- Type for using multiple dispatch to get the moving blocks bootstrap"
struct BootMoving <: BootMethod ; end
"BootNoOverlap <- Type for using multiple dispatch to get the non-overlapping blocks bootstrap"
struct BootNoOverlap <: BootMethod ; end
"BootCircular <- Type for using multiple dispatch to get the circular blocks bootstrap"
struct BootCircular <: BootMethod ; end
# "BootTapered <- Type for using multiple dispatch to get the tapered block bootstrap"
# struct BootTapered{Tkf<:KernelFunctionMethod} <: BootMethod
# 	kernelfunction::Tkf
# 	function BootTapered(kf::Tkf) where {Tkf}
# 		new{Tkf}(kf)
# 	end
# end
# BootTapered()::BootTapered = BootTapered(KernelTrap())
Base.show(io::IO, bm::BootDummy) = print(io, "Dummy bootstrap method")
Base.show(io::IO, bm::BootIID) = print(io, "IID bootstrap")
Base.show(io::IO, bm::BootStationary) = print(io, "Stationary bootstrap")
Base.show(io::IO, bm::BootMoving) = print(io, "Moving block bootstrap")
Base.show(io::IO, bm::BootNoOverlap) = print(io, "Non-overlapping block bootstrap")
Base.show(io::IO, bm::BootCircular) = print(io, "Circular block bootstrap")
#Base.show(io::IO, bm::BootTapered) = print(io, "Tapered block bootstrap")
boot_method_dict_input = Pair{Symbol, BootMethod}[
	:iid => BootIID(),
	:efron => BootIID(),
	:stationary => BootStationary(),
	:movingblock => BootMoving(),
	:moving => BootMoving(),
	:nonoverlappingblock => BootNoOverlap(),
	:nooverlap => BootNoOverlap(),
	:circularblock => BootCircular(),
	:circular => BootCircular()]::Vector{Pair{Symbol,BootMethod}}
"BOOT_METHOD_DICT <- Dictionary for converting string or symbol inputs into bootstrap methods"
const BOOT_METHOD_DICT = Dict{Union{Symbol,String},BootMethod}()::Dict{Union{Symbol,String},BootMethod}
sizehint!(BOOT_METHOD_DICT, 2*length(boot_method_dict_input))
for bm in boot_method_dict_input ; BOOT_METHOD_DICT[bm[1]] = bm[2] ; end
for bm in boot_method_dict_input ; BOOT_METHOD_DICT[string(bm[1])] = bm[2] ; end

#bandwidth method types
"P2003 <- Type for using multiple dispatch to get the bandwidth selection procedure of Politis (2003)"
struct P2003 <: BandwidthMethod ; end
Base.show(io::IO, bw::P2003) = print(io, "Bandwidth selection of Politis (2003)")
bandwidth_method_dict_input = Pair{Symbol,BandwidthMethod}[
	:politis2003 => P2003(),
	:p2003 => P2003()
]::Vector{Pair{Symbol,BandwidthMethod}}
"BANDWIDTH_METHOD_DICT <- Dictionary for converting string or symbol inputs into bandwidth methods"
const BANDWIDTH_METHOD_DICT = Dict{Union{Symbol,String},BandwidthMethod}()::Dict{Union{Symbol,String},BandwidthMethod}
sizehint!(BANDWIDTH_METHOD_DICT, 2*length(bandwidth_method_dict_input))
for bw in bandwidth_method_dict_input ; BANDWIDTH_METHOD_DICT[bw[1]] = bw[2] ; end
for bw in bandwidth_method_dict_input ; BANDWIDTH_METHOD_DICT[string(bw[1])] = bw[2] ; end

#blocklength method types
"BLDummy <- Dummy type for block length method"
struct BLDummy <: BlockLengthMethod ; end
# "BLPP2002 <- Type for using multiple dispatch to get the block length selection procedure of Paparoditis and Politis (2002)"
# struct BLPP2002{Tbw<:BandwidthMethod,Tkf<:KernelFunctionMethod} <: BlockLengthMethod
# 	bandwidthmethod::Tbw
# 	kernelfunction::Tkf
# end
"BLPPW2009 <- Type for using multiple dispatch to get the block length selection procedure of Patton, Politis, and White (2009)"
struct BLPPW2009{Tbw<:BandwidthMethod} <: BlockLengthMethod
	bandwidthmethod::Tbw
end
BLPPW2009()::BLPPW2009{P2003} = BLPPW2009(P2003())
# Base.show(io::IO, bl::BLPP2002) = print(io, "Block length selection of Paparoditis and Politis (2002)")
Base.show(io::IO, bl::BLPPW2009) = print(io, "Block length selection of Patton, Politis, and White (2009)")
blocklength_method_dict_input = Pair{Symbol,BlockLengthMethod}[
	:ppw2009 => BLPPW2009(P2003())]::Vector{Pair{Symbol,BlockLengthMethod}}
"BLOCKLENGTH_METHOD_DICT <- Dictionary for converting string or symbol inputs into bandwidth methods"
const BLOCKLENGTH_METHOD_DICT = Dict{Union{Symbol,String},BlockLengthMethod}()::Dict{Union{Symbol,String},BlockLengthMethod}
sizehint!(BLOCKLENGTH_METHOD_DICT, 2*length(blocklength_method_dict_input))
for bl in blocklength_method_dict_input ; BLOCKLENGTH_METHOD_DICT[bl[1]] = bl[2] ; end
for bl in blocklength_method_dict_input ; BLOCKLENGTH_METHOD_DICT[string(bl[1])] = bl[2] ; end

"bootmethod_to_blocklengthmethod <- Convert a bootstrap method to the most appropriate block length method"
(bootmethod_to_blocklengthmethod(bm::T)::BLPPW2009{P2003}) where {T<:Union{BootStationary,BootMoving,BootCircular,BootNoOverlap,BootIID}} = BLPPW2009(P2003())
#(bootmethod_to_blocklengthmethod(bm::BootTapered)::BLPP2002{P2003,KernelDummy}) = BLPP2002(P2003(), KernelDummy())


"""
	BootInput

Core type that defines all parameters needed to perform a bootstrap procedure. The
vast majority of users should use the keyword argument constructor that has the
method signature:

    BootInput(data ; kwargs...)

where data is the dataset to be bootstrapped, and kwargs denotes a set of keyword arguments
(defined below) that are used for every exported function in the DependentBootstrap package.
The following keyword arguments and default values follow: \n
    - blocklength         <- Block length for bootstrapping procedure. Default value is 0. Set to <= 0 to auto-estimate the optimal block length from the dataset. Float64 inputs allowed.
	- numresample         <- Number of times to resample the input dataset. Default value is the module constant NUM_RESAMPLE, currently set to 1000.
    - bootmethod          <- Bootstrapping methodology to use. Default value is the Symbol :stationary (for the stationary bootstrap).
	- blocklengthmethod   <- Block length selection procedure to use if user wishes to auto-estimate the block length. Default value is the Symbol :ppw2009 (use the method described in Patton, Politis, and White (2009)).
	- flevel1             <- A function that converts the input dataset to the estimator that the user wishes to bootstrap. Default value is the sample mean.
	- flevel2             <- A function that converts a vector of estimators constructed by flevel1 into a distributional parameter. Default value is sample variance.
	- numobsperresample   <- Number of observations to be drawn (with replacement) per resample. The default value is the number of observations in the dataset (the vast majority of users will want this default value).
	- fblocklengthcombine <- A function for converting a Vector{Float64} of estimated blocklengths to a single Float64 blocklength estimate. Default value is median.

The constructor will attempt to convert all provided keyword arguments to appropriate types,
and will notify the user via an error if a supplied keyword argument is not valid.

Note that the bootmethod and blocklengthmethod keyword arguments will accept both
Symbol and String inputs, and will convert them to BootMethod and BlockLengthMethod
types internally. To see a list of acceptable Symbol or String values for the bootmethod and
blocklengthmethod keyword arguments, use: \n
	- collect(keys(DependentBootstrap.BOOT_METHOD_DICT))
	- collect(keys(DependentBootstrap.BLOCKLENGTH_METHOD_DICT))

respectively. A small proportion of users may need the fine-grained control that comes
from constructing BootMethod and BlockLengthMethod types explicitly and then providing
them to the keyword constructor. These users should use ?BootMethod and ?BlockLengthMethod
at the REPL for more info.

BootInput is not mutable, but the type is near instantaneous to construct, so if a user wishes
to amend a BootInput it is recommended to just construct another one. A special constructor
is provided to facilitate this process that has the method definition: \n
	- BootInput(data, bootinput::BootInput ; kwargs...)

where the new BootInput draws its fields from the keyword arguments that are provided, and then
the input BootInput for any keyword arguments that are not provided.

Note that all exported functions in the DependentBootstrap package exhibit the
method signature: \n
    - exported_func(data ; kwargs...)

which in practice just wraps the keyword argument constructor for a BootInput, and
then calls the method signature: \n
     - exported_func(data, bootinput::BootInput)
"""
struct BootInput{Tbm<:BootMethod,Tbl<:BlockLengthMethod,Tf1<:Function,Tf2<:Function,Tfc<:Function}
	numobs::Int
    blocklength::Float64
    numresample::Int
	bootmethod::Tbm
	blocklengthmethod::Tbl
	flevel1::Tf1
	flevel2::Tf2
    numobsperresample::Int
	fblocklengthcombine::Tfc
	function BootInput(numobs::Int, blocklength::Float64, numresample::Int, bootmethod::Tbm, blmethod::Tbl,
                       flevel1::Tf1, flevel2::Tf2, numobsperresample::Int, fblocklengthcombine::Tfc) where {Tbm<:BootMethod,Tbl<:BlockLengthMethod,Tf1<:Function,Tf2<:Function,Tfc<:Function}
		typeof(bootmethod) <: BootDummy && error("bootmethod is set to BootDummy. You should not have been able to accidentally reach this point. Please file an issue.")
		numobs < 2 && error("Number of observations input to BootInput must be 2 or greater: $(numobs)")
		(isnan(blocklength) || isinf(blocklength)) && error("Invalid blocklength: $(blocklength)")
		numresample < 1 && error("Number of resamples must be strictly positive: $(numresample)")
		numobsperresample < 1 && error("Number of observations per resample must be strictly positive: $(numobsperresample)")
		blocklength <= 0.0 && error("Block length must be strictly positive. BootInput outer constructors that include the dataset as the first argument will automatically estimate the block length if you do not specify it")
		Tbm == BootDummy && error("Do not use BootDummy as an input type. It is for internal module use only")
		Tbl == BLDummy && error("Do not use BLDummy as an input type. It is for internal module use only")
		new{Tbm,Tbl,Tf1,Tf2,Tfc}(numobs, blocklength, numresample, bootmethod, blmethod, flevel1, flevel2, numobsperresample, fblocklengthcombine)
	end
end
#BootInput empty constructor
BootInput() = BootInput(2, 1.0, 1, BootIID(), BLPPW2009(), mean, var, 1, identity)
#BootInput constructor that ensures blocklength is auto-detected if need be
function BootInput(data, numobs::Int, blocklength::Number, numresample::Int, bootmethod::Tbm, blocklengthmethod::Tbl, flevel1::Tf1,
                   flevel2::Tf2, numobsperresample::Int, fblocklengthcombine::Tfc) where {Tbm<:BootMethod,Tbl<:BlockLengthMethod,Tf1<:Function,Tf2<:Function,Tfc<:Function}
	blocklength = Float64(blocklength)
    Tbm <: BootIID && (blocklength = 1.0)
    blocklength <= 0.0 && (blocklength = optblocklength(data, blocklengthmethod, bootmethod, fblocklengthcombine))
	return BootInput(numobs, blocklength, numresample, bootmethod, blocklengthmethod, flevel1, flevel2, numobsperresample, fblocklengthcombine)
end
#BootInput constructors that use keyword arguments
function BootInput(data ; blocklength=0, numresample=NUM_RESAMPLE, bootmethod=:stationary, blocklengthmethod=:dummy,
				   flevel1=mean, flevel2=var, numobsperresample=num_obs(data), fblocklengthcombine=median)
	numobs = num_obs(data)
	blocklength = bootinput_get_blocklength(blocklength)
	numresample = bootinput_get_numresample(numresample)
	bootmethod = bootinput_get_bootmethod(bootmethod)
	blocklengthmethod = bootinput_get_blocklengthmethod(blocklengthmethod, bootmethod)
	flevel1 = bootinput_get_flevel1(flevel1)
	flevel2 = bootinput_get_flevel2(flevel2)
	numobsperresample = bootinput_get_numobsperresample(numobsperresample)
	fblocklengthcombine = bootinput_get_fblocklengthcombine(fblocklengthcombine)
	typeof(bootmethod) <: BootIID && (blocklength = 1.0)
	blocklength <= 0.0 && (blocklength = optblocklength(data, blocklengthmethod, bootmethod, fblocklengthcombine))
	return BootInput(numobs, blocklength, numresample, bootmethod, blocklengthmethod, flevel1, flevel2, numobsperresample, fblocklengthcombine)
end
#Constructor for building a new BootInput using keyword arguments where possible, and, failing that, the fields
#of an existing BootInput
_db_bi_dummy_f() = error("This function is designed to never be called. It is used as a default value for keyword arguments so as to check whether the user has specified them")
function BootInput(data, bi::BootInput ; blocklength=-9, numresample=-9, bootmethod=:z, blocklengthmethod=:z,
				   flevel1=_db_bi_dummy_f, flevel2=_db_bi_dummy_f, numobsperresample=-9,
				   fblocklengthcombine=_db_bi_dummy_f)
    blocklength == -9 && (blocklength = bi.blocklength)
	numresample == -9 && (numresample = bi.numresample)
	bootmethod == :z && (bootmethod = bi.bootmethod)
	blocklengthmethod == :z && (blocklengthmethod = bi.blocklengthmethod)
	flevel1 == _db_bi_dummy_f && (flevel1 = bi.flevel1)
	flevel2 == _db_bi_dummy_f && (flevel2 = bi.flevel2)
	numobsperresample == -9 && (numobsperresample = bi.numobsperresample)
	fblocklengthcombine == _db_bi_dummy_f && (fblocklengthcombine = bi.fblocklengthcombine)
	return BootInput(data, blocklength=blocklength, numresample=numresample, bootmethod=bootmethod, blocklengthmethod=blocklengthmethod,
					 flevel1=flevel1, flevel2=flevel2, numobsperresample=numobsperresample, fblocklengthcombine=fblocklengthcombine)
end
bootinput_get_blocklength(x::Number)::Float64 = Float64(x)
bootinput_get_blocklength(x) = error("Invalid type for blocklength input. Use a subtype of Number, eg Int or Float64.")
bootinput_get_numresample(x::Number)::Int = Int(x)
bootinput_get_numresample(x) = error("Invalid type for numresample input. Use a subtype of Number, eg Int or Float64.")
function bootinput_get_bootmethod(bootmethod::T) where {T<:Union{Symbol,String}}
	bm = get(BOOT_METHOD_DICT, bootmethod, BootDummy())
	typeof(bm) <: BootDummy && error("No matching entry found in dictionary for bootmethod input $(bootmethod). Please use collect(keys(DependentBootstrap.BOOT_METHOD_DICT)) at the REPL to see a list of valid keyword arguments.")
	return bm
end
bootinput_get_bootmethod(x) = typeof(x) <: BootMethod ? x : error("Invalid type for bootmethod input. Use Symbol, String, or a subtype of BootMethod")
function bootinput_get_blocklengthmethod(blocklengthmethod::T, bootmethod::Tbm) where {T<:Union{Symbol,String},Tbm<:BootMethod}
	string(blocklengthmethod) == "dummy" && return bootmethod_to_blocklengthmethod(bootmethod)
	blm = get(BLOCKLENGTH_METHOD_DICT, blocklengthmethod, BLDummy())
	typeof(blm) <: BLDummy && error("No matching entry found in dictionary for blocklengthmethod input $(blocklengthmethod). Please use collect(keys(DependentBootstrap.BLOCKLENGTH_METHOD_DICT)) at the REPL to see a list of valid keyword arguments.")
	return blm
end
bootinput_get_blocklengthmethod(x, bootmethod) = typeof(x) <: BlockLengthMethod ? x : error("Invalid type for bootmethod input. Use Symbol, String, or a subtype of BootMethod")
(bootinput_get_flevel1(f::Tf)::Tf) where {Tf<:Function} = f
bootinput_get_flevel1(f) = error("Invalid type for flevel1 input. Use a subtype of Function")
(bootinput_get_flevel2(f::Tf)::Tf) where {Tf<:Function} = f
bootinput_get_flevel2(f) = error("Invalid type for flevel2 input. Use a subtype of Function")
bootinput_get_numobsperresample(x::Number)::Int = Int(x)
bootinput_get_numobsperresample(x) = error("Invalid type for numobsperresample input. Use a subtype of Number, eg Int or Float64.")
(bootinput_get_fblocklengthcombine(f::Tf)::Tf) where {Tf<:Function} = f
bootinput_get_fblocklengthcombine(f) = error("Invalid type for fblocklengthcombine input. Use a subtype of Function")
function Base.show(io::IO, x::BootInput)
	println(io, "Dependent bootstrap input:")
	println(io, "    Number of observations in dataset = $(x.numobs)")
    println(io, "    Current block length = $(x.blocklength)")
	println(io, "    Number of resamples = $(x.numresample)")
	println(io, "    Bootstrap method = $(x.bootmethod)")
	println(io, "    Block length method = $(x.blocklengthmethod)")
	println(io, "    Level 1 statistic = $(x.flevel1)")
	println(io, "    Level 2 statistic = $(x.flevel2)")
	println(io, "    Number of observations per resample = $(x.numobsperresample)")
	println(io, "    Block length combine function = $(x.fblocklengthcombine)")
end
