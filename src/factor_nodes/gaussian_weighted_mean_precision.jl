export NormalWMP, prod!, convert

"""
Description:

    A Gaussian with weighted-mean-precision parameterization:

    f(out,xi,w) = 𝒩(out|xi,w) = (2π)^{-D/2} |w|^{1/2} exp(-1/2 xi' w^{-1} xi) exp(-1/2 xi' w xi + out' xi)

Interfaces:

    1. out
    2. xi (weighted mean, w*m)
    3. w (precision)

"""

struct NormalWMP <: FLNormal
    xi::Float64
    w::Float64
end

function convert(::Type{NormalWMP}, dist::NormalMV)
    w = cholinv(dist.V)
    xi = w*dist.m

    return NormalWMP(xi, w)
end

function convert(::Type{NormalMV}, dist::NormalWMP)
    v = cholinv(dist.w)
    m = v*dist.xi

    return NormalMV(m, v)
end

weightedMean(x::NormalWMP) = x.xi
precision(x::NormalWMP) = x.w

mean(x::NormalWMP) = cholinv(x.w)*x.xi
var(x::NormalWMP) =  1.0/x.w

function prod!(
    x::D1,
    y::D2,
    z::NormalWMP=NormalWMP(0.0, 1.0)) where {D1<:FLNormal, D2<:FLNormal}

    return NormalWMP(weightedMean(x) + weightedMean(y), precision(x) + precision(y))
end

# mutable struct GaussianWeightedMeanPrecision <: Gaussian
#     id::Symbol
#     interfaces::Vector{Interface}
#     i::Dict{Symbol,Interface}

#     function GaussianWeightedMeanPrecision(out, xi, w; id=generateId(GaussianWeightedMeanPrecision))
#         @ensureVariables(out, xi, w)
#         self = new(id, Array{Interface}(undef, 3), Dict{Symbol,Interface}())
#         addNode!(currentGraph(), self)
#         self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
#         self.i[:xi] = self.interfaces[2] = associate!(Interface(self), xi)
#         self.i[:w] = self.interfaces[3] = associate!(Interface(self), w)

#         return self
#     end
# end

# slug(::Type{GaussianWeightedMeanPrecision}) = "𝒩"

# format(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = "$(slug(GaussianWeightedMeanPrecision))(xi=$(format(dist.params[:xi])), w=$(format(dist.params[:w])))"

# ProbabilityDistribution(::Type{Univariate}, ::Type{GaussianWeightedMeanPrecision}; xi=0.0, w=1.0) = ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
# ProbabilityDistribution(::Type{GaussianWeightedMeanPrecision}; xi::Number=0.0, w::Number=1.0) = ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))
# ProbabilityDistribution(::Type{Multivariate}, ::Type{GaussianWeightedMeanPrecision}; xi=[0.0], w=transpose([1.0])) = ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}(Dict(:xi=>xi, :w=>w))

# dims(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = length(dist.params[:xi])

# vague(::Type{GaussianWeightedMeanPrecision}) = ProbabilityDistribution(Univariate, GaussianWeightedMeanPrecision, xi=0.0, w=tiny)
# vague(::Type{GaussianWeightedMeanPrecision}, dims::Int64) = ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(dims), w=tiny*diageye(dims))
# vague(::Type{GaussianWeightedMeanPrecision}, dims::Tuple{Int64}) = ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=zeros(dims), w=tiny*diageye(dims[1]))

# unsafeMean(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = cholinv(dist.params[:w])*dist.params[:xi]

# unsafeMode(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = cholinv(dist.params[:w])*dist.params[:xi]

# unsafeVar(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}) = 1.0/dist.params[:w] # unsafe variance
# unsafeVar(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}) = diag(cholinv(dist.params[:w]))

# unsafeCov(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = cholinv(dist.params[:w])

# function unsafeMeanCov(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
#     v = cholinv(dist.params[:w])
#     return (v*dist.params[:xi], v)
# end

# unsafeWeightedMean(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = deepcopy(dist.params[:xi]) # unsafe weighted mean

# unsafePrecision(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = deepcopy(dist.params[:w]) # unsafe precision

# unsafeWeightedMeanPrecision(dist::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType = (deepcopy(dist.params[:xi]), deepcopy(dist.params[:w]))

# isProper(dist::ProbabilityDistribution{Univariate, GaussianWeightedMeanPrecision}) = (floatmin(Float64) < dist.params[:w] < floatmax(Float64))
# isProper(dist::ProbabilityDistribution{Multivariate, GaussianWeightedMeanPrecision}) = isRoundedPosDef(dist.params[:w])

# function ==(t::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, u::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}) where V<:VariateType
#     (t === u) && return true
#     isApproxEqual(t.params[:xi], u.params[:xi]) && isApproxEqual(t.params[:w], u.params[:w])
# end

# function ==(t::ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, u::ProbabilityDistribution{V, F}) where {V<:VariateType, F<:Gaussian}
#     (t === u) && return true
#     isApproxEqual(t.params[:xi], unsafeWeightedMean(u)) && isApproxEqual(t.params[:w], unsafePrecision(u))
# end