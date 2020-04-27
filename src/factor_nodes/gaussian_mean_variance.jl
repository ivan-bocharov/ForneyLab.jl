export NormalMV, NormalWMP, NormalMP, prod!, convert

"""
Description:

    A Gaussian with mean-variance parameterization:

    f(out,m,v) = 𝒩(out|m,v) = (2π)^{-D/2} |v|^{-1/2} exp(-1/2 (out - m)' v^{-1} (out - m))

Interfaces:

    1. out
    2. m (mean)
    3. v (covariance)

Construction:

    GaussianMeanVariance(out, m, v, id=:some_id)
"""
abstract type FLNormal <: ContinuousUnivariateDistribution end

struct NormalMV <: FLNormal
    m::Float64
    V::Float64
end

struct NormalWMP <: FLNormal
    xi::Float64
    w::Float64
end

struct NormalMP <: FLNormal
    m::Float64
    T::Float64
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



# mutable struct GaussianMeanVariance <: Gaussian
#     id::Symbol
#     interfaces::Vector{Interface}
#     i::Dict{Symbol,Interface}

#     function GaussianMeanVariance(out, m, v; id=generateId(GaussianMeanVariance))
#         @ensureVariables(out, m, v)
#         self = new(id, Array{Interface}(undef, 3), Dict{Symbol,Interface}())
#         addNode!(currentGraph(), self)
#         self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
#         self.i[:m] = self.interfaces[2] = associate!(Interface(self), m)
#         self.i[:v] = self.interfaces[3] = associate!(Interface(self), v)

#         return self
#     end
# end

# slug(::Type{GaussianMeanVariance}) = "𝒩"

# format(dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType = "$(slug(GaussianMeanVariance))(m=$(format(dist.params[:m])), v=$(format(dist.params[:v])))"

# ProbabilityDistribution(::Type{Univariate}, ::Type{GaussianMeanVariance}; m=0.0, v=1.0) = ProbabilityDistribution{Univariate, GaussianMeanVariance}(Dict(:m=>m, :v=>v))
# ProbabilityDistribution(::Type{GaussianMeanVariance}; m::Number=0.0, v::Number=1.0) = ProbabilityDistribution{Univariate, GaussianMeanVariance}(Dict(:m=>m, :v=>v))
# ProbabilityDistribution(::Type{Multivariate}, ::Type{GaussianMeanVariance}; m=[0.0], v=mat(1.0)) = ProbabilityDistribution{Multivariate, GaussianMeanVariance}(Dict(:m=>m, :v=>v))

# dims(dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType = length(dist.params[:m])

# vague(::Type{GaussianMeanVariance}) = ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=huge)
# vague(::Type{GaussianMeanVariance}, dims::Int64) = ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=zeros(dims), v=huge*diageye(dims))
# vague(::Type{GaussianMeanVariance}, dims::Tuple{Int64}) = ProbabilityDistribution(Multivariate, GaussianMeanVariance, m=zeros(dims), v=huge*diageye(dims[1]))

# unsafeMean(dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType = deepcopy(dist.params[:m]) # unsafe mean

# unsafeMode(dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType = deepcopy(dist.params[:m]) # unsafe mode

# unsafeVar(dist::ProbabilityDistribution{Univariate, GaussianMeanVariance}) = dist.params[:v] # unsafe variance
# unsafeVar(dist::ProbabilityDistribution{Multivariate, GaussianMeanVariance}) = diag(dist.params[:v])

# unsafeCov(dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType = deepcopy(dist.params[:v]) # unsafe covariance

# unsafeMeanCov(dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType = (deepcopy(dist.params[:m]), deepcopy(dist.params[:v]))

weightedMean(dist::NormalMV) = cholinv(dist.V)*dist.m

precision(dist::NormalMV) = cholinv(dist.V)

# logPdf(dist::ProbabilityDistribution{Univariate, GaussianMeanVariance},x) = -0.5*(log(2pi)+log(dist.params[:v]) + (x-dist.params[:m])^2/dist.params[:v])
# logPdf(dist::ProbabilityDistribution{Multivariate, GaussianMeanVariance},x) = -0.5*(dims(dist)*log(2pi) + log(det(dist.params[:v])) + transpose(x-dist.params[:m])*inv(dist.params[:v])*(x-dist.params[:m]))

# @symmetrical function prod!(
#     x::ProbabilityDistribution{Univariate, Function},
#     y::ProbabilityDistribution{Univariate, F},
#     z::ProbabilityDistribution{Univariate, GaussianMeanVariance}=ProbabilityDistribution(Univariate, GaussianMeanVariance, m=0.0, v=1.0)) where {F<:Gaussian}

#     # The product of a log-pdf and Gaussian distribution is computed by importance sampling
#     y = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, y)
#     samples = y.params[:m] .+ sqrt(y.params[:v]).*randn(1000)

#     p = exp.((x.params[:log_pdf]).(samples))
#     Z = sum(p)
#     mean = sum(p./Z.*samples)
#     var = sum(p./Z.*(samples .- mean).^2)

#     z.params[:m] = mean
#     z.params[:v] = var
#     return z
# end

# function prod!(
#     x::ProbabilityDistribution{Univariate, Function},
#     y::ProbabilityDistribution{Univariate, Function},
#     z::ProbabilityDistribution{Univariate, Function}=ProbabilityDistribution(Univariate, Function, log_pdf=(s)->s))

#     z.params[:log_pdf] = ((s) -> x.params[:log_pdf](s) + y.params[:log_pdf](s))

#     return z
# end

# # Converting from m,v to xi,w would require two separate inversions of the covariance matrix;
# # this function ensures only a single inversion is performed
# function unsafeWeightedMeanPrecision(dist::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType
#     W = cholinv(dist.params[:v])
#     return (W*dist.params[:m], W)
# end

# isProper(dist::ProbabilityDistribution{Univariate, GaussianMeanVariance}) = (floatmin(Float64) < dist.params[:v] < floatmax(Float64))
# isProper(dist::ProbabilityDistribution{Multivariate, GaussianMeanVariance}) = isRoundedPosDef(dist.params[:v])

# function ==(t::ProbabilityDistribution{V, GaussianMeanVariance}, u::ProbabilityDistribution{V, GaussianMeanVariance}) where V<:VariateType
#     (t === u) && return true
#     isApproxEqual(t.params[:m], u.params[:m]) && isApproxEqual(t.params[:v], u.params[:v])
# end

# function ==(t::ProbabilityDistribution{V, GaussianMeanVariance}, u::ProbabilityDistribution{V, F}) where {V<:VariateType, F<:Gaussian}
#     (t === u) && return true
#     isApproxEqual(t.params[:m], unsafeMean(u)) && isApproxEqual(t.params[:v], unsafeCov(u))
# end

# # Average energy functional
# function averageEnergy(::Type{GaussianMeanVariance}, marg_out::ProbabilityDistribution{Univariate}, marg_mean::ProbabilityDistribution{Univariate}, marg_var::ProbabilityDistribution{Univariate})
#     (m_mean, v_mean) = unsafeMeanCov(marg_mean)
#     (m_out, v_out) = unsafeMeanCov(marg_out)

#     0.5*log(2*pi) +
#     0.5*unsafeLogMean(marg_var) +
#     0.5*unsafeInverseMean(marg_var)*(v_out + v_mean + (m_out - m_mean)^2)
# end

# function averageEnergy(::Type{GaussianMeanVariance}, marg_out::ProbabilityDistribution{Multivariate}, marg_mean::ProbabilityDistribution{Multivariate}, marg_var::ProbabilityDistribution{MatrixVariate})
#     (m_mean, v_mean) = unsafeMeanCov(marg_mean)
#     (m_out, v_out) = unsafeMeanCov(marg_out)

#     0.5*dims(marg_out)*log(2*pi) +
#     0.5*unsafeDetLogMean(marg_var) +
#     0.5*tr( unsafeInverseMean(marg_var)*(v_out + v_mean + (m_out - m_mean)*(m_out - m_mean)'))
# end

# function averageEnergy(::Type{GaussianMeanVariance}, marg_out_mean::ProbabilityDistribution{Multivariate, F}, marg_var::ProbabilityDistribution{Univariate}) where F<:Gaussian
#     (m, V) = unsafeMeanCov(marg_out_mean)

#     0.5*log(2*pi) +
#     0.5*unsafeLogMean(marg_var) +
#     0.5*unsafeInverseMean(marg_var)*( V[1,1] - V[1,2] - V[2,1] + V[2,2] + (m[1] - m[2])^2 )
# end

# function averageEnergy(::Type{GaussianMeanVariance}, marg_out_mean::ProbabilityDistribution{Multivariate, F}, marg_var::ProbabilityDistribution{MatrixVariate}) where F<:Gaussian
#     (m, V) = unsafeMeanCov(marg_out_mean)
#     d = Int64(dims(marg_out_mean)/2)

#     0.5*d*log(2*pi) +
#     0.5*unsafeDetLogMean(marg_var) +
#     0.5*tr( unsafeInverseMean(marg_var)*( V[1:d,1:d] - V[1:d,d+1:end] - V[d+1:end,1:d] + V[d+1:end,d+1:end] + (m[1:d] - m[d+1:end])*(m[1:d] - m[d+1:end])' ) )
# end
