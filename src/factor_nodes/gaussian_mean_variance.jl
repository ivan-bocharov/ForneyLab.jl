export NormalMV
using Random:AbstractRNG
"""
Description:

    A Gaussian with mean-variance parameterization:

    f(out,m,v) = 𝒩(out|m,v) = (2π)^{-D/2} |v|^{-1/2} exp(-1/2 (out - m)' v^{-1} (out - m))

Interfaces:

    1. out
    2. m (mean)
    3. v (covariance)


"""
abstract type FLNormal <: ContinuousUnivariateDistribution end

struct NormalMV <: FLNormal
    m::Float64
    V::Float64
    function NormalMV(m::Float64, V::Float64; check_args=false)
        if check_args
            V > 0.0 || error("Variance has to be positive")
        end
        return new(m, V)
    end    
end


mean(dist::NormalMV) = dist.m
var(dist::NormalMV) = dist.V
weightedMean(dist::NormalMV) = cholinv(dist.V)*dist.m
precision(dist::NormalMV) = cholinv(dist.V)

function rand(::AbstractRNG, d::NormalMV)
    return d.m+sqrt(d.V)*randn()
end


function pdf(d::NormalMV, x::Real)
    return exp(-0.5*(d.m-x)^2/(2*d.V))/(sqrt(2*pi*d.V))
end

function logpdf(d::NormalMV, x::Real)
    return log(pdf(d, x))
end 

function cdf(d::NormalMV, x::Real)
    return 0.5*(1+erf(x-d.m/(sqrt(2*d.V))))
end


function quantile(d::NormalMV, q::Real)
    return d.m +sqrt(2*d.V)/erf(2*q-1)
end


minimum(d::NormalMV) = -Inf

maximum(d::UnivariateDistribution) = Inf

insupport(d::UnivariateDistribution, x::Real) = true