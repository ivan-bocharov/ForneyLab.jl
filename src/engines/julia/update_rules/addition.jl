export  
# ruleSPAdditionOutNGG,
ruleSPAdditionOutNGP,
# ruleSPAdditionOutNPG,
# ruleSPAdditionOutNPP,
# ruleSPAdditionIn1GNG,
ruleSPAdditionIn1PNG,
ruleSPAdditionIn1GNP
# ruleSPAdditionIn1PNP,
# ruleSPAdditionIn2GGN,
# ruleSPAdditionIn2PGN,
# ruleSPAdditionIn2GPN,
# ruleSPAdditionIn2PPN,
# ruleMAdditionNGG

# function ruleSPAdditionOutNGG(
#     msg_out::Nothing,
#     msg_in1::Message{F1, V},
#     msg_in2::Message{F2, V}) where {F1<:Gaussian, F2<:Gaussian, V<:Union{Univariate, Multivariate}}

#     d_in1 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in1.dist)
#     d_in2 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in2.dist)

#     Message(V, GaussianMeanVariance, m=d_in1.params[:m] + d_in2.params[:m], v=d_in1.params[:v] + d_in2.params[:v])
# end

# function ruleSPAdditionIn2GGN(
#     msg_out::Message{F1, V},
#     msg_in1::Message{F2, V},
#     msg_in2::Nothing) where {F1<:Gaussian, F2<:Gaussian, V<:Union{Univariate, Multivariate}}

#     d_in1 = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_in1.dist)
#     d_out = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_out.dist)

#     Message(V, GaussianMeanVariance, m=d_out.params[:m] - d_in1.params[:m], v=d_out.params[:v] + d_in1.params[:v])
# end

# function ruleSPAdditionIn1GNG(
#     msg_out::Message{F1, V},
#     ::Nothing, 
#     msg_in2::Message{F2, V}) where {F1<:Gaussian, F2<:Gaussian, V<:Union{Univariate, Multivariate}}

#     ruleSPAdditionIn2GGN(msg_out, msg_in2, nothing)
# end

function ruleSPAdditionOutNGP(  
    msg_out::Nothing,
    msg_in1::Message{F},
    msg_in2::Message{PointMass}) where {F<:FLNormal}

    d_in1 = convert(NormalMV, msg_in1.dist)

    Message(NormalMV(d_in1.m + msg_in2.dist.m, d_in1.V))
end

# function ruleSPAdditionOutNPG(
#     ::Nothing, 
#     msg_in1::Message{PointMass, V}, 
#     msg_in2::Message{F, V}) where {F<:Gaussian, V<:Union{Univariate, Multivariate}}

#     ruleSPAdditionOutNGP(nothing, msg_in2, msg_in1)
# end

function ruleSPAdditionIn1PNG(  
    msg_out::Message{PointMass},
    msg_in1::Nothing,
    msg_in2::Message{NormalMV})

    Message(NormalMV(msg_out.dist.m - msg_in2.dist.m, msg_in2.dist.V))
end

# function ruleSPAdditionIn2PGN(
#     msg_out::Message{PointMass, V}, 
#     msg_in1::Message{F, V}, 
#     msg_in2::Nothing) where {F<:Gaussian, V<:Union{Univariate, Multivariate}}
    
#     ruleSPAdditionIn1PNG(msg_out, nothing, msg_in1)
# end

function ruleSPAdditionIn1GNP(  
    msg_out::Message{F},
    msg_in1::Nothing,
    msg_in2::Message{PointMass}) where {F<:FLNormal}

    d_out = convert(NormalMV, msg_out.dist)

    Message(NormalMV(d_out.m - msg_in2.dist.m, d_out.V))
end

# function ruleSPAdditionIn2GPN(
#     msg_out::Message{F, V}, 
#     msg_in1::Message{PointMass, V}, 
#     msg_in2::Nothing) where {F<:Gaussian, V<:Union{Univariate, Multivariate}}

#     ruleSPAdditionIn1GNP(msg_out, nothing, msg_in1)
# end

# function ruleSPAdditionOutNPP(
#     msg_out::Nothing, 
#     msg_in1::Message{PointMass, V}, 
#     msg_in2::Message{PointMass, V}) where V<:Union{Univariate, Multivariate}

#     Message(V, PointMass, m=msg_in1.dist.params[:m] + msg_in2.dist.params[:m])
# end

# function ruleSPAdditionIn2PPN(
#     msg_out::Message{PointMass, V}, 
#     msg_in1::Message{PointMass, V}, 
#     msg_in2::Nothing) where V<:Union{Univariate, Multivariate}

#     Message(V, PointMass, m=msg_out.dist.params[:m] - msg_in1.dist.params[:m])
# end

# function ruleSPAdditionIn1PNP(
#     msg_out::Message{PointMass, V}, 
#     msg_in1::Nothing, 
#     msg_in2::Message{PointMass, V}) where V<:Union{Univariate, Multivariate}

#     Message(V, PointMass, m=msg_out.dist.params[:m] - msg_in2.dist.params[:m])
# end

# function ruleMAdditionNGG(
#     msg_out::Message{<:Gaussian, V}, 
#     msg_in1::Message{<:Gaussian, V}, 
#     msg_in2::Message{<:Gaussian, V}) where V<:Union{Univariate, Multivariate}

#     d_out = convert(ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, msg_out.dist)
#     d_in1 = convert(ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, msg_in1.dist)
#     d_in2 = convert(ProbabilityDistribution{V, GaussianWeightedMeanPrecision}, msg_in2.dist)

#     xi_out = d_out.params[:xi]
#     W_out = d_out.params[:w]
#     xi_in1 = d_in1.params[:xi]
#     W_in1 = d_in1.params[:w]
#     xi_in2 = d_in2.params[:xi]
#     W_in2 = d_in2.params[:w]

#     return ProbabilityDistribution(Multivariate, GaussianWeightedMeanPrecision, xi=[xi_in1+xi_out; xi_in2+xi_out], w=[W_in1+W_out W_out; W_out W_in2+W_out])
# end
