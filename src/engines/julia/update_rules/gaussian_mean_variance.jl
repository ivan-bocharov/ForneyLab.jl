export
ruleSPNormalMVOutNPP,
ruleSPNormalOutNPP
# ruleSPGaussianMeanVarianceMPNP,
# ruleSPGaussianMeanVarianceOutNGP, 
# ruleSPGaussianMeanVarianceMGNP, 
# ruleSPGaussianMeanVarianceVGGN, 
# ruleSPGaussianMeanVarianceVPGN,
# ruleVBGaussianMeanVarianceM,
# ruleVBGaussianMeanVarianceOut

ruleSPNormalMVOutNPP(msg_out::Nothing,
                                    msg_mean::Message{PointMass},
                                    msg_var::Message{PointMass}) =
    Message(NormalMV(deepcopy(msg_mean.dist.m), deepcopy(msg_var.dist.m)))

ruleSPGaussianMeanVarianceMPNP(msg_out::Message{PointMass}, msg_mean::Nothing, msg_var::Message{PointMass}) =
    ruleSPGaussianMeanVarianceOutNPP(msg_mean, msg_out, msg_var)

ruleSPNormalOutNPP(msg_out::Nothing,
    msg_mean::Message{PointMass},
    msg_sigma::Message{PointMass}) =
Message(NormalMV(deepcopy(msg_mean.dist.m), deepcopy(msg_sigma.dist.m)^2))

# function ruleSPGaussianMeanVarianceOutNGP(  msg_out::Nothing,
#                                             msg_mean::Message{F, V},
#                                             msg_var::Message{PointMass}) where {F<:Gaussian, V<:VariateType}

#     d_mean = convert(ProbabilityDistribution{V, GaussianMeanVariance}, msg_mean.dist)

#     Message(V, GaussianMeanVariance, m=d_mean.params[:m], v=d_mean.params[:v] + msg_var.dist.params[:m])
# end

# ruleSPGaussianMeanVarianceMGNP(msg_out::Message{F}, msg_mean::Nothing, msg_var::Message{PointMass}) where F<:Gaussian = 
#     ruleSPGaussianMeanVarianceOutNGP(msg_mean, msg_out, msg_var)

# function ruleSPGaussianMeanVarianceVGGN(msg_out::Message{F1, Univariate},
#                                         msg_mean::Message{F2, Univariate},
#                                         msg_var::Nothing) where {F1<:Gaussian, F2<:Gaussian}

#     d_out  = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_out.dist)
#     d_mean = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_mean.dist)

#     Message(Univariate, Function, log_pdf=(x)-> -0.5*log(d_out.params[:v] + d_mean.params[:v] + x) - 1/(2*x)*(d_out.params[:m] - d_mean.params[:m])^2)
# end

# function ruleSPGaussianMeanVarianceVPGN(msg_out::Message{PointMass, Univariate},
#                                         msg_mean::Message{F, Univariate},
#                                         msg_var::Nothing) where F<:Gaussian

#     d_mean = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance}, msg_mean.dist)

#     Message(Univariate, Function, log_pdf=(x)-> -0.5*log(d_mean.params[:v] + x) - 1/(2*x)*(msg_out.dist.params[:m] - d_mean.params[:m])^2)
# end

# ruleVBGaussianMeanVarianceM(dist_out::ProbabilityDistribution{V},
#                             dist_mean::Any,
#                             dist_var::ProbabilityDistribution) where V<:VariateType =
#     Message(V, GaussianMeanVariance, m=unsafeMean(dist_out), v=unsafeMean(dist_var))

# ruleVBGaussianMeanVarianceOut(  dist_out::Any,
#                                 dist_mean::ProbabilityDistribution{V},
#                                 dist_var::ProbabilityDistribution) where V<:VariateType =
#     Message(V, GaussianMeanVariance, m=unsafeMean(dist_mean), v=unsafeMean(dist_var))