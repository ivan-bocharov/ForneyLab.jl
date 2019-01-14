module PoissonTest

using Test
using ForneyLab
using Herring

import ForneyLab: outboundType, isApplicable, unsafeMean, unsafeVar, slug, isProper, FactorNode, SoftFactor, Interface, FactorGraph
import Herring: VBPoissonOut, VBPoissonL, SPPoissonOutVP, SPPoissonLPV

@testset "Poisson ProbabilityDistribution construction" begin
    @test ProbabilityDistribution(Univariate, Poisson, l=2.0) == ProbabilityDistribution{Univariate, Poisson}(Dict(:l=>2.0))
    @test ProbabilityDistribution(Poisson, l=2.0) == ProbabilityDistribution{Univariate, Poisson}(Dict(:l=>2.0))
    @test ProbabilityDistribution(Poisson) == ProbabilityDistribution{Univariate, Poisson}(Dict(:l=>1.0))
    @test_throws MethodError ProbabilityDistribution(Multivariate, Poisson, l=2.0)
end

@testset "Poisson Message construction" begin
    @test Message(Univariate, Poisson, l=2.0) == Message{Poisson, Univariate}(ProbabilityDistribution{Univariate, Poisson}(Dict(:l=>2.0))    )
    @test Message(Poisson, l=2.0) == Message{Poisson, Univariate}(ProbabilityDistribution{Univariate, Poisson}(Dict(:l=>2.0))    )
    @test Message(Poisson) == Message{Poisson, Univariate}(ProbabilityDistribution{Univariate, Poisson}(Dict(:l=>1.0)))
    @test_throws MethodError Message(Multivariate, Poisson, l=2.0)
end

@testset "dims" begin
    @test dims(ProbabilityDistribution(Poisson)) == 1
end

@testset "slug" begin
    @test slug(Poisson) == "Poisson"
end

@testset "vague" begin
    @test vague(Poisson) == ProbabilityDistribution(Poisson, l=huge)
end

@testset "isProper" begin
    @test isProper(ProbabilityDistribution(Poisson, l=1.0))
    @test !isProper(ProbabilityDistribution(Poisson, l=0.0))
    @test !isProper(ProbabilityDistribution(Poisson, l=-1.0))
end

@testset "unsafe mean and variance" begin
    @test unsafeMean(ProbabilityDistribution(Poisson, l=2.0)) == 2.0
    @test unsafeVar(ProbabilityDistribution(Poisson, l=2.0)) == 2.0
end


#-------------
# Update rules
#-------------

@testset "SPPoissonOutVP" begin
    @test SPPoissonOutVP <: SumProductRule{Poisson}
    @test outboundType(SPPoissonOutVP) == Message{Poisson}
    @test isApplicable(SPPoissonOutVP, [Nothing, Message{PointMass}])

    @test ruleSPPoissonOutVP(nothing, Message(Univariate, PointMass, m=2.0)) == Message(Poisson, l=2.0)
end

@testset "SPPoissonLPV" begin
    @test SPPoissonLPV <: SumProductRule{Poisson}
    @test outboundType(SPPoissonLPV) == Message{Gamma}
    @test isApplicable(SPPoissonLPV, [Message{PointMass}, Nothing])

    @test ruleSPPoissonLPV(Message(Univariate, PointMass, m=2.0), nothing) == Message(Gamma, a=3.0, b=1.0)
end

@testset "VBPoissonOut" begin
    @test VBPoissonOut <: NaiveVariationalRule{Poisson}
    @test outboundType(VBPoissonOut) == Message{Poisson}
    @test isApplicable(VBPoissonOut, [Nothing, ProbabilityDistribution])

    @test ruleVBPoissonOut(nothing, ProbabilityDistribution(Gamma, a=2.0, b=3.0)) == Message(Poisson, l=0.5)
end

@testset "VBPoissonL" begin
    @test VBPoissonL <: NaiveVariationalRule{Poisson}
    @test outboundType(VBPoissonL) == Message{Gamma}
    @test isApplicable(VBPoissonL, [ProbabilityDistribution, Nothing])

    @test ruleVBPoissonL(ProbabilityDistribution(Poisson, l=2.0), nothing) == Message(Gamma, a=3.0, b=1.0)
    @test ruleVBPoissonL(ProbabilityDistribution(Univariate, PointMass, m=2.0), nothing) == Message(Gamma, a=3.0, b=1.0)
end

@testset "averageEnergy and differentialEntropy" begin
    # TODO: Eventually these need to be implemented, but they require additional work for the Poisson distribution
    # @test isapprox(differentialEntropy(ProbabilityDistribution(Poisson, l=2.0)), averageEnergy(Poisson, ProbabilityDistribution(Poisson, l=2.0), ProbabilityDistribution(Univariate, PointMass, m=2.0)))
    @test true == false
end


#------------
# Integration
#------------

@testset "Poisson node construction" begin
    g = FactorGraph()

    test_node = Poisson(Variable(), 1.0)

    # Node should be of correct type
    @test isa(test_node, SoftFactor)

    # Node fields should be of correct types
    @test isa(test_node.id, Symbol)
    @test isa(test_node.interfaces, Vector{Interface})
    @test isa(test_node.i, Dict)

    # Node constructor should automatically assign an id
    @test !isempty(string(test_node.id))

    # Node constructor should assign interfaces to itself
    for iface in test_node.interfaces
        @test ===(iface.node, test_node)
    end

    # Node constructor should add node to graph
    @test ===(g.nodes[test_node.id], test_node)
end

@testset "Parameter estimation" begin
    g = FactorGraph()

    # Construct model
    @RV x ~ Gamma(1.0, 1.0)
    y = Vector{Variable}(undef, 3)
    for k=1:3
        @RV y[k] ~ Poisson(x)
        placeholder(y[k], :y, index=k)
    end

    # Construct algorithm
    algo = sumProductAlgorithm(x)

    # Load algorithm
    eval(Meta.parse(algo))

    # Execute algorithm
    data = Dict(:y => [1.0, 2.0, 3.0])
    marginals = step!(data)

    @test marginals[:x] == ProbabilityDistribution(Gamma, a=7.0, b=4.0)
end

end # module
