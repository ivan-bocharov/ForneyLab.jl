mutable struct SPClamp <: SumProductRule{Clamp} end
outboundType(::Type{SPClamp}) = Message{PointMass}
isApplicable(::Type{SPClamp}, input_types::Vector{<:Type}) = true