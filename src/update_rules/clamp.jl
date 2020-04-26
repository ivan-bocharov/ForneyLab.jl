mutable struct SPClamp <: SumProductRule{Clamp} end
outboundType(::Type{SPClamp}) = Message{PointMass}
isApplicable(::Type{SPClamp{T}}, input_types::Vector{<:Type}) where T<:VariateType = true