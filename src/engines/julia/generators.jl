export algorithmSourceCode

"""
Generate Julia code for message passing
and optional free energy evaluation
"""
function algorithmSourceCode(algo::InferenceAlgorithm; free_energy=false)
    source_code_blocks = Vector{Expr}(undef, length(algo.posterior_factorization.posterior_factors))
    
    for (i, pf) in algo.posterior_factorization
        source_code_blocks[i] = posteriorFactorSourceCode(pf)
    end

    println(source_code_blocks)

    # if free_energy
    #     algo_code *= freeEnergySourceCode(algo)
    #     algo_code *= "\n\n"
    # end

    # algo_code *= "end # block"
end

"""
Generate Julia code for free energy evaluation
"""
function freeEnergySourceCode(algo::InferenceAlgorithm)
    algo.posterior_factorization.free_energy_flag || error("Required quantities for free energy evaluation are not computed by the algorithm. Make sure to flag free_energy=true upon algorithm construction to schedule computation of required quantities.")

    fe_code  = "function freeEnergy$(algo.id)(data::Dict, marginals::Dict)\n\n"
    fe_code *= "F = 0.0\n\n"
    fe_code *= energiesSourceCode(algo.average_energies)
    fe_code *= "\n"
    fe_code *= entropiesSourceCode(algo.entropies)
    fe_code *= "\nreturn F\n\n"
    fe_code *= "end"

    return fe_code
end

"""
Generate Julia code for message passing on a single posterior factor
"""
function posteriorFactorSourceCode(pf::PosteriorFactor)
    pf_code = Vector{Expr}(undef)

    if pf.optimize
        push!(pf_code, optimizeSourceCode(pf))
    end

    if pf.initialize
        push!(pf_code, initializationSourceCode(pf))
    end

    n_entries = length(pf.schedule)
    
    mp_code = quote
        function step$(pf.algorithm_id)$(pf.id)!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=Array{Message}(undef, $n_entries))
            $(scheduleSourceCode(pf.schedule))
            $(marginalTableSourceCode(pf.marginal_table))
            return marginals
        end
    end

    return pf_code
end

"""
Generate template code for optimize block
"""
function optimizeSourceCode(pf::PosteriorFactor)
    # optim_code =  "# You have created an algorithm that requires updates for (a) clamped parameter(s).\n"
    # optim_code *= "# This algorithm requires the definition of a custom `optimize!` function that updates the parameter value(s)\n"
    # optim_code *= "# by altering the `data` dictionary in-place. The custom `optimize!` function may be based on the mockup below:\n\n"
    # optim_code *= "# function optimize$(pf.algorithm_id)$(pf.id)!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=init$(pf.algorithm_id)$(pf.id)())\n"
    # optim_code *= "# \t...\n"
    # optim_code *= "# \treturn data\n"
    # optim_code *= "# end"

    # return optim_code
    return Expr(:block)
end

"""
Generate code for initialization block (if required)
"""
function initializationSourceCode(pf::PosteriorFactor)
    # init_code = "function init$(pf.algorithm_id)$(pf.id)()\n\n"

    # n_messages = length(pf.schedule)

    # init_code *= "messages = Array{Message}(undef, $n_messages)\n\n"

    # for entry in pf.schedule
    #     if entry.initialize
    #         init_code *= "messages[$(entry.schedule_index)] = Message($(vagueSourceCode(entry)))\n"
    #     end
    # end
    # init_code *= "\nreturn messages\n\n"
    # init_code *= "end"

    # return init_code
    return Expr(:block)
end

"""
Generate code for evaluating the average energy
"""
function energiesSourceCode(average_energies::Vector)
    energies_code = ""
    for energy in average_energies
        count_code = countingNumberSourceCode(energy[:counting_number])
        node_code = removePrefix(energy[:node])
        inbounds_code = inboundsSourceCode(energy[:inbounds])
        energies_code *= "F += $(count_code)averageEnergy($node_code, $inbounds_code)\n"
    end

    return energies_code
end

"""
Generate code for evaluating the entropy
"""
function entropiesSourceCode(entropies::Vector)
    entropies_code = ""
    for entropy in entropies
        count_code = countingNumberSourceCode(entropy[:counting_number])
        inbound_code = inboundSourceCode(entropy[:inbound])
        entropies_code *= "F -= $(count_code)differentialEntropy($inbound_code)\n"
    end

    return entropies_code
end

"""
Generate code for counting number
"""
function countingNumberSourceCode(cnt::Int64)
    if cnt == 1
        count_code = ""
    else
        count_code = "$(cnt)*"
    end

    return count_code
end

"""
Construct code for message updates
"""
function scheduleSourceCode(schedule::Schedule)
    schedule_code = Vector{Expr}(undef, length(schedule))
    
    for (i, entry) in enumerate(schedule)
        rule_code = removePrefix(entry.message_update_rule)
        inbounds_code = inboundsSourceCode(entry.inbounds)
        schedule_code[i] = 
            quote
                messages[$(entry.schedule_index)] = rule$(rule_code)($inbounds_code)
            end
    end

    return schedule_code
end

"""
Generate code for marginal updates
"""
function marginalTableSourceCode(table::MarginalTable)
    table_code = Vector{Expr}(undef, length(table))
    
    for (i, entry) in enumerate(table)
        table_code[i] = marginalEntrySourceCode(entry)
    end
    
    return Expr(:block, table_code...)
end

# Should this be dispatched on the rule type?
function marginalEntrySourceCode(entry::MarginalEntry)
    if entry.marginal_update_rule == Nothing
        inbound = entry.inbounds[1]
        return quote
            marginals[:$(entry.marginal_id)] = messages[$(inbound.schedule_index)].dist
        end 
    elseif entry.marginal_update_rule == Product
        inbound1 = entry.inbounds[1]
        inbound2 = entry.inbounds[2]
        return quote
            marginals[:$(entry.marginal_id)] = messages[$(inbound1.schedule_index)].dist * messages[$(inbound2.schedule_index)].dist
        end
    else
        rule_code = removePrefix(entry.marginal_update_rule)
        inbounds_code = inboundsSourceCode(entry.inbounds)
        return quote
            marginals[:$(entry.marginal_id)] =  rule$(rule_code)($inbounds_code)
        end
    end
end

"""
Generate code for vague initializations
"""
function vagueSourceCode(entry::ScheduleEntry)
    family_code = removePrefix(entry.family)
    dims = entry.dimensionality
    if dims == ()
        vague_code = "vague($family_code)"
    else
        vague_code = "vague($family_code, $dims)"
    end

    return vague_code
end

"""
Generate code for message/marginal/energy/entropy computation inbounds
"""
function inboundsSourceCode(inbounds::Vector)
    inbounds_code = String[]
    for inbound in inbounds
        push!(inbounds_code, inboundSourceCode(inbound))
    end
    return join(inbounds_code, ", ")
end

"""
Generate code for a single inbound (overloaded for specific inbound type)
"""
inboundSourceCode(inbound::Nothing) = "nothing"
inboundSourceCode(inbound::ScheduleEntry) = "messages[$(inbound.schedule_index)]"
inboundSourceCode(inbound::MarginalEntry) = "marginals[:$(inbound.marginal_id)]"
function inboundSourceCode(inbound::Dict) # Custom inbound
    keyword_flag = true # Default includes keyword in custom argument
    if haskey(inbound, :keyword)
        keyword_flag = inbound[:keyword]
    end

    inbound_code = ""
    for (key, val) in inbound
        if key != :keyword
            if keyword_flag
                inbound_code = "$(removePrefix(key))=$(removePrefix(val))"
            else
                inbound_code = removePrefix(val)
            end
        end
    end

    return inbound_code
end
function inboundSourceCode(inbound::Clamp{V}) where V<:VariateType # Buffer or value inbound
    dist_or_msg_code = removePrefix(inbound.dist_or_msg)
    variate_type_code = removePrefix(V)

    inbound_code = "$dist_or_msg_code($variate_type_code, PointMass, m="
    if isdefined(inbound, :buffer_id)
        # Inbound is read from buffer
        inbound_code *= "data[:$(inbound.buffer_id)]"
        if isdefined(inbound, :buffer_index) && (inbound.buffer_index > 0)
            inbound_code *= "[$(inbound.buffer_index)]"
        end
    else
        # Inbound is read from clamp node value
        inbound_code *= valueSourceCode(inbound.value)
    end
    inbound_code *= ")"

    return inbound_code
end

inboundSourceCode(inbound::Number) = string(inbound)

"""
Convert a value to parseable Julia code
"""
valueSourceCode(val::Union{Vector, Number}) = string(val)
valueSourceCode(val::Diagonal) = "Diagonal($(val.diag))"
function valueSourceCode(val::AbstractMatrix)
    # If the dimensionality in at least one direction is 1, a matrix needs to be
    # constructed explicitly; otherwise the output Julia code will construct a vector.
    d = size(val)
    if d == (1,1)
        val_code = "mat($(val[1]))" # Shorthand notation for a (1,1) matrix reshape
    elseif 1 in d
        val_code = "reshape($(vec(val)), $d)" # Explicit reshape
    else
        val_code = string(val)
    end

    return val_code
end

"""
Remove module prefixes from types and functions
"""
removePrefix(arg::Any) = arg # Do not remove prefix in general
removePrefix(num::Number) = string(num)
removePrefix(type::Type) = split(string(type), '.')[end]
removePrefix(func::Function) = split(string(func), '.')[end]
