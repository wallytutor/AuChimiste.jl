# -*- coding: utf-8 -*-

export NASAThermo
export ShomateThermo
export thermo_factory
export CompiledThermoFunctions

#######################################################################
# GENERAL DATA
#######################################################################

macro thermocoefs(data)
    return quote
        m = reduce(hcat, $(esc(data)))
        SMatrix{size(m)...}(m)
    end
end

"""
Generic storage of thermodynamic data with arbitrary sizes. This structure
is not associated to any specific thermodynamic model/representation.
"""
struct ThermoData{K, N, M}
    "Data matrix with `K` coefficients for `N` temperature ranges."
    params::SMatrix{K, N, Float64}

    "Temperature bounds for `N=M-1` intervals associated to data."
    bounds::NTuple{M, Float64}

    function ThermoData(params::SMatrix{K, N, Float64},
                        bounds::NTuple{M, Float64}) where {K, N, M}
        if M-1 != N
            error("""\
                Bounds size ($M) must be one more than the number of \
                provided coefficient ranges ($N). Please check your data.
                """)
        end

        if !issorted(bounds)
            error("""\
                Tuple of ranges `bounds` must be sorted: got $(bounds).
                """)
        end

        return new{K, N, M}(params, bounds)
    end

    function ThermoData(params::Vector{Vector{Float64}}, bounds::Vector{Float64})
        return ThermoData(@thermocoefs(params), Tuple(bounds))
    end
end

function heaviside(T, r)
    return 0.5 * (sign(T - r) + 1)
 end

#######################################################################
# TYPES
#######################################################################

"""
Stores data for NASA-`K` parametrization with `N` temperature ranges.
"""
struct NASAThermo{K, N} <: ThermodynamicModelData
    data::ThermoData{K, N}
    h_ref::Float64
    s_ref::Float64
    
    function NASAThermo(data::ThermoData{K, N}) where {K, N}
        return new{K, N}(data, data.params[end-1, 1], data.params[end, 1])
    end

    function NASAThermo(data::Vector{Vector{Float64}}, bounds::Vector{Float64})
        return NASAThermo(ThermoData(data, bounds))
    end
end

"""
Stores data for Shomate parametrization with `N` temperature ranges.
"""
struct ShomateThermo{K, N} <: ThermodynamicModelData
end

# struct EinsteinThermo <: ThermodynamicModelData end

#######################################################################
# INTERNALS
#######################################################################

function factory_symbolic(m::NASAThermo{K, N}) where {K, N}
    @variables T

    function thermo_nasa(c)
        eval_cp = specific_heat_nasa(T, c)
        eval_hm = enthalpy_nasa(T, c)
        eval_sm = entropy_nasa(T, c)
        return eval_cp, eval_hm, eval_sm
    end

    jumps = m.data.bounds[1:N]
    coefs = m.data.params

    funs = thermo_nasa(SVector{K}(coefs[1:end, 1]))
    fun_cp, fun_hm, fun_sm = funs
    
    for k in range(2, N)
        δ = heaviside(T, jumps[k])

        funs = thermo_nasa(SVector{K}(coefs[1:end, k]))
        new_cp, new_hm, new_sm = funs
        
        Δcp = new_cp - fun_cp
        Δhm = new_hm - fun_hm
        Δsm = new_sm - fun_hm

        fun_cp += δ * Δcp
        fun_hm += δ * Δhm
        fun_sm += δ * Δsm
    end

    fun_cp = simplify(fun_cp; expand = false)
    fun_hm = simplify(fun_hm; expand = false)
    fun_sm = simplify(fun_sm; expand = false)

    return fun_cp, fun_hm, fun_sm
end

function factory_symbolic(m::ShomateThermo)
    error("not implemented")
end

function factory_numeric(m::NASAThermo{K, N}) where {K, N}
    error("not implemented")
end

function factory_numeric(m::ShomateThermo)
    error("not implemented")
end

#######################################################################
# IMPLEMENTATIONS
#######################################################################

function specific_heat_nasa(T, c::SVector{7, Float64})
    f = c[4] + T * c[5]
    f = c[3] + T * f
    f = c[2] + T * f
    f = c[1] + T * f
    return f
end

function specific_heat_nasa(T, c::SVector{9, Float64})
    error("not implemented")
end

function specific_heat_shomate(T, c::SVector{7, Float64})
    error("not implemented")
end

function enthalpy_nasa(T, c::SVector{7, Float64})
    f = c[4] / 4 + T * c[5] / 5
    f = c[3] / 3 + T * f
    f = c[2] / 2 + T * f
    f = c[1] / 1 + T * f
    return c[6] + T * f
end

function enthalpy_nasa(T, c::SVector{9, Float64})
    error("not implemented")
end

function enthalpy_shomate(T, c::SVector{7, Float64})
    error("not implemented")
end

function entropy_nasa(T, c::SVector{7, Float64})
    f = c[4] / 3 + T * c[5] / 4
    f = c[3] / 2 + T * f
    f = c[2] / 1 + T * f
    return c[7] + c[1] * log(T) + T * f
end

function entropy_nasa(T, c::SVector{9, Float64})
    error("not implemented")
end

function entropy_shomate(T, c::SVector{7, Float64})
    error("not implemented")
end

#######################################################################
# FACTORY
#######################################################################

function thermo_models()
    return Dict(
        :NASA7 => NASAThermo,
        :NASA9 => NASAThermo,
    )
end

function get_thermo_model(name::String)
    symb = Symbol(name)
    models = thermo_models()

    if !haskey(models, symb)
        error("""\
        Unknown thermodynamic model $(name); model name must be \
        among the following: $(keys(models))
        """)
    end
        
    return models[symb]
end

function thermo_factory(m::ThermodynamicModelData; how = :symbolic)
    how == :symbolic && return factory_symbolic(m)
    how == :numeric  && return factory_numeric(m)

    # XXX: for now this seems better as error handling than trying
    # something as getfield(Module, Symbol("nasa7_$(how)"))(m)
    error("""\
        Unknown factory method $(how); currently supported values are \
        `:symbolic` and `:numeric`.
        """)
end

function thermo_factory(model::String, data, bounds; how = :symbolic)
    data_builder = get_thermo_model(model)
    thermo_data = data_builder(data, bounds)
    return thermo_factory(thermo_data; how)
end

function compile_function(f, expression)
    return build_function(f, Symbolics.get_variables(f); expression)
end

struct CompiledThermoFunctions
    specific_heat::Function
    enthalpy::Function
    entropy::Function

    function CompiledThermoFunctions(funcs; expression = Val{false})
        return new(compile_function.(GAS_CONSTANT .* funcs, expression)...)
    end
end
