# -*- coding: utf-8 -*-

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
Model equations are as provided by [Shomate1954](@cite).
"""
struct ShomateThermo{K, N} <: ThermodynamicModelData
    data::ThermoData{K, N}
    h_ref::Float64
    s_ref::Float64

    function ShomateThermo(data::ThermoData{K, N}) where {K, N}
        return new{K, N}(data, data.params[end-1, 1], data.params[end, 1])
    end

    function ShomateThermo(data::Vector{Vector{Float64}}, bounds::Vector{Float64})
        return ShomateThermo(ThermoData(data, bounds))
    end
end

# struct MaierKelleyThermo <: ThermodynamicModelData end
# struct EinsteinThermo <: ThermodynamicModelData end

#######################################################################
# INTERNALS
#######################################################################

function factory_symbolic(data::ThermoData{K, N, M}, properties) where {K, N, M}
    @variables T

    jumps = data.bounds[1:N]
    coefs = data.params

    funs = properties(T, SVector{K}(coefs[1:end, 1]))
    fun_cp, fun_hm, fun_sm = funs
    
    for k in range(2, N)
        δ = heaviside(T, jumps[k])

        funs = properties(T, SVector{K}(coefs[1:end, k]))
        new_cp, new_hm, new_sm = funs
        
        Δcp = simplify(new_cp - fun_cp; expand = true)
        Δhm = simplify(new_hm - fun_hm; expand = true)
        Δsm = simplify(new_sm - fun_hm; expand = true)

        fun_cp += δ * Δcp
        fun_hm += δ * Δhm
        fun_sm += δ * Δsm
    end

    # XXX: this is producing identically zero results in some cases
    # for Shomate models of specific heat and entropy. Why? Probably
    # due to some overflow due to the T/1000 factor. Keep this note
    # as this is a weak point of the implementation. It was also
    # observed that this breaks the management of the heaviside
    # intended behavior.
    # fun_cp = simplify(fun_cp; expand = true)
    # fun_hm = simplify(fun_hm; expand = true)
    # fun_sm = simplify(fun_sm; expand = true)

    return fun_cp, fun_hm, fun_sm
end

function factory_symbolic(m::NASAThermo{K, N}) where {K, N}
    return factory_symbolic(m.data, properties_nasa)
end

function factory_symbolic(m::ShomateThermo{K, N}) where {K, N}
    return factory_symbolic(m.data, properties_shomate)
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

const Nasa7Coefs   = SVector{7, Float64}
const Nasa9Coefs   = SVector{9, Float64}
const ShomateCoefs = SVector{8, Float64}

function specific_heat_nasa(T, c::Nasa7Coefs)
    f = c[4] + T * c[5]
    f = c[3] + T * f
    f = c[2] + T * f
    f = c[1] + T * f
    return f
end

function specific_heat_nasa(T, c::Nasa9Coefs)
    error("not implemented")
end

function specific_heat_shomate(t, c::ShomateCoefs)
    p = c[1] + t * (c[2] + t * (c[3] + t * c[4]))
    return p + c[5] / t^2
end

function enthalpy_nasa(T, c::Nasa7Coefs)
    f = c[4] / 4 + T * c[5] / 5
    f = c[3] / 3 + T * f
    f = c[2] / 2 + T * f
    f = c[1] / 1 + T * f
    f = c[6] + T * f
    return f
end

function enthalpy_nasa(T, c::Nasa9Coefs)
    error("not implemented")
end

function enthalpy_shomate(t, c::ShomateCoefs)
    p = t * (c[1] + t * (c[2]/2 + t * (c[3]/3 + t * c[4]/4)))
    return p - c[5] / t + c[6] - c[8]
end

function entropy_nasa(T, c::Nasa7Coefs)
    f = c[4] / 3 + T * c[5] / 4
    f = c[3] / 2 + T * f
    f = c[2] / 1 + T * f
    f = c[7] + c[1] * log(T) + T * f
    return f
end

function entropy_nasa(T, c::Nasa9Coefs)
    error("not implemented")
end

function entropy_shomate(t, c::ShomateCoefs)
    p = log(t) * c[1] + t * (c[2] + t * (c[3]/2 + t * c[4]/3))
    return p - c[5] / (2 * t^2) + c[7]
end

function properties_nasa(T, c)
    eval_cp = GAS_CONSTANT * specific_heat_nasa(T, c)
    eval_hm = GAS_CONSTANT * enthalpy_nasa(T, c)
    eval_sm = GAS_CONSTANT * entropy_nasa(T, c)
    return eval_cp, eval_hm, eval_sm
end

function properties_shomate(T, c)
    eval_cp = specific_heat_shomate(T/1000, c)
    eval_hm = enthalpy_shomate(T/1000, c)
    eval_sm = entropy_shomate(T/1000, c)
    return eval_cp, eval_hm, eval_sm
end

#######################################################################
# FACTORY
#######################################################################

function thermo_models()
    return Dict(
        :NASA7 => NASAThermo,
        :NASA9 => NASAThermo,
        :SHOMATE => ShomateThermo,
        # :MAIERKELLEY => MaierKelleyThermo,
        # :EINSTEIN => EinsteinThermo,
    )
end

function get_thermo_model(name::String)
    name = uppercase(name)
    symb = Symbol(name)
    models = thermo_models()

    if !haskey(models, symb)
        error("""\
        Unknown thermodynamic model `$(name)`; model name must be \
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
        return new(compile_function.(funcs, expression)...)
    end
end

function CompiledThermoFunctions(model::String, data, bounds; 
        how = :symbolic, expression = Val{false})
    funcs = thermo_factory(model, data, bounds; how)
    return new(compile_function.(funcs, expression)...)
end
