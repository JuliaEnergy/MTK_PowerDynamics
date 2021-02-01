using ModelingToolkit
const MTK = ModelingToolkit

struct SystemBase
    D
    t
    function SystemBase()
        MTK.@parameters t
        MTK.@derivatives D'~t
        new(D, t)
    end       
end
  
const standard_base = SystemBase()

function open_loop(;sb::SystemBase=standard_base)
    MTK.@variables x(sb.t), u(sb.t), y(sb.t)
    MTK.@parameters a, b, c, d
    return MTK.ODESystem([sb.D(x) ~ a * x + b * u, y ~ c * x  + d * u,], sb.t, [x, u, y], [a, b, c, d], pins=[u], observed=[ ], name=:open_loop)
end

# Strangely, MTK errors during alias_elimination when I add the input to y as in y ~ c * x + d * u.
# Also, the same point errors when I put the eq in "observed as it might make sense. :-(
    
# ERROR: TypeError: non-boolean (SymbolicUtils.Term{Bool}) used in boolean context
# Stacktrace:
#  [1] (::ModelingToolkit.var"#39#40"{Array{Any,2},Int64})(::Int64) at /home/paul/.julia/packages/ModelingToolkit/hkIWj/src/solve.jl:22
#  [2] iterate at ./generator.jl:47 [inlined]
#  [3] _collect(::UnitRange{Int64}, ::Base.Generator{UnitRange{Int64},ModelingToolkit.var"#39#40"{Array{Any,2},Int64}}, ::Base.EltypeUnknown, ::Base.HasShape{1}) at ./array.jl:699
#  [4] collect_similar at ./array.jl:628 [inlined]
#  [5] map at ./abstractarray.jl:2162 [inlined]
#  [6] sym_lu(::Array{Any,2}) at /home/paul/.julia/packages/ModelingToolkit/hkIWj/src/solve.jl:22
#  [7] _solve(::Array{Expression,2}, ::Array{Operation,1}) at /home/paul/.julia/packages/ModelingToolkit/hkIWj/src/solve.jl:75
#  [8] solve_for(::Array{Equation,1}, ::Array{Operation,1}) at /home/paul/.julia/packages/ModelingToolkit/hkIWj/src/solve.jl:69
#  [9] alias_elimination(::ODESystem) at /home/paul/.julia/packages/ModelingToolkit/hkIWj/src/systems/reduction.jl:71
#  [10] top-level scope at REPL[178]:1

function P_control(;sb::SystemBase=standard_base)
    MTK.@variables u(sb.t), y(sb.t)
    MTK.@parameters k_P
    return MTK.ODESystem(Array{Equation, 1}(), sb.t, [u, y], [k_P, ], pins=[y], observed=[u ~ k_P * y, ],  name=:P_control)
end

ol = open_loop()
pc = P_control()

# this does not work, pins/observed seemingly need to be states since 
# this is producing Equation(u(t), u(t)), Equation(y(t), y(t))
# connections = [
#     pins(ol)[1](standard_base.t) ~  MTK.lhss(MTK.observed(pc))[1], 
#     pins(pc)[1](standard_base.t) ~ MTK.lhss(MTK.observed(ol))[1]
#     ]

connections = [
    ol.u ~ pc.u,
    pc. y ~ ol.y
]

connected = MTK.ODESystem(Array{Equation, 1}(),standard_base.t,[],[], observed=connections, systems=[ol, pc])

flattened_system = MTK.flatten(connected) # Turns the composite system into one

@show flattened_system.states 
@show flattened_system.pins # should be empty
@show flattened_system.observed
@show flattened_system.ps # a, b, k_P  --> what happened to c?

aliased_flattened_system = MTK.alias_elimination(flattened_system) # Connects the pins and the observations, and also removes other equalities.

# the output does not make sense, the u is set to 0.

@show aliased_flattened_system.eqs 
@show aliased_flattened_system.states # only x
@show aliased_flattened_system.pins # 
@show aliased_flattened_system.observed
@show aliased_flattened_system.ps 