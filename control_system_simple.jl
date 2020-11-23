## Setup
using ModelingToolkit
const MTK = ModelingToolkit

@parameters t
@derivatives D'~t

## Define an open loop control system:

@variables x(t), u(t), y(t)
@parameters a, b, c, d
ol = ODESystem([D(x) ~ a * x + b * u], t, pins=[u], observed=[y ~ c * x + d * u], name=:ol)

## Define a proportional control law

@variables u_c(t), y_c(t)
@parameters k_P
pc = ODESystem(Array{Equation, 1}(), t, [y_c], [k_P], pins=[y_c], observed=[u_c ~ k_P * y_c],  name=:pc)
# We need to explicitly declare y_c here...

## Connect them
connections = [
    ol.u ~ pc.u_c,
    pc.y_c ~ ol.y
]
# ... or I get the complaint that pc.y_c does not exist here.

connected = MTK.ODESystem(Array{Equation, 1}(),t, observed=connections, systems=[ol, pc])

## Get the flattened system

flattened_system = MTK.flatten(connected) # Turns the composite system into one

@show flattened_system.states 
@show flattened_system.pins # should be empty
@show flattened_system.observed
@show flattened_system.ps # a, b --> what happened to c, k_P ?

## Actually reduce the connections

aliased_flattened_system = MTK.alias_elimination(flattened_system) # Connects the pins and the observations, and also removes other equalities.

# the output does not make sense, the u is set to 0.

@show aliased_flattened_system.eqs 
@show aliased_flattened_system.states # only x
@show aliased_flattened_system.pins # 
@show aliased_flattened_system.observed
@show aliased_flattened_system.ps 