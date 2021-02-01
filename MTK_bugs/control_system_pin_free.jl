## Setup
using ModelingToolkit
const MTK = ModelingToolkit

@parameters t
@derivatives D'~t

## Define an open loop control system:

@variables x(t), u(t), y(t)
@parameters a, b, c, d
ol = ODESystem([D(x) ~ a * x + b * u, y ~ c * x], t, name=:ol)

## Define a proportional control law

@variables u_c(t), y_c(t)
@parameters k_P
pc = ODESystem([u_c ~ k_P * y_c], t,  name=:pc)

## Connect them

connections = [
    ol.u ~ pc.u_c,
    pc.y_c ~ ol.y
]

connected = MTK.ODESystem(connections,t,systems=[ol, pc])

## Get the flattened system

flattened_system = MTK.flatten(connected) # Turns the composite system into one

@show flattened_system.states 
@show flattened_system.pins # should be empty
@show flattened_system.observed
@show flattened_system.ps # a, b, k_P  --> what happened to c?

## Actually reduce the connections

aliased_flattened_system = MTK.alias_elimination(flattened_system) # Connects the pins and the observations, and also removes other equalities.

@show aliased_flattened_system.eqs 
@show aliased_flattened_system.states # only x
@show aliased_flattened_system.pins
@show aliased_flattened_system.observed
@show aliased_flattened_system.ps 