## Setup
using ModelingToolkit
const MTK = ModelingToolkit

@parameters t
@derivatives D'~t

## Define an open loop control system:

@variables x(t), u(t), y(t)
@parameters a, b, c, d
ol = ODESystem([D(x) ~ a * x + b * u], t, pins=[u], observed=[y ~ c * x], name=:ol)

## Define a proportional control law

@variables u_c(t), y_c(t)
@parameters k_P
pc = ODESystem(Equation[], t, [y_c], [k_P], pins=[y_c], observed=[u_c ~ k_P * y_c],  name=:pc)

## Connect them
connections = [
    ol.u ~ pc.u_c,
    pc.y_c ~ ol.y
]

connected = MTK.ODESystem(Equation[],t, observed=connections, systems=[ol, pc])

## Get the flattened system

flattened_system = MTK.flatten(connected) # Turns the composite system into one
aliased_flattened_system = MTK.alias_elimination(flattened_system)

##

@show aliased_flattened_system.ps # [ol₊b, ol₊a, pc₊k_P] ol.c is missing

of = ODEFunction(aliased_flattened_system)
of(ones(1), ones(5), 0.) # ERROR: UndefVarError: ol₊c not defined
