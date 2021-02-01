# Stack Overflow on reducing system:

## Setup
using ModelingToolkit
const MTK = ModelingToolkit

@parameters t
@derivatives D'~t

## Define an open loop control system which feeds the input back into the output:

@variables x(t), u(t), y(t)
@parameters a, b, c, d
ol = ODESystem([D(x) ~ a * x + b * u], t, pins=[u], observed=[y ~ c * x + d * u], name=:ol)

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

## Try to connect them:

flattened_system = MTK.flatten(connected)
aliased_flattened_system = MTK.alias_elimination(flattened_system) # ERROR: StackOverflow
