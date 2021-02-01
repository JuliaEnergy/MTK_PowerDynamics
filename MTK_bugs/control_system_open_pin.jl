# This example shows open pins not being propagated by flatten system

## Setup
using ModelingToolkit
const MTK = ModelingToolkit

@parameters t
@derivatives D'~t

## Define an open loop control system:

@variables x(t), u(t), u1(t), y(t)
@parameters a, b, c, d
ol = ODESystem([D(x) ~ a * x + b * u - u1], t, pins=[u, u1], observed=[y ~ c * x], name=:ol)

# We will keep u1 open

## Define a proportional control law

@variables u_c(t), y_c(t)
@parameters k_P
pc = ODESystem(Equation[], t, [y_c], [k_P], pins=[y_c], observed=[u_c ~ k_P * y_c],  name=:pc)
# We need to explicitly declare y_c here...

## Connect them
connections = [
    ol.u ~ pc.u_c,
    pc.y_c ~ ol.y
]
# ... or I get the complaint that pc.y_c does not exist here. This is an example
# of not being able to connect to pins unless they are states
# (they are not namespaced unless they are states either it seems)

connected = MTK.ODESystem(Equation[],t, observed=connections, systems=[ol, pc])

## Get the flattened system

flattened_system = MTK.flatten(connected) # Turns the composite system into one

@show flattened_system.eqs # same as ol.eqs, the connections do not occur here...
@show flattened_system.pins # Num[] should be ol₊u1(t)
@show flattened_system.observed # Num[] should be ol₊u1(t)

## Actually reduce the connections

aliased_flattened_system = MTK.alias_elimination(flattened_system) # ERROR: StackOverflow

@show aliased_flattened_system.eqs # same as ol.eqs, the connections do not occur here...
@show aliased_flattened_system.pins # Num[] should be ol₊u1(t)
@show aliased_flattened_system.ps # Num[] should be ol₊u1(t)

