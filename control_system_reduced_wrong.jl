# This example shows all equations reduced incorrectly by alias removal:

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
pc = ODESystem(Array{Equation, 1}(), t, [y_c], [k_P], pins=[y_c], observed=[u_c ~ k_P * y_c],  name=:pc)
# We need to explicitly declare y_c here...

## Connect them
connections = [
    ol.u ~ pc.u_c,
    pc.y_c ~ ol.y
]
# ... or I get the complaint that pc.y_c does not exist here. This is an example
# of not being able to connect to pins unless they are states
# (they are not namespaced unless they are states either it seems)

connected = MTK.ODESystem(Array{Equation, 1}(),t, observed=connections, systems=[ol, pc])

## Get the flattened system

flattened_system = MTK.flatten(connected) # Turns the composite system into one

@show flattened_system.eqs # same as ol.eqs, would expect the connections to occur here...
@show flattened_system.states # [ol₊x(t), ol₊u(t)]
@show flattened_system.pins # Num[]
@show flattened_system.observed # all the connections occur here...
@show flattened_system.ps # [ol₊b, ol₊a] only the parameters in the eqs not those in the connections show up here,
# these are not sufficient to run the connected system.

## Actually reduce the connections

aliased_flattened_system = MTK.alias_elimination(flattened_system)

@show aliased_flattened_system.eqs # The equation is incorrect, an Infinity occurs here
@show aliased_flattened_system.states # [ol₊x(t)]
@show aliased_flattened_system.pins # Num[]
@show aliased_flattened_system.observed # Lots of inifinities in these equations as well...

##

@show aliased_flattened_system.ps # [ol₊b, ol₊a] and now this leads to problems:

of = ODEFunction(aliased_flattened_system)
of(ones(1), ones(5), 0.) # ERROR: UndefVarError: ol₊c not defined

##

of = ODEFunction(ODESystem(aliased_flattened_system.eqs))
of(ones(1), ones(5), 0.) # Inf.