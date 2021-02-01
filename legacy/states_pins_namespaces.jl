## Setup
import ModelingToolkit: ODESystem, @parameters, @derivatives, @variables, Differential

@parameters t
@derivatives D'~t

## Define an open loop control system:

@variables x(t), u(t), y(t)
@parameters a, b, c, d
ol_s = ODESystem([D(x) ~ a * x + b * u], t, pins=[u], observed=[y ~ c * x], name=:ol)

@show ol_s.states # [x(t), u(t)]

# So this suggests that pins are states.

## If I explicitly exclude the pins from the states, I don't get namespacing:

ol = ODESystem([D(x) ~ a * x + b * u], t, [x], [a,b], pins=[u], name=:ol)
ol.u # Error: Variable u does not exist

## And then I can not connect to the pin:

@variables y_c(t)
cc = ODESystem([y_c ~ 5], t,  name=:cc)

## Connecting like this:

connections = [
    u ~ cc.y_c,
]

connected = MTK.ODESystem(Equation[],t, observed=connections, systems=[ol, cc])
flattened_system = MTK.flatten(connected)

## Fails to produce the right equations:

flattened_system.eqs[1].rhs # contains the namespaced ol.u(t)
flattened_system.observed[1].lhs # contains the un-namespaced u(t)

## On the other hand, the pins are states version works:

connections_s = [
    ol_s.u ~ cc.y_c,
]

connected_s = MTK.ODESystem(Equation[],t, observed=connections_s, systems=[ol_s, cc])
flattened_system_s = MTK.flatten(connected_s)

@show flattened_system_s.eqs[1].rhs # contains the namespaced ol.u
@show flattened_system_s.observed[1].lhs # contains the namespaced ol.u
