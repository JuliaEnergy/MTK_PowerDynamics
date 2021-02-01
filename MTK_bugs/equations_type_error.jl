## Setup
using ModelingToolkit
const MTK = ModelingToolkit

@parameters t
@derivatives D'~t

## Define a system

@variables x(t), u(t), y(t)
@parameters a, b, c, d
ol = ODESystem([D(x) ~ a * x + b * u, y ~ c * x], t, name=:ol)

## Define another system that is purely observing

@variables u_c(t), y_c(t)
@parameters k_P
pc = ODESystem(Equation[], t, pins=[y_c], observed = [u_c ~ k_P * y_c], name=:pc)

## Connect them

connections = [
    ol.u ~ pc.u_c,
    y_c ~ ol.y
]

connected = MTK.ODESystem(connections, t, systems=[ol, pc])

##

# equations(connected) has eltype Any:
@show eltype(equations(connected))
@show eltype(parameters(connected)) # So does parameters
@show eltype(parameters(pc)) # So does parameters of pc

# hence flatten doesn't work:
MTK.flatten(connected) # MethodError