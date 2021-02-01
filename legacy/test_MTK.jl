#as in https://github.com/SciML/ModelingToolkit.jl/blob/master/test/reduction.jl

using ModelingToolkit
const MTK = ModelingToolkit
include("src/components.jl")



active_droop = ODESystem(virtual_intertia_active_power_droop_control,t,pins=[p_m],name=:active_droop)#[u_ϕ],[K_P,P_ref,ω_ref],
freq_integrator = ODESystem(integrator_frequency,t,pins=[ω],name=:freq_inegrator)
active_power_filter = ODESystem(filter_active_power,t,pins=[p],name=:active_power_filter)
grid = ODESystem(slack,t,[u_r,u_i,ϕ],[i_r,i_i],pins=[ϕ],observed=aliases,name=:grid)



connections = [active_droop.u_ϕ ~ freq_integrator.ω,
            active_power_filter.p_m ~ active_droop.p_m,
            active_power_filter.p ~ grid.p,
            freq_integrator.ϕ ~ grid.ϕ]

connected = ODESystem(Array{Equation, 1}(), t, [], [], observed = connections, systems=[active_droop,freq_integrator,active_power_filter,grid])

flattened_system = ModelingToolkit.flatten(connected) # Turns the composite system into one
afs = alias_elimination(flattened_system) # Connects the pins and the observations, and also removes other equalities.


# TODO: check why this is not working
#flattened_system = ModelingToolkit.flatten(connected)
#aliased_flattened_system = alias_elimination(connected)
u0 = [active_droop.p_m => 1.0,
      active_droop.u_ϕ => 0.0,
      freq_integrator.ϕ => 0.0,
      freq_integrator.ω => 0.0,
      active_power_filter.p => 1.0,# TODO: this should be elimtaed by flattening and alias alias_elimination
      active_power_filter.p_m => 1.0,
      grid.p => 1.0,
      grid.ϕ => 0.
      ]

p  = [active_droop.P_ref => 1.0,
      active_droop.ω_ref => 0.0,
      active_droop.K_P => 1.,
      a => 0.0]

tspan = (0.0,100.0)
# TODO: build_function bauen, was NetworkDynamics aufruft und VertexFunctions baut
prob = ODEProblem(connected,u0,tspan,p)
sol = solve(prob,Rodas5())

using Plots; plot(sol,vars=(a,lorenz1.x,lorenz2.z))
