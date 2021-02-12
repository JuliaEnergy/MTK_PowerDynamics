# First attempt at VSC in BlockSystems (nee IOSystems)

# Pkg.add(url="https://github.com/hexaeder/IOSystems_prototype")
using BlockSystems
using ModelingToolkit

## Common time

@parameters t 
D = Differential(t)

## Low pass filter

@parameters τ input(t) 
@variables filtered(t)

lpf = IOBlock([D(filtered) ~ 1/τ * (- filtered + input)], [input], [filtered], name=:lp_filter)

## integrator

@parameters x(t)
@variables int(t)
    
integrator = IOBlock([D(int) ~ x], [x], [int], name=:integrator)

## Droop control

@parameters K u_ref x_ref x(t)
@variables u(t)

droop_control = IOBlock([
    u ~ - K * (x - x_ref) + u_ref # output is the droop voltage v
    ], [x], [u], name = :droop)

##

p_filter = IOBlock(lpf, name = :active_power_filter)
q_filter = IOBlock(lpf, name = :reactive_power_filter)
p_droop = IOBlock(droop_control, name = :active_power_droop)
q_droop = IOBlock(droop_control, name = :reactive_power_droop)
f_integrator = IOBlock(integrator, name = :frequency_integrator)

##
@variables ϕ(t) v(t)
@parameters p(t) q(t)

gfi = IOSystem([f_integrator.x => p_droop.u,
          p_droop.x => p_filter.filtered,
          q_droop.x => q_filter.filtered],
          [p_filter, q_filter, p_droop, q_droop, f_integrator],
          name = :GridForming,
          namespace_map = [p_filter.input => p, q_filter.input => q, f_integrator.int => ϕ, q_droop.u => v],
          outputs = [ϕ, v])

connected_gfi = connect_system(gfi)

##

gen = generate_io_function(connected_gfi, f_states=[v, ϕ], f_inputs=[p, q])

##

using OrdinaryDiffEq
pq(t) = [1.0, 1.0]
odefun(du, u, p, t) = gen.f_ip(du, u, pq(t), p, t)

p = [0.5, 1.0, 1.0, 1.0, 0.5, 1.0, 1.0, 1.0]
u0 = [0.0, 0.0, 0.0, 0.0]
tspan = (0.0, 30.0)
prob = ODEProblem(odefun, u0, tspan, p) # Massmatrix missing here
sol = solve(prob, Rodas4())

##

using Plots
plot(sol)

## Infinite grid (WIP)

@parameters Y v(t) ϕ(t)
@variables p(t) q(t)

infinite_grid = IOBlock([
  p ~ real(conj(Y) * v * exp(1im * ϕ) * (v * exp(- 1im * ϕ) - 1)),
  p ~ real(conj(Y) * v * exp(1im * ϕ) * (v * exp(- 1im * ϕ) - 1))], [v, ϕ], [p, q], name = :droop)

##

