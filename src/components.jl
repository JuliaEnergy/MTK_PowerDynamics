using ModelingToolkit


# grid definition
@variables t, u_r(t), u_i(t),ω(t),ϕ(t),p(t)
@parameters i_r, i_i
@derivatives D'~t

slack = [
    D(u_r) ~ -ω*1*sin(ϕ),
    D(u_i) ~  ω*1*cos(ϕ)
    ]
    
aliases = [p ~ u_r*i_r + u_i*i_i]


# outer droop control for virtual intertia one 
@parameters t K_P P_ref  ω_ref 
@variables  u_ϕ(t) p_m(t)
@derivatives D'~t

virtual_intertia_active_power_droop_control = [
    u_ϕ ~ -K_P*(p_m-P_ref) + ω_ref, # output is the droop frequency ω
    ]

# outer droop control for virtual intertia 
@parameters t K_Q Q_ref  V_ref
@variables u_V q_m(t)
@derivatives D'~t

virtual_intertia_reactive_power_droop_control = [
    u_V ~ -K_Q*(q_m-Q_ref) + V_ref # output is the droop voltage v
    ]


    
# component active power filter 
@parameters t τ_P
@variables p(t) p_m(t)
@derivatives D'~t

filter_active_power = [
    D(p_m) ~ 1/τ_P* (-p_m-p),
    ]

# component reactive power filter 
@parameters t τ_Q
@variables q(t) q_m(t)
@derivatives D'~t
    
filter_reactive_power = [
    D(q_m) ~ 1/τ_P* (-q_m-q),
    ]

# component integral controller for voltage 
@parameters t τ_V
@variables v(t) u_V(t)
@derivatives D'~t
    
integral_controller_voltage = [
    D(v) ~ 1/τ_V* (-v+u_V),
    ]

# integrator for frequency 
@parameters t
@variables ϕ(t) ω(t)
@derivatives D'~t
    
integrator_frequency = [
    D(ϕ) ~ ω
    ]