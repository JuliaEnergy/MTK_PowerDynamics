## Setup

# A general note for handling equations and terms: == is overloaded and does not return bool,
# Many of the normal ways of handling arrays don't work. Use isequal instead. Set based manipulations
# seem fine. More generally you can check here https://docs.julialang.org/en/v1/base/collections/ which functions use isequal instead of ==

using ModelingToolkit
const MTK = ModelingToolkit

##

# Right now we don't define IOSystem as a subtype of MTK.AbstractSystem as we don't
# define all neccessary fields (leading to StackOverflow on field access).

"""
IOSystems model systems of the form:

dx/dt = f(x,i)

o = g(x)

"""
mutable struct IOSystem # variables are all that occur in eqs, inputs and outputs are subsets of the variables.
    eqs::Vector{Equation}
    variables
    inputs
    outputs
    function IOSystem(eqs; inputs, outputs)
        os = ODESystem(eqs, name=:sys) # We use ODESystem to analyze the equations for us.
        # We could also construct from two sets of eqs: dynamics (f) and output (g).
        
        state_vars = setdiff(os.states, inputs ∪ outputs)
        # Make sure the outputs and inputs actually occur in the equations
        @assert Set(inputs) ⊆ Set(os.states) "Inputs need to occur in the equations"
        @assert Set(outputs) ⊆ Set(os.states) "Outputs need to occur in the equations"
        @assert isempty(Set(inputs) ∩ Set(outputs)) "Outputs and Inputs need to be disjoint"
        @assert isequal( Set(inputs) ∪ Set(outputs) ∪ Set(state_vars), Set(os.states) ) "This assert should be redundant and never be reached"
        

        # We require output variables to occur as left hand side in the equations
        lhs_eqs = [eq.lhs for eq in eqs]
        rhs_eqs = [eq.rhs for eq in eqs]
        @assert Set(outputs) ⊆ Set(lhs_eqs) "Outputs need to be the left hand side of an equation."
        
        # For now we neglect the constraints we talk about here....

        # for eq in eqs
        #     if Set([eq.lhs]) ⊆ Set(outputs)
        #         @assert isempty( MTK.vars(eq.rhs) ∩ Set(inputs) ) "Inputs may not appear in output equations"
        #     end
        # end

        # We make sure the right hand side of output equations does not contain inputs.
        # Otherwise we can not eliminate inputs simply by plugging in the outputs, but need to potentially solve equations.
        # One can always work around this by introducing an intermediate variable. That will show up as a constraint equation
        # in the dynamics.
        # There might be additional constraints we need for our IO System equations, like we don't want any outputs to have two equations.

        new(eqs, state_vars, inputs, outputs)
    end
end

##

@parameters t
@derivatives D'~t

@variables x(t), u(t), y(t)
@parameters a, b, c, d

##

IOSystem([D(x) ~ a * x + b * u, y ~ c * x], inputs = [u], outputs = [y])

##

# A note on types in MTK, explicit type signatures are tricky because there are a lot of types floating around:

os = ODESystem([D(x) ~ a * x + b * u, y ~ c * x], name=:sys1)
@show typeof(x) # typeof(x) = Num
@show typeof(x.val) # typeof(x.val) = Term{Real}
@show typeof(os.states[1]) # typeof(os.states[1]) = Term{Real}

@show isequal(x, os.states[1]) # true
@show Num <: Term{Real} # false
@show Term{Real} <: Num # false

##


# A note on namespacing in MTK:
os = ODESystem([D(x) ~ a * x + b * u, y ~ c * x + d * u], name=:sys1)

@show os.y # sys1₊y(t)
@show os.eqs[2].lhs # y(t)
@show os.states[2] # y(t)

@variables z(t)
eqs2 = vcat(os.eqs, [D(z) ~ os.y]) # This looks natural for adding an equation to a system but is incorrect:
os2 = ODESystem(eqs2)
@show os2.states[[2,5]] # [y(t), sys1₊y(t)] so the namespaced access with os.y and the variable in the equation end up being different.

# Accessing variables through the dot notation is namespaced, and combining them namespaces them, too.
# There are a number of undocumented helperfunctions (like namespace_equations that gets namespace equations)
# that do what you would expect them to.

eqs3 = MTK.namespace_equations(os)
os3 = ODESystem(vcat(eqs3, [D(z) ~ os.y]), name=:newsys)
@show os3.states # this is a namespaced 4 equation array, but now...
try 
    os3.y 
catch err @assert err isa ErrorException end # ... this is an error, and this:
os3.sys1₊y # newsys₊sys1₊y(t)
# is awkward and error prone.

# We might want to explore this space a little. As there might be many u and y variables flying around, namespacing is important.
# We probably don't want it to be automatic. I am leaning towards adding optional namespacing to connect, and checking
# for duplicate variable names when connecting systems.

##

function connect(connections, io_systems::Array{IOSystem,1})
    # The connections are pairs o => i that tell us which output to connect to which input.

    # For now simply concatenate everything. TODO: make sure things are suitably disjoint!
    all_inputs = vcat((ios.inputs for ios in io_systems)...)
    all_outputs = vcat((ios.outputs for ios in io_systems)...)

    @assert allunique(vcat(all_inputs, all_outputs)) "All outputs and all inputs in the systems to connect need to be disjoint."

    output_set = Set(all_outputs)
    input_set = Set(all_inputs)

    connected_inputs = [conn.second for conn in connections]
    connected_outputs = [conn.first for conn in connections]
    
    @assert allunique(connected_inputs) "Can not connect two or more connections to same input"

    @assert Set(connected_inputs) ⊆ Set(all_inputs) "Must connect to inputs"
    @assert Set(connected_outputs) ⊆ Set(all_outputs) "Must connect from outputs"
    
    eqs = vcat((ios.eqs for ios in io_systems)...)

    # as == is overloaded for variables one has to use isequal(x, y) to compare Terms instead.
    # This is the core substitution: We substitute the input that we are connecting to (conn.second)
    # with the right hand side of the output (conn.first) we're plugging into the input.
    sub_rules = [conn.second => eqs[findfirst(x -> isequal(x.lhs, conn.first), eqs)].rhs  for conn in connections]
    
    # Because we do not enforce the "no inputs in the outputs" equations rule, this can cause a stack overflow.
    # If each substitution introduces a new term that needs to be substituted we get a runaway recursion.
    new_eqs = [eq.lhs ~ substitute(eq.rhs, sub_rules) for eq in eqs]
    
    open_inputs = setdiff(Set(all_inputs), Set(connected_inputs))

    IOSystem(new_eqs; inputs = open_inputs, outputs = all_outputs)
end

## Define an open loop control system:

@variables x(t), u(t), y(t)
@parameters a, b, c, d
ol = IOSystem([D(x) ~ a * x + b * u, y ~ c * x]; inputs = [u], outputs = [y])

## Define a proportional control law

@variables u_c(t), y_c(t)
@parameters k_P
# With the "no inputs in the outputs rule we could not write simple replacements like this:
pc = IOSystem([u_c ~ k_P * y_c]; inputs = [y_c], outputs = [u_c])

## Connect them:

closed_loop = connect([u_c => u, y => y_c], [ol, pc])


##

function build_io_functions(io_sys::IOSystem; first_outputs = [], first_inputs = [], first_dynamic_states = []) # The outputs, inputs and dynamic_states parameters are merely to enforce ordering.
    outputs = vcat(first_outputs, io_sys.outputs) |> unique # We add the variables we want first to the start of the array and then we make unique, which is order preserving.
    inputs = vcat(first_inputs, io_sys.inputs) |> unique

    @assert Set(first_dynamic_states) ⊆ Set(io_sys.variables) "dynamic_states parameter must contain variables that occur in the equations"
    @assert Set(io_sys.outputs) == Set(outputs) && length(io_sys.outputs) == length(outputs) "parameter outputs needs to contain ouptuts of the IO System"
    @assert Set(io_sys.inputs) == Set(inputs) && length(io_sys.inputs) == length(inputs) "parameter inputs needs to contain inputs of the IO System"
    
    out_eqs_idx = [findfirst(x -> isequal(x.lhs, o), io_sys.eqs) for o in outputs]

    output_equations = io_sys.eqs[out_eqs_idx]
    dynamic_equations_unordered = setdiff(io_sys.eqs, output_equations)

    os = ODESystem(dynamic_equations_unordered)
    
    not_outputs = setdiff(os.states, io_sys.outputs)
    dynamic_states_unordered = setdiff(not_outputs, io_sys.inputs)

    dynamic_states = vcat(first_dynamic_states, dynamic_states_unordered) |> unique # We add the ordered states to the start of the array and then we make unique, which is order preserving.
    
    # Then we make sure we didn't add any new variables accidentally here:
    @assert length(dynamic_states) == length(dynamic_states_unordered) "dynamic_states parameter must contain variables that occur in the dynamic equations"

    dyn_eqs_idx = [findfirst(x -> Set([s]) ⊆ MTK.vars(x.lhs), dynamic_equations_unordered) for s in dynamic_states]
    
    dynamic_equations = dynamic_equations_unordered[dyn_eqs_idx]

    dynamic_formulas = [eq.rhs for eq in dynamic_equations] # We assume that the left hand side is d/dx. Otherwise we need to add mass matrix support and use the make_lhs_0 function to get the right formulas.
    @assert all( [ Set([state]) ⊆ MTK.vars(eq.lhs) for (state, eq) in zip(dynamic_states, dynamic_equations)] )

    f_oop, f_ip = build_function(dynamic_formulas, dynamic_states, io_sys.inputs, os.ps, os.iv; expression = Val{false})
    # g_oop, g_ip = build_functions(output_equations, dynamic_states, , os.iv)
    f_oop, f_ip, dynamic_states
end

# We need an IOFunctions object that contains a compiled f, a compiled g, mass matrix for f, symbols and jacobians.

#=
struct IOFunction
    f_ip
    f_oop
    g_ip
    g_oop
    mass_matrix
    dynamic_symbols
    input_symbols
    output_symbols
    parameter_symbols
    f_jac
    g_jac
end
=#

##

f_oop, f_ip, ds = build_io_functions(ol)

f_oop([1], [2], [3,4], 0.)

##

using NetworkDynamics

struct NetworkLayer
    vertex_vars
    edge_vars 
    vertex_s_vars # We need these separately so we can have two inputs of the edge functions
    vertex_d_vars # 
    aggregate_vars
    aggregate # The function that aggregates
end


@variables ϕ(t), ϕ_s(t), ϕ_d(t), f(t), total_flow(t)

@inline function flow_sum_1(e_s, e_d)
    e_sum = 0.
    @inbounds for e in e_s
        e_sum -= e[1]
    end
    @inbounds for e in e_d
        e_sum += e[1]
    end
    e_sum
end

KuramotoNetwork = NetworkLayer([ϕ], [f], [ϕ_s], [ϕ_d], [total_flow], flow_sum_1)

##

function build_nd_vertex(ios::IOSystem, nl::NetworkLayer)
    @assert length(nl.aggregate_vars) == length(ios.inputs)
    @assert length(nl.vertex_vars) <= length(ios.variables)


    @assert all(isequal.(ios.inputs, nl.aggregate_vars)) "Vertex IO system needs to have aggregate variables as inputs"

    agg = nl.aggregate

    f_oop, f_ip, ds = build_io_functions(ios; first_inputs = nl.aggregate_vars, first_dynamic_states = nl.vertex_vars)
    f_vert(dx,x,es,ed,p,t) = f_ip(dx, x, agg(es,ed), p, t)

    @assert all(isequal.(ds[1:length(nl.vertex_vars)], nl.vertex_vars)) "The first IO System variables need to match the vertex_variables"

    ODEVertex(f! = f_vert, dim = length(ds); sym = MTK.tosymbol.(ds))
end

##

@variables ω(t)

KuramotoVertex = IOSystem(
    [D(ω) ~ 1. - total_flow,
    D(ϕ) ~ ω];
    inputs = [total_flow], outputs = [])

nd_k_vertex = build_nd_vertex(KuramotoVertex, KuramotoNetwork)

using LightGraphs

function diffusionedge!(e, v_s, v_d, p, t)
    e[1] = v_s[1] - v_d[1]
    nothing
end

nd_diffusion_edge = StaticEdge(f! = diffusionedge!, dim = 1)
nd = network_dynamics(nd_k_vertex, nd_diffusion_edge, barabasi_albert(10,2))

x0 = randn(20) # random initial conditions
dx0 = similar(x0)
nd(dx0, x0, nothing, 0.)

##

using OrdinaryDiffEq

ode_prob = ODEProblem(nd, x0, (0., 4.))
sol = solve(ode_prob, Tsit5());

##

KuramotoEdge = IOSystem(
    [f ~ sin(ϕ_s - ϕ_d)];
    inputs = [], outputs = [])

function build_static_nd_edge(io_sys::IOSystem, nl::NetworkLayer)
    # Very rudimentary
    edge_eq_idx = [findfirst(x -> isequal(x.lhs, e_var), io_sys.eqs) for e_var in nl.edge_vars]
    @assert length(edge_eq_idx) == length(nl.edge_vars)
    edge_eq = io_sys.eqs[edge_eq_idx]

    os = ODESystem(edge_eq)

    Set(os.states) ⊆ ( Set(nl.vertex_s_vars) ∪ Set(nl.vertex_d_vars) ∪ Set(nl.edge_vars) )

    edge_formulas = [eq.rhs for eq in edge_eq]

    # So MTK includes a "fill zeros" at the start of the inplace function. That throws an error. It's unclear to me why
    # This might be our fault due to the eltype of edge_views being weird in some way...
    f_oop, f_ip = build_function(edge_formulas, nl.vertex_s_vars, nl.vertex_d_vars, os.ps, os.iv; expression = Val{false}, fillzeros = false)
    StaticEdge(f! = f_ip, dim = length(nl.edge_vars), sym = MTK.tosymbol.(nl.edge_vars))
end

nd_k_edge = build_static_nd_edge(KuramotoEdge, KuramotoNetwork)

a = [0.]
nd_k_edge.f!(a, [0., 1.], [1., 1.], nothing, 0.)
@assert a[1] == sin(-1.)

##

nd = network_dynamics(nd_k_vertex, nd_k_edge, barabasi_albert(10,2))

x0 = randn(20) # random initial conditions
dx0 = similar(x0)
nd(dx0, x0, nothing, 0.)


##