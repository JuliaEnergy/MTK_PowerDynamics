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
        # We could also construct from two sets of eqs: dynamics and observed.
        
        state_vars = setdiff(os.states, inputs ∪ outputs)
        # Make sure the outputs and inputs actually occur in the equations
        @assert Set(inputs) ⊆ Set(os.states) "Inputs need to occur in the equations"
        @assert Set(outputs) ⊆ Set(os.states) "Outputs need to occur in the equations"
        @assert isempty(Set(inputs) ∩ Set(outputs)) "Outputs and Inputs need to be disjoint"
        @assert isequal( Set(inputs) ∪ Set(outputs) ∪ Set(state_vars), Set(os.states) )
        

        # We require output variables to occur as left hand side in the equations
        lhs_variables = [eq.lhs for eq in eqs]
        rhs_variables = [eq.rhs for eq in eqs]
        @assert Set(outputs) ⊆ Set(lhs_variables) "Outputs need to be the left hand side of an equation."
      

        # Todo: Make sure the right hand side of output equations does not contain inputs.
        # Otherwise we can not eliminate inputs simply by plugging in the outputs, but need to potentially solve equations.
        # One can always work around this by introducing an intermediate variable. That will show up as a constraint equation
        # in the dynamics.
        # There might be additional constraints we need for our IO System equations, like we don't want any outputs to have two equations.

        # for eq in eqs
        #     if Set(output_vars) ⊆ MTK.vars(eq.lhs) # output eqn
        #         @assert isempty( intersect(Set(input_vars), MTK.vars(eq.rhs)) ) "Make sure the right hand side of output equations does not contain inputs."
        #     end
        # end

        new(eqs, state_vars, inputs, outputs)
    end
end

##

@parameters t
@derivatives D'~t

@variables x(t), u(t), y(t)
@parameters a, b, c, d

ol = IOSystem([D(x) ~ a * x + b * u, y ~ c * x + d * u], inputs = [u], outputs = [y])


##

# A note on types in MTK, explicit type signatures are tricky because there are a lot of types floating around:

os = ODESystem([D(x) ~ a * x + b * u, y ~ c * x + d * u], name=:sys1)
@show typeof(x) # typeof(x) = Num
@show typeof(x.val) # typeof(x) = Term{Real}
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

    @assert allunique(vcat(all_inputs, all_outputs))

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
pc = IOSystem([u_c ~ k_P * y_c]; inputs = [y_c], outputs = [u_c])

## Connect them:

closed_loop = connect([u_c => u, y => y_c], [ol, pc])


##

function build_io_functions(io_sys::IOSystem)
    out_eqs_idx = [findfirst(x -> isequal(x.lhs, o), io_sys.eqs) for o in io_sys.outputs]

    output_equations = io_sys.eqs[out_eqs_idx]
    dynamic_equations = setdiff(io_sys.eqs, output_equations)
    
    os = ODESystem(dynamic_equations)
    
    not_outputs = setdiff(os.states, io_sys.outputs)
    dynamic_states = setdiff(not_outputs, io_sys.inputs)
    dynamic_formulas = [eq.rhs for eq in dynamic_equations] # Here we assume that the left hand side is d/dx.

    # There is a subtlty here with ordering. The equations in dynamic_formulas and the dynamic_states need to match up
    # so the left hand side ofthe equations matches the order of the states. Maybe this is already guaranteed as a consequence of the way ODESystem and setdiff work.
    # the function MTK.vars can extract variables from Differentials on the left hand side.

    f_oop, f_ip = build_function(dynamic_formulas, dynamic_states, io_sys.inputs, os.ps, os.iv; expression = Val{false})
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

f_oop, f_ip = build_io_functions(ol)

f_oop([1], [2], [3,4], 0.)

##

using NetworkDynamics

struct NetworkLayer
    vertex_vars
    edge_vars
    aggregate_vars
    aggregate # The function that aggregates
end

function build_nd_vertex(ios::IOSystem, nl::NetworkLayer)
    @assert length(nl.aggregate_vars) == length(ios.inputs)
    @assert length(nl.vertex_vars) <= length(ios.variables)

    # So ordering is crucial to get all of this right. We need to check and double check our ordering guarantees
    # for building io functions.

    @assert all(isequal.(ios.inputs, nl.aggregate_vars)) "Vertex IO system needs to have aggregate variables as inputs"
    @assert all(isequal.(ios.variables[1:length(nl.vertex_vars)], nl.vertex_vars)) "The first IO System variables need to match the vertex_variables"

    const agg = nl.aggregate

    f_oop, f_ip = build_io_functions(ios)
    f_vert(dx,x,es,ed,p,t) = f_ip(dx, x, agg(es,ed), p, t)

    ODEVertex(f_vert, dim, mass_matrix, sym)
end

##