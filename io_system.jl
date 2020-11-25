## Setup

# A general note for handling equations and terms: == is overloaded and does not return bool,
# Many of the normal ways of handling arrays don't work. Use isequal instead. Set based manipulations
# seem fine.

using ModelingToolkit
const MTK = ModelingToolkit

"""
IOSystems model systems of the form:

dx/dt = f(x,i)

o = g(x)

"""
struct IOSystem # variables are all that occur in eqs, inputs and outputs are subsets of the variables.
    eqs
    variables
    inputs
    outputs
    function IOSystem(eqs; inputs, outputs)
        os = ODESystem(eqs, name=:sys) # We use ODESystem to analyze the equations for us.
        # We could also construct from two sets of eqs: dynamics and observed.
        
        # Make sure the outputs and inputs actually occur in the equations
        @assert Set(inputs) ⊆ Set(os.states) "Inputs need to occur in the equations"
        @assert Set(outputs) ⊆ Set(os.states) "Outputs need to occur in the equations"

        # We require output variables to occur as left hand side in the equations
        lhs_variables = Set([eq.lhs for eq in eqs])
        @assert Set(outputs) ⊆ Set(lhs_variables) "Outputs need to be the left hand side of an equation."

        @assert isempty(Set(inputs) ∩ Set(outputs)) "Outputs and Inputs need to be disjoint"

        # Todo: Make sure the right hand side of output equations does not contain inputs.
        # There might be additional constraints we need for our IO System equations, like we don't want any outputs to have two equations.
        
        new(eqs, os.states, inputs, outputs)
    end
end

##

@parameters t
@derivatives D'~t

@variables x(t), u(t), y(t)
@parameters a, b, c

ol = IOSystem([D(x) ~ a * x + b * u, y ~ c * x], inputs = [u], outputs = [y])

##

function connect(connections, io_systems::Array{IOSystem,1})
    # The connections are pairs o => i that tell us which output to connect to which input.

    # For now simply concatenate everything. TODO: make sure things are suitably disjoint!
    all_inputs = vcat((ios.inputs for ios in io_systems)...)
    all_outputs = vcat((ios.outputs for ios in io_systems)...)
    output_set = Set(all_outputs)
    
    connected_inputs = [conn.second for conn in connections]
    connected_outputs = [conn.first for conn in connections]
    
    @assert Set(connected_inputs) ⊆ Set(all_inputs) "Must connect to inputs"
    @assert Set(connected_outputs) ⊆ Set(all_outputs) "Must connect from outputs"
    
    eqs = vcat((ios.eqs for ios in io_systems)...)

    # as == is overloaded for variables one has to use isequal(x, y) to compare Terms instead.
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
    dynamic_formulas = [eq.rhs for eq in dynamic_equations]
    g_oop, g_ip = build_function(dynamic_formulas, dynamic_states, io_sys.inputs, os.ps, os.iv; expression = Val{false})
end

##

g_oop, g_ip = build_io_functions(ol)

g_oop([1], [2], [3,4], 0.)
