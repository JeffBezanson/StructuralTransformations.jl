module StructuralTransformations

using ModelingToolkit
using ModelingToolkit: ODESystem, var_from_nested_derivative, Differential, states, equations, vars, Symbolic, symtype

function lower_varname_dae(var::Term, idv, order)
    order == 0 && return var
    name = Symbol(nameof(var.op), :ˍˍ, string(idv)^order)
    return Sym{symtype(var.op)}(name)(var.args[1])
end

# V-nodes `[x_1, x_2, x_3, ..., dx_1, dx_2, ..., y_1, y_2, ...]` where `x`s are
# differential variables and `y`s are algebraic variables.
function get_vnodes(sys)
    dxvars = []
    eqs = equations(sys)
    edges = map(_->Int[], 1:length(eqs))
    for (i, eq) in enumerate(eqs)
        if eq.lhs isa Symbolic
            # Make sure that the LHS is a first order derivative of a var.
            @assert eq.lhs.op isa Differential
            @assert !(eq.lhs.args[1] isa Differential) # first order

            push!(dxvars, eq.lhs)
            # For efficiency we note down the diff edges here
            push!(edges[i], length(dxvars))
        end
    end

    xvars = (first ∘ var_from_nested_derivative).(dxvars)
    algvars  = setdiff(states(sys), xvars)
    return xvars, dxvars, edges, algvars
end

function sys2bigraph(sys)
    xvars, dxvars, edges, algvars = get_vnodes(sys)
    xvar_offset = length(xvars)
    algvar_offset = 2xvar_offset
    for edge in edges
        isempty(edge) || (edge .+= xvar_offset)
    end

    eqs = equations(sys)
    for (i, eq) in enumerate(eqs)
        # T or D(x):
        # We assume no derivatives appear on the RHS at this point
        vs = vars(eq.rhs)
        for v in vs
            for (j, target_v) in enumerate(xvars)
                if isequal(v, target_v)
                    push!(edges[i], j)
                end
            end
            for (j, target_v) in enumerate(algvars)
                if isequal(v, target_v)
                    push!(edges[i], j+algvar_offset)
                end
            end
        end
    end

    fullvars = [xvars; dxvars; algvars] # full list of variables
    vars_asso = Int[(1:xvar_offset) .+ xvar_offset; zeros(Int, length(fullvars) - xvar_offset)] # variable association list
    return edges, fullvars, vars_asso
end

function match_equation!(edges, i, assign, active, vcolor=falses(length(active)), ecolor=falses(length(edges)))
    # `edge[active]` are active edges
    # i: equations
    # j: variables
    # assign: assign[j] == i means (i-j) is assigned
    #
    # color the equation
    ecolor[i] = true
    # if a V-node j exists s.t. edge (i-j) exists and assign[j] == 0
    for j in edges[i]
        if active[j] && assign[j] == 0
            assign[j] = i
            return true
        end
    end
    # for every j such that edge (i-j) exists and j is uncolored
    for j in edges[i]
        (active[j] && !vcolor[j]) || continue
        # color the variable
        vcolor[j] = true
        if match_equation!(edges, assign[j], assign, active, vcolor, ecolor)
            assign[j] = i
            return true
        end
    end
    return false
end

function matching(edges, nvars, active=trues(nvars))
    assign = zeros(Int, nvars)
    for i in 1:length(edges)
        match_equation!(edges, i, assign, active)
    end
    return assign
end

# Naive subtree matching, we can make the dict have levels
# Going forward we should look into storing the depth in Terms
function walk_and_substitute(expr, substitution_dict)
    if haskey(substitution_dict, expr)
        return substitution_dict[expr]
    elseif expr isa Term
        return Term(walk_and_substitute(expr.op, substitution_dict),
                         map(x->walk_and_substitute(x, substitution_dict), expr.args))
    elseif expr isa Equation
        error("TODO: can rhs contain equations?")
    else
        return expr
    end
end

function pantelides_reassemble(sys, vars, vars_asso, eqs_asso, assign)
    # Step 1: write derivative equations
    in_eqs = equations(sys)
    out_eqs = Vector{Any}(undef, length(eqs_asso))
    fill!(out_eqs, nothing)
    out_eqs[1:length(in_eqs)] .= in_eqs

    d_dict = Dict(zip(vars, 1:length(vars)))
    D = ModelingToolkit.Differential(sys.iv)
    for (i, e) in enumerate(eqs_asso)
        if e === 0
            continue
        end
        # LHS variable is looked up from vars_asso
        # the vars_asso[i]-th variable is the differentiated version of var at i
        eq = out_eqs[i]
        lhs = if !(eq.lhs isa Symbolic)
            0
        elseif eq.lhs isa Term && eq.lhs.op isa Differential
            # look up the variable that represents D(lhs)
            @assert !(eq.lhs.args[1] isa Differential) # first order
            i = get(d_dict, eq.lhs.args[1], nothing)
            if i !== nothing
                D(vars[vars_asso[i]])
            else
                D(eq.lhs)
            end
        else
            D(eq.lhs)
        end
        rhs = ModelingToolkit.expand_derivatives(D(eq.rhs))
        out_eqs[e] = lhs ~ rhs
        substitution_dict = Dict(x.lhs => x.rhs for x in out_eqs if x !== nothing && x.lhs isa Symbolic)
        sub_rhs = walk_and_substitute(rhs, substitution_dict)
        out_eqs[e] = lhs ~ sub_rhs
    end
    @assert !any(x->x === nothing, out_eqs)

    final_vars = filter(x->!(x.op isa Differential), vars)
    final_eqs = map(identity, out_eqs[sort(filter(!iszero, assign))])

    # remove clashing equations (from order lowering vs index reduction)
    return ODESystem(final_eqs, sys.iv, final_vars, sys.ps)
end

function pantelides(sys::ODESystem; kwargs...)
    edges, fullvars, vars_asso = sys2bigraph(sys)
    return pantelides!(edges, fullvars, vars_asso, sys.iv; kwargs...)
end

function pantelides!(edges, vars, vars_asso, iv; maxiter = 8)
    neqs = length(edges)
    nvars = length(vars)
    assign = zeros(Int, nvars)
    eqs_asso = fill(0, neqs)
    neqs′ = neqs
    D = Differential(iv)
    for k in 1:neqs′
        i = k
        pathfound = false
        # In practice, `maxiter=8000` should never be reached, otherwise, the
        # index would be on the order of thousands.
        for iii in 1:maxiter
            # run matching on (dx, y) variables
            #
            # the derivatives and algebraic variables are zeros in the variable
            # association list
            active = vars_asso .== 0
            vcolor = falses(nvars)
            ecolor = falses(neqs)
            pathfound = match_equation!(edges, i, assign, active, vcolor, ecolor)
            pathfound && break # terminating condition
            # for every colored V-node j
            for j in eachindex(vcolor); vcolor[j] || continue
                # introduce a new variable
                nvars += 1
                # introduce a differential variable (x)
                lhsj = vars[j]
                vj, order = var_from_nested_derivative(lhsj)
                _newvarj = lower_varname_dae(vj, iv, order)
                newvarj = _newvarj
                vars[j] = newvarj
                vars_asso[j] = nvars
                # introduce a derivative variable (dx)
                push!(vars, D(newvarj))
                push!(vars_asso, 0)
                push!(assign, 0)
            end

            # for every colored E-node l
            for l in eachindex(ecolor); ecolor[l] || continue
                neqs += 1
                # create new E-node
                push!(edges, copy(edges[l]))
                # create edges from E-node `neqs` to all V-nodes `j` and
                # `vars_asso[j]` s.t. edge `(l-j)` exists
                for j in edges[l]
                    if !(vars_asso[j] in edges[neqs])
                        push!(edges[neqs], vars_asso[j])
                    end
                end
                push!(eqs_asso, 0)
                eqs_asso[l] = neqs
            end

            # for every colored V-node j
            for j in eachindex(vcolor); vcolor[j] || continue
                assign[vars_asso[j]] = eqs_asso[assign[j]]
            end
            i = eqs_asso[i]
        end # for _ in 1:maxiter
        pathfound || error("maxiter=$maxiter reached! File a bug report if your system has a reasonable index (<100), and you are using the default `maxiter`. Try to increase the maxiter by `pantelides(sys::ODESystem; maxiter=1_000_000)` if your system has an incredibly high index and it is truly extremely large.")
    end # for k in 1:neqs′
    return edges, assign, vars_asso, eqs_asso, vars
end

function dae_index_lowering(sys::ODESystem; kwargs...)
    edges, fullvars, vars_asso = sys2bigraph(sys)
    edges, assign, vars_asso, eqs_asso, vars = pantelides!(edges, fullvars, vars_asso, sys.iv; kwargs...)
    return pantelides_reassemble(sys, vars, vars_asso, eqs_asso, assign)
end

# Debug use only
print_bigraph(sys, vars, edges) = print_bigraph(stdout, sys, vars, edges)
function print_bigraph(io::IO, sys, vars, edges)
    println(io, "Equations:")
    eqs = equations(sys)
    foreach(x->println(io, x), [i => eqs[i] for i in 1:length(eqs)])
    for (i, edge) in enumerate(edges)
        println(io, "\nEq $i has:")
        print(io, '[')
        for e in edge
            print(io, "$(vars[e]), ")
        end
        print(io, ']')
    end
    return nothing
end


end
