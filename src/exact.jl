export exactpruning!

# Exact pruning
"""
    exactpruning!(man::AbstractCutPruner, optimizer_constructor;
                  ub=Inf, lb=-Inf, epsilon=1e-5)


Remove dominated cuts in CutPruner `man`.

We use a LP solver to determine whether a cut is dominated or not.

# Arguments
* `man::AbstractCutPruner`
    Cut pruner where to remove cuts
* `optimizer_constructor`
    Optimizer used to solve LP
* `ub::Union{Float64, Vector{Float64}}`
    State x upper bound
* `lb::Union{Float64, Vector{Float64}}`
    State x lower bound
* `epsilon::Float64`
    Pruning's tolerance
"""
function exactpruning!(man::AbstractCutPruner, optimizer_constructor; ub = Inf, lb = -Inf, epsilon = 1e-5)
    K = getdominated(man.A, man.b, man.islb, man.isfun, optimizer_constructor, lb, ub, epsilon)
    removecuts!(man, K)
end

"""Return dominated cuts."""
function getdominated(A, b, islb, isfun, optimizer_constructor, lb, ub, epsilon)
    red = Int[]
    if size(A, 1) == 1
        return red
    end
    for i in 1:size(A, 1)
        if isdominated(A, b, islb, isfun, i, optimizer_constructor, lb, ub, epsilon)
            push!(red, i)
        end
    end
    red
end

"""State whether a cut is dominated with a tolerance epsilon."""
function isdominated(A, b, islb, isfun, k, optimizer_constructor, lb, ub, epsilon)
    # we use JuMP to solve the test
    # For instance, if islb & isfun, the LP becomes:
    # min - y
    #     x ∈ R^n, y ∈ R
    #     y <= (A[k, :] - A[i, :])*x + b[k] - b[i]     ∀ i != k

    # we formulate an equivalent problem as
    # min c'*z
    # s.t H*z <= h

    # get problem's dimension
    ncuts, nx = size(A)
    # allocate arrays
    h = zeros(ncuts - 1)
    H = zeros(ncuts - 1, nx + 1)
    c = zeros(nx + 1)
    c[1] = (islb) ? -1 : 1

    λk = @view A[k, :]

    ic = 0
    for ix in 1:ncuts
        if ix != k
            ic += 1
            δb = b[k] - b[ix]
            @inbounds h[ic] = ((isfun & islb) || (~isfun & ~islb)) ? δb : -δb
            @inbounds H[ic, 1] = (islb) ? 1 : -1

            for jx in 1:nx
                dl = A[ix, jx] - λk[jx]
                @inbounds H[ic, jx + 1] = (islb) ? dl : -dl
            end
        end
    end

    # solve the LP with JuMP
    model = Model(optimizer_constructor)
    z = @variable(model, [1:nx+1])
    set_lower_bound.(z[2:end], lb)
    set_upper_bound.(z[2:end], ub)

    @constraint(model, H * z .≤ h)
    @objective(model, Min, c ⋅ z)
    optimize!(model)

    status = termination_status(model)

    if status == MOI.INFEASIBLE
        return true
    elseif status == MOI.DUAL_INFEASIBLE
        return false
    elseif status == MOI.OPTIMAL
        res = objective_value(model)
        return (islb) ? -res < epsilon : res > -epsilon
    else
        error("Solver returned status $status: $(raw_status(model)).")
    end
end
