export LevelOneCutPruner, LevelOnePruningAlgo


type LevelOnePruningAlgo <: AbstractCutPruningAlgo
    # maximum number of cuts
    maxncuts::Int
    function LevelOnePruningAlgo(maxncuts::Int)
        new(maxncuts)
    end
end


"""
$(TYPEDEF)

Removes the cuts with lower trust where the trust is: nused / nwith + bonus
where the cut has been used `nused` times amoung `nwith` optimization done with it.
We say that the cut was used if its dual value is nonzero.
It has a bonus equal to `mycutbonus` if the cut was generated using a trial given by the problem using this cut.
If `nwidth` is zero, `nused/nwith` is replaced by `newcuttrust`.
"""
type LevelOneCutPruner{N, T} <: AbstractCutPruner{N, T}
    # used to generate cuts
    isfun::Bool
    islb::Bool
    A::AbstractMatrix{T}
    b::AbstractVector{T}

    maxncuts::Int

    trust::Vector{Float64}
    ids::Vector{Int} # small id means old
    id::Int # current id

    #set of states where cut k is active
    territories::Vector{Vector{Tuple{Int64, T}}}
    nstates::Int
    states::Array{T, 2}

    # tolerance to check redundancy between two cuts
    TOL_EPS::Float64

    function LevelOneCutPruner(sense::Symbol, maxncuts::Int, tol=1e-6)
        isfun, islb = gettype(sense)
        new(isfun, islb, spzeros(T, 0, N), T[], maxncuts, Tuple{Int64, T}[], Int[], 0, [], 0, zeros(T, 0, N), tol)
    end
end

(::Type{CutPruner{N, T}}){N, T}(algo::LevelOnePruningAlgo, sense::Symbol) = LevelOneCutPruner{N, T}(sense, algo.maxncuts)

getnreplaced(man::LevelOneCutPruner, R, ncur, nnew, mycut) = nnew, length(R)

"""Update territories with cuts previously computed during backward pass.

$(SIGNATURES)

# Arguments
* `man::LevelOneCutPruner`
* `position::Array{T, 2}`
    New visited positions
"""
function updatestats!{T}(man::LevelOneCutPruner, position::Array{T, 2})
    # get number of new positions to analyse:
    nx = size(position, 1)

    for i in 1:nx
        addstate!(man, position[i, :])
    end

    updatetrust!(man)
end


"""Add a new state to test and accordingly update territories of each cut.

$(SIGNATURES)

"""
function addstate!(man::LevelOneCutPruner, x::Vector)
    # update number of states
    man.nstates += 1
    # Add `x` to the list of visited state:
    man.states = vcat(man.states, x')

    giveterritory!(man, man.nstates, x)
end

function giveterritory!(man::LevelOneCutPruner, ix::Int, x::Vector=man.states[ix, :])
    # Get cut which is active at point `x`:
    bcost, bcuts = optimalcut(man, x)
    # Add `x` with index nstates  to the territory of cut with index `bcuts`:
    push!(man.territories[bcuts], (ix, bcost))
end

"""Find active cut at point `xf`.

$(SIGNATURES)

# Arguments
* `man::LevelOneCutPruner`:
    CutPruner
* `xf::Vector{Float64}`:

# Return
`bestcost::Float64`
    Value of supporting cut at point `xf`
`bestcut::Int64`
    Index of supporting cut at point `xf`
"""
function optimalcut{T}(man::LevelOneCutPruner,
                       xf::Vector{T})
    bestcost = -Inf::Float64
    bestcut = -1
    dimstates = length(xf)
    nc = ncuts(man)

    @inbounds for i in 1:nc
        cost = cutvalue(man, i, xf)
        if cost > bestcost
            bestcost = cost
            bestcut = i
        end
    end
    return bestcost, bestcut
end


"""Update territories (i.e. the set of tested states where
    a given cut is active) considering new cut given by index `indcut`.

$(SIGNATURES)

# Arguments
* `man::LevelOneCutPruner`:
* `indcut::Int64`:
    new cut index
"""
function updateterritory!(man::LevelOneCutPruner, indcut::Int64)
    @assert length(man.territories) == ncuts(man)
    for k in 1:ncuts(man)
        if k == indcut
            continue
        end
        todelete = []
        for (num, (ix, cost)) in enumerate(man.territories[k])
            x = man.states[ix, :]

            costnewcut = cutvalue(man, indcut, x)

            if costnewcut > cost
                push!(todelete, num)
                push!(man.territories[indcut], (ix, costnewcut))
            end
        end
        deleteat!(man.territories[k], todelete)
    end
end


"""
Get value of cut with index `indc` at point `x`.

$(SIGNATURES)

# Arguments
- `man::LevelOneCutPruner`
    Approximation of the value function as linear cuts
- `indc::Int`
    Index of cut to consider
- `x::Vector`
    Coordinates of state

# Return
`cost::Float64`
    Value of cut `indc` at point `x`.
    If `man` is a polyhedral function, then it is the value of the cut at `x`,
    otherwise, it is the distance between `x` and the cut.
    As a rule of thumb, the higher the `cutvalue` is, the less it is redundant.
"""
function cutvalue{T}(man::LevelOneCutPruner, indc::Int, x::Vector{T})
    β = man.b[indc]
    a = @view man.A[indc, :]
    ax = dot(a, x)
    cost = isfun(man) ? ax + β : (β - ax) / norm(a, 2)
    islb(man) ? cost : -cost
end

flength(a)::Float64 = length(a)

function updatetrust!(man)
    @assert length(man.territories) == ncuts(man)
    if ncuts(man) == length(man.trust)
        # Avoid new allocation. Avoiding this allocation is the whole point of
        # doint updatetrust! instead of doing trust = nothing
        for i in 1:ncuts(man)
            man.trust[i] = length(man.territories[i])
        end
    else
        man.trust = flength.(man.territories)
    end
end

function replacecuts!{N, T}(man::LevelOneCutPruner{N, T}, K::AbstractVector{Int}, A, b, mycut::AbstractVector{Bool})
    @assert length(man.territories) == ncuts(man)
    # FIXME If K is 1:ncuts, then checkconsistency will be true and trust will not be recomputed by gettrust
    _replacecuts!(man, K, A, b)
    # Do not do view here since will will modify the entries
    freeterritories = man.territories[K]
    for k in K
        man.territories[k] = Tuple{Int64, T}[]
    end
    for k in K
        updateterritory!(man, k)
    end
    for freet in freeterritories
        for (ik, _) in freet
            giveterritory!(man, ik)
        end
    end
    updatetrust!(man)
    @assert length(man.territories) == ncuts(man)
end


"""Push new cut in CutPruner `man`."""
function appendcuts!{N, T}(man::LevelOneCutPruner{N, T}, A, b, mycut::AbstractVector{Bool})
    @assert length(man.territories) == ncuts(man)
    oldncuts = ncuts(man)
    _appendcuts!(man, A, b)
    nnew = length(b)
    man.territories = vcat(man.territories, [Tuple{Int64, T}[] for _ in 1:nnew])
    for k in oldncuts+(1:nnew)
        updateterritory!(man, k)
    end
    updatetrust!(man)
    @assert length(man.territories) == ncuts(man)
end
