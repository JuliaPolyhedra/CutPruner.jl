################################################################################
# Implement abstract type of CutPruner
################################################################################
export AbstractCutPruningAlgo
export CutPruner, AbstractCutPruner
# High-level functions
export addcuts!, ncuts
# Low-level functions
export appendcuts!, replacecuts!, keeponlycuts!, removecuts!

abstract AbstractCutPruningAlgo

"""
A cut pruner maintains a matrix `A` and a vector `b` such that
represents `size(A, 1)` (` == length(b)`) cuts.
Let `a_i` be `A[i,:]` and `β_i` be `b[i]`, the meaning of the cut depends on the sense.
Cuts (A, b) defines the half-space satisfying:
Ax >= b if islb
Ax <= b otherwise
If `sense` is
* `:Min`, then the cut pruner represents the concave polyhedral function `min ⟨a_i, x⟩ + β_i`;
* `:Max`, then the cut pruner represents the convex polyhedral function `max ⟨a_i, x⟩ + β_i`;
* `:≤`, then the cut pruner represents the polyhedra defined by the intersection of the half-space `⟨a_i, x⟩ ≤ β_i`;
* `:≥`, then the cut pruner represents the polyhedra defined by the intersection of the half-space `⟨a_i, x⟩ ≥ β_i`.

Internally, instead of `sense`, the booleans `isfun` and `islb` are stored.
The mapping between `sense` and these two booleans is given by the following table

| `sense` | `isfun` | `islb` |
| ------- | ------- | ------ |
| Min     | true    | false  |
| Max     | true    | true   |
| ≤       | false   | false  |
| ≥       | false   | true   |

"""
abstract AbstractCutPruner{N, T}

function gettype(sense::Symbol)
    if sense == :Min
        true, false
    elseif sense == :Max
        true, true
    elseif sense == :≤
        false, false
    elseif sense == :≥
        false, true
    else
        throw(ArgumentError("Invalid value `$sense' for sense. It should be :Min, :Max, :≤ or :≥."))
    end
end
function getsense(isfun::Bool, islb::Bool)
    if isfun
        islb ? :Max : :Min
    else
        islb ? :≥ : :≤
    end
end
getsense(pruner::AbstractCutPruner) = getsense(isfun(pruner), islb(pruner))
isfun(pruner::AbstractCutPruner) = pruner.isfun
islb(pruner::AbstractCutPruner) = pruner.islb

immutable CutPruner{N, T} end

"""Return whether the CutPruner `man` has any cut."""
function Base.isempty(man::AbstractCutPruner)
    return isempty(man.b)
end

"""Return number of cuts in CutPruner `man`."""
function ncuts(man::AbstractCutPruner)
    return length(man.b)
end

# COMPARISON
hastrust(man::AbstractCutPruner) = true
"""Get current `trust` of CutPruner `man`."""
gettrust(man::AbstractCutPruner) = man.trust

function _indmin(a::Vector, tiebreaker::Vector)
    imin = 1
    for i in 2:length(a)
        if a[i] < a[imin] || (a[i] == a[imin] && tiebreaker[i] < tiebreaker[imin])
            imin = i
        end
    end
    imin
end

"""Remove `num` cuts in CutPruner `man`."""
function choosecutstoremove(man::AbstractCutPruner, num::Int)
    # MergeSort is stable so in case of equality, the oldest cut loose
    # However PartialQuickSort is a lot faster

    trust = gettrust(man)
    if num == 1
        [_indmin(trust, man.ids)]
    else
        # /!\ PartialQuickSort is unstable, here it does not matter
        function _lt(i, j)
            # If cuts have same trust, remove oldest cut
            if trust[i] == trust[j]
                man.ids[i] < man.ids[j]
            # Else, remove cuts with lowest trust
            else
                trust[i] < trust[j]
            end
        end
        # Return index of `num` cuts with lowest trusts
        sort(1:length(trust), alg=PartialQuickSort(num), lt=_lt)[1:num]
    end
end

"""Test if cut `i` is better than `newcuttrust`."""
isbetter(man::AbstractCutPruner, i::Int, mycut::Bool) = gettrust(man)[i] > initialtrust(man, mycut)

# CHANGE

function getnreplaced(man::AbstractCutPruner, R, ncur, nnew, nmycut)
    # Check if some new cuts should be ignored
    take = man.maxncuts - ncur
    # get number of cuts to remove
    nreplaced = length(R)
    # Start:
    # |      | R      j|
    # | take |         |
    # |  mycut  |      |
    # End:
    # |      | R   j   |
    # | take       |   |
    # |  mycut  |      |
    while take + length(R) - nreplaced < nnew
        # I first try to see if the nmycut cuts generated by me can be taken
        # because they have better trust than the nnew-nmycut others
        if isbetter(man, R[nreplaced], take < nmycut)
            nreplaced -= 1
        else
            take += 1
        end
    end
    take, nreplaced
end

# Add cuts Ax >= b
# If mycut then the cut has been added because of one of my trials
function addcuts!{N, T}(man::AbstractCutPruner{N, T},
                     A::AbstractMatrix{T},
                     b::AbstractVector{T},
                     mycut::AbstractVector{Bool})
    # get current number of cuts:
    ncur = ncuts(man)
    nincumbents = size(A, 1)

    # check redundancy
    redundants = checkredundancy(man.A, man.b, A, b, man.isfun, man.islb, man.TOL_EPS)
    if !isempty(redundants)
        tokeep = setdiff(collect(1:nincumbents), redundants)

        # if all cuts are redundants, then do nothing:
        if length(tokeep) == 0
            return zeros(Int, nincumbents)
        end
        A = A[tokeep, :]
        b = b[tokeep]
        mycut = mycut[tokeep]
    end

    # get number of new cuts in A:
    nnew = size(A, 1)
    @assert length(mycut) == length(b) == nnew
    @assert nnew > 0

    if man.maxncuts == -1 || ncur + nnew <= man.maxncuts
        # If enough room, just append cuts
        status = ncur + (1:nnew)
        appendcuts!(man, A, b, mycut)
    else
        # Otherwise, we need need to remove some cuts

        # get indexes of cuts with lowest trusts:
        R = choosecutstoremove(man, ncur + nnew - man.maxncuts)

        nmycut = sum(mycut)
        take, nreplaced = getnreplaced(man, R, ncur, nnew, nmycut)

        # Cuts that will be replaced
        R = @view R[1:nreplaced]
        # Nowe we split A, b into
        # * A, b   : cuts to be pushed
        # * Ar, br : cuts to replace old cuts
        # * _, _   : cuts removed
        if nreplaced == nnew
            status = R
            Ar = A
            br = b
            mycutr = mycut
            A = similar(A, 0, size(A, 2))
            b = similar(b, 0)
            mycut = similar(mycut, 0)
        else
            status = zeros(Int, nnew)
            if take < nnew
                # Remove ignored cuts
                takemycut = min(take, nmycut)
                takenotmycut = take - takemycut
                takeit = zeros(Bool, nnew)
                for i in 1:nnew
                    ok = false
                    if mycut[i]
                        if takemycut > 0
                            takemycut -= 1
                            ok = true
                        end
                    else
                        if takenotmycut > 0
                            takenotmycut -= 1
                            ok = true
                        end
                    end
                    if ok
                        takeit[i] = true
                    end
                end

                takeit = find(takeit)
            else
                takeit = collect(1:nnew)
            end
            @assert take == length(takeit)
            replaced = takeit[1:nreplaced]
            pushed = takeit[(nreplaced+1):end]

            status[replaced] = R
            Ar = @view A[replaced,:]
            br = @view b[replaced]
            mycutr = @view mycut[replaced]
            status[pushed] = ncur + (1:length(pushed))
            appendcuts!(man, (@view A[pushed,:]), (@view b[pushed]), (@view mycut[pushed]))
        end
        if nreplaced > 0
            man.A[R, :] = Ar
            man.b[R] = br
            replacecuts!(man, R, Ar, br, mycutr)
        end
    end

    status
end

# Low-level functions
function _keeponlycuts!(man::AbstractCutPruner, K::AbstractVector{Int})
    man.A = man.A[K, :]
    man.b = man.b[K]
    man.ids = man.ids[K]
    if hastrust(man)
        man.trust = gettrust(man)[K]
    end
end

"""Keep only cuts whose indexes are in Vector `K` in CutPruner `man`. If `K` is not sorted, the cuts will change their order accordingly."""
function keeponlycuts!(man::AbstractCutPruner, K::AbstractVector{Int})
    _keeponlycuts!(man, K)
end

"""Remove cuts whose indexes are in Vector `K` in CutPruner `man`."""
function removecuts!(man::AbstractCutPruner, K::AbstractVector{Int})
    keeponlycuts!(man, setdiff(1:ncuts(man), K))
end

function _replacecuts!(man::AbstractCutPruner, K::AbstractVector{Int}, A, b)
    man.A[K, :] = A
    man.b[K] = b
    man.ids[K] = newids(man, length(K))
end

"""Replace cuts at indexes in `K` by cuts in (A, b, mycut) in CutPruner `man`."""
function replacecuts!(man::AbstractCutPruner, K::AbstractVector{Int}, A, b, mycut::AbstractVector{Bool})
    _replacecuts!(man, K, A, b)
    man.trust[K] = initialtrusts(man, mycut)
end

function _appendcuts!(man::AbstractCutPruner, A, b)
    man.A = [man.A; A]
    man.b = [man.b; b]
    append!(man.ids, newids(man, length(b)))
end

"""Append cuts (A, b, mycut) in CutPruner `man`."""
function appendcuts!(man::AbstractCutPruner, A, b, mycut::AbstractVector{Bool})
    _appendcuts!(man, A, b)
    append!(man.trust, initialtrusts(man, mycut))
end

# Unexported utilities

function newids(man::AbstractCutPruner, n::Int)
    (man.id+1):(man.id += n)
end

"""Get a Vector of Float64 specifying the initial trusts of `mycut`."""
function initialtrusts(man::AbstractCutPruner, mycut::AbstractVector{Bool})
    Float64[initialtrust(man, mc) for mc in mycut]
end

function normalizedcut{T}(A::AbstractMatrix{T}, b::AbstractVector{T}, k::Int, isfun::Bool, tol::Float64)
    a = @view A[k, :]
    β = b[k]
    na = norm(a, 2)
    if isfun || na < tol
        a, β
    else
        a / na, β / na
    end
end

"""
Check redundant cuts. Return index of redundant cuts in `Anew`.

$(SIGNATURES)

"""
function checkredundancy{T}(A::AbstractMatrix{T}, b::AbstractVector{T},
                            Anew::AbstractMatrix{T}, bnew::AbstractVector{T},
                            isfun::Bool, islb::Bool, tol::Float64)
    # index of redundants cuts
    redundants = Int[]
    # number of new lines
    nnew = size(Anew, 1)

    for kk in 1:nnew
        a, β = normalizedcut(Anew, bnew, kk, isfun, tol)
        chk, indk = isinside(A, b, a, isfun, tol)
        if chk
            ared, βred = normalizedcut(A, b, indk, isfun, tol)
            if islb ? β <= βred+tol : β+tol >= βred
                push!(redundants, kk)
            end
        end
    end

    redundants
end


"""Check if `λ` is a line of matrix `A`. `λ` might not have the same `eltype` as `A` and `b` as it might have been scaled by `normalizecut`."""
function isinside{T}(A::AbstractMatrix{T}, b::AbstractVector{T}, λ::AbstractVector, isfun::Bool, tol::Float64)
    nlines = size(A, 1)

    check = false
    k = 0
    while ~check && k < nlines
        k += 1
        a, β = normalizedcut(A, b, k, isfun, tol)
        check = norm(a - λ, Inf) < tol
    end
    check, k
end
