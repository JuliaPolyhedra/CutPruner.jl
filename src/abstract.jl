################################################################################
# Implement abstract type of CutPruner
################################################################################
export AbstractCutPruner, addcuts!, start!, isstarted

abstract AbstractCutPruner{S}

"""Return number of cuts in CutPruner `man`."""
function ncuts(man::AbstractCutPruner)
    return length(get(man.cuts_de))
end

function isfeasibilitycut(man::AbstractCutPruner, cut)
    if length(man.σs) < length(man.ρs)
        cut in man.σs
    else
        !(cut in man.ρs)
    end
end

function init!(man::AbstractCutPruner, mycut_d::Vector{Bool}, mycut_e::Vector{Bool})
    mycut = [mycut_d; mycut_e]
    man.trust = Float64[initialtrust(man, mc) for mc in mycut]
    man.ids = newids(man, length(mycut))
end

function start!{S}(man::AbstractCutPruner{S}, ncols::Integer)
    start!(man, Matrix{S}(0, ncols), Matrix{S}(0, ncols), S[], S[], Bool[], Bool[])
end

"""Start CutPruner `man`."""
function start!{S}(man::AbstractCutPruner{S},
                   cuts_D::AbstractMatrix{S},
                   cuts_E::AbstractMatrix{S},
                   cuts_d::AbstractVector{S},
                   cuts_e::AbstractVector{S},
                   mycut_d::AbstractVector{Bool},
                   mycut_e::AbstractVector{Bool})
    man.nσ = length(cuts_d)
    man.nρ = length(cuts_e)
    man.cuts_DE = [cuts_D; cuts_E]
    man.cuts_de = [cuts_d; cuts_e]
    man.σs = collect(1:man.nσ)
    man.ρs = collect(man.nσ+(1:man.nρ))
    init!(man, mycut_d, mycut_e)
end

"""State if CutPruner `man` was initialized."""
function isstarted(man::AbstractCutPruner)
    @assert isnull(man.cuts_DE) == isnull(man.cuts_de)
    !isnull(man.cuts_DE)
end

# COMPARISON
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

# Add cut ax >= β
# FIXME: Ax >= b ?
# If fc then it is a feasibility cut, otherwise it is an optimality cut
# If mycut then the cut has been added because of one of my trials
function addcuts!{S}(man::AbstractCutPruner{S},
                     A::AbstractMatrix{S},
                     b::AbstractVector{S},
                     isfc::Bool,
                     mycut::Vector{Bool})
    # get number of new cuts in A:
    nnew = size(A, 1)
    @assert length(mycut) == length(b) == nnew
    @assert nnew > 0
    status = Symbol[:Pushed for i in 1:nnew]

    # If not enough room, need to remove some cuts
    if man.maxncuts != -1 && ncuts(man)+nnew > man.maxncuts
        # get indexes of cuts with lowest trusts:
        J = choosecutstoremove(man, ncuts(man) + nnew - man.maxncuts)

        # Check if some new cuts should be ignored
        take = man.maxncuts - ncuts(man)
        # get number of cuts to remove
        j = length(J)
        nmycut = sum(mycut)
        while j > 0 && take < nnew
            # I first try to see if the nmycut cuts generated by me can be taken
            # because they have better trust than the nnew-nmycut others
            if isbetter(man, J[j], take < nmycut)
                j -= 1
            else
                take += 1
            end
        end
        J = @view J[1:j]
        if take < size(A, 1) || length(J) < length(b)
            if take < size(A, 1)
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
                    else
                        status[i] = :Ignored
                    end
                end

                takeit = find(takeit)
            else
                takeit = collect(1:nnew)
            end
            nreplaced = min(length(J), length(takeit))
            replaced = takeit[1:nreplaced]
            pushed = takeit[(nreplaced+1):end]

            status[replaced] = :Replaced
            Ar = A[replaced,:]
            br = b[replaced]
            mycutr = mycut[replaced]
            A = A[pushed,:]
            b = b[pushed]
            mycut = mycut[pushed]
        else
            Ar = A
            br = b
            mycutr = mycut
            A = similar(A, 0, size(A, 2))
            b = similar(b, 0)
            mycut = similar(mycut, 0)
            status[:] = :Replaced
        end
        if !isempty(br)
            if length(br) == length(J)
                js = J
                J = []
            else
                @assert length(br) < length(J)
                js = J[(length(J) - length(br) + 1):end]
            end
            get(man.cuts_DE)[js,:] = Ar
            get(man.cuts_de)[js] = br
            replacecuts!(man, js, mycutr)
            cutadded = true
            needupdate_σsρs = reduce(|, false, Bool[isfc $ isfeasibilitycut(man, j) for j in js])
        else
            cutadded = false
            needupdate_σsρs = false
        end

        if !isempty(J) || needupdate_σsρs
            keep = ones(Bool, ncuts(man))
            keep[J] = false
            K = find(keep)
            isσcut = zeros(Bool, ncuts(man))
            isσcut[man.σs] = true
            if cutadded
                isσcut[js] = isfc
            end
            isσcut = isσcut[K]
            man.σs = (1:length(isσcut))[isσcut]
            man.ρs = (1:length(isσcut))[!isσcut]
            man.nσ = length(man.σs)
            man.nρ = length(man.ρs)
        end

        if !isempty(J)
            # TODO: dry these two lines in keeponly to avoid side effect?
            man.cuts_DE = get(man.cuts_DE)[K,:]
            man.cuts_de = get(man.cuts_de)[K]
            keeponly!(man, K)
        end
    elseif !isempty(b)
        # Just append cuts
        if isfc
            append!(man.σs, man.nσ + man.nρ + (1:nnew))
            man.nσ += nnew
        else
            append!(man.ρs, man.nσ + man.nρ + (1:nnew))
            man.nρ += nnew
        end
        man.cuts_DE = [get(man.cuts_DE); A]
        man.cuts_de = [get(man.cuts_de); b]
        pushcuts!(man, mycut)
    end

    status
end

"""Keep only cuts whose indexes are in Vector `K`."""
function keeponly!(man::AbstractCutPruner, K::Vector{Int})
    man.trust = man.trust[K]
end

"""Get a Vector of Float64 specifying the initial trusts of `mycut`."""
function initialtrusts(man::AbstractCutPruner, mycut::Vector{Bool})
    Float64[initialtrust(man, mc) for mc in mycut]
end

"""Reset trust of cuts with indexes in `js`."""
function replacecuts!(man::AbstractCutPruner, js::AbstractVector{Int}, mycut::Vector{Bool})
    man.trust[js] = initialtrusts(man, mycut)
    man.ids[js] = newids(man, length(js))
end

function pushcuts!(man::AbstractCutPruner, mycut::Vector{Bool})
    append!(man.trust, initialtrusts(man, mycut))
    append!(man.ids, newids(man, length(mycut)))
end

function newids(man::AbstractCutPruner, n::Int)
    (man.id+1):(man.id += n)
end
