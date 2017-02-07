export DecayCutPruner

"""
$(TYPEDEF)

Removes the cuts with lower trust where the trust is initially `newcuttrust + bonus` and is updated using `trust -> λ * trust + used` after each optimization done with it.
The value `used` is 1 if the cut was used and 0 otherwise.
It has a bonus equal to `mycutbonus` if the cut was generated using a trial given by the problem using this cut.
We say that the cut was used if its dual value is nonzero.
"""
type DecayCutPruner{S} <: AbstractCutPruner{S}
    # used to generate cuts
    cuts_DE::Nullable{AbstractMatrix{S}}
    cuts_de::Nullable{AbstractVector{S}}

    nσ::Int
    nρ::Int
    σs::Vector{Int}
    ρs::Vector{Int}

    maxncuts::Int

    trust::Vector{Float64}
    ids::Vector{Int}
    id::Int

    λ::Float64
    newcuttrust::Float64
    mycutbonus::Float64

    function DecayCutPruner(maxncuts::Int, λ=0.9, newcuttrust=0.8, mycutbonus=1)#newcuttrust=(1/(1/0.9-1))/2, mycutbonus=(1/(1/0.9-1))/2)
        new(nothing, nothing, 0, 0, Int[], Int[], maxncuts, Float64[], Int[], 0, λ, newcuttrust, mycutbonus)
    end
end

DecayCutPruner(maxncuts::Int, λ=0.9, newcuttrust=(1/(1/0.9-1))/2, mycutbonus=(1/(1/0.9-1))/2) = DecayCutPruner{Float64}(maxncuts, λ, newcuttrust, mycutbonus)

function clone{S}(man::DecayCutPruner{S})
    DecayCutPruner{S}(man.maxncuts, man.λ, man.newcuttrust, man.mycutbonus)
end

# COMPARISON
function updatestats!(man::DecayCutPruner, σρ)
    if ncuts(man) > 0
        man.trust *= man.λ
        man.trust[σρ .> 1e-6] += 1
    end
end

function initialtrust(man::DecayCutPruner, mycut)
    if mycut
        man.newcuttrust + man.mycutbonus
    else
        man.newcuttrust
    end
end

function isbetter(man::DecayCutPruner, i::Int, mycut::Bool)
    if mycut
        # If the cut has been generated, that means it is useful
        false
    else
        # The new cut has initial trust initialtrust(man, false)
        # but it is a bit disadvantaged since it is new so
        # as we advantage the new cut if mycut == true,
        # we advantage this cut by taking initialtrust(man, true)
        # with true instead of false
        man.trust[i] > initialtrust(man, mycut)
    end
end
