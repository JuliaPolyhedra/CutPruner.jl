
using CutPruners, Clp
using Base.Test

for algo in [AvgCutPruningAlgo(20), DecayCutPruningAlgo(20), DeMatosPruningAlgo(20)]
    @testset "Exact pruning" begin
        pruner = CutPruner{2, Int}(algo, :Max)
        # add 10 cuts in a row
        for i in 1:10
            addcuts!(pruner, [1 0], [i], [true])
        end

        exactpruning!(pruner, ClpSolver())
        # normally the exact pruning saves only the last cuts
        @test pruner.b == [10]
        # ... and add another set of cuts ...
        for i in 1:10
            addcuts!(pruner, [2 0], [i+10], [true])
        end

        exactpruning!(pruner, ClpSolver())
        # perform the same test again
        @test pruner.b == [20]
    end
end