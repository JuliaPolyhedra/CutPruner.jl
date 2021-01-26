using CutPruners

include("solvers.jl")

isempty(lp_solvers) && warn("Exact Pruning tests not run!")
@testset "Exact Pruning with $solver" for solver in lp_solvers
    optimizer_constructor = optimizer_with_attributes(solver.Optimizer, MOI.Silent() => true)
    for algo in [AvgCutPruningAlgo(20), DecayCutPruningAlgo(20), DeMatosPruningAlgo(20)]
        @testset "Exact pruning" begin
            pruner = CutPruner{2, Int}(algo, :Max)

            # test pruning with one cut
            addcuts!(pruner, [1 0], [0], [true])
            exactpruning!(pruner, optimizer_constructor)
            @test pruner.b == [0]

            # add 10 cuts in a row
            for i in 1:10
                addcuts!(pruner, [1 0], [i], [true])
            end

            exactpruning!(pruner, optimizer_constructor)
            # normally the exact pruning saves only the last cuts
            @test pruner.b == [10]
            # ... and add another set of cuts ...
            for i in 1:10
                addcuts!(pruner, [2 0], [i+10], [true])
            end

            exactpruning!(pruner, optimizer_constructor)
            # perform the same test again
            @test pruner.b == [10, 20]

            pruner = CutPruner{1, Int}(algo, :Max)

            addcuts!(pruner, [1]', [0], [true])
            addcuts!(pruner, [-1]', [1], [true])
            exactpruning!(pruner, optimizer_constructor, lb = 0.5)
            @test pruner.b == [0]

            addcuts!(pruner, [-1]', [1], [true])
            exactpruning!(pruner, optimizer_constructor, lb = [0.5])
            @test pruner.b == [0]

            addcuts!(pruner, [-1]', [1], [true])
            exactpruning!(pruner, optimizer_constructor, ub = 0.4)
            @test pruner.b == [1]

            addcuts!(pruner, [1]', [0], [true])
            exactpruning!(pruner, optimizer_constructor, ub = [0.4])
            @test pruner.b == [1]

        end
    end
    for algo in [AvgCutPruningAlgo(20), DecayCutPruningAlgo(20), DeMatosPruningAlgo(20)]
        @testset "Exact pruning" begin
            pruner = CutPruner{2, Int}(algo, :Min)
            # add 10 cuts in a row
            for i in 1:10
                addcuts!(pruner, [1 0], [i], [true])
            end

            exactpruning!(pruner, optimizer_constructor)
            # normally the exact pruning saves only the last cuts
            @test pruner.b == [1]
            # ... and add another set of cuts ...
            for i in 1:10
                addcuts!(pruner, [2 0], [i+10], [true])
            end

            exactpruning!(pruner, optimizer_constructor)
            # perform the same test again
            @test pruner.b == [1, 11]
        end
    end
end
