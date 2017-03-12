@testset "Avg cut pruning" begin
    algo = AvgCutPruningAlgo(2)
    pruner = CutPruner{2, Int}(algo, :≤)
    @test 1:1 == addcuts!(pruner, [1 0], [1], [true])
    @test 2:2 == addcuts!(pruner, [0 1], [1], [true])
    CutPruners.updatestats!(pruner, [1, 0])
    @test [2, 0, 0] == addcuts!(pruner, [1 1; -1 -1; 0 1], [1, 1, 2], [true, false, true])
    @test pruner.A == [1 0; 1 1]
    @test pruner.b == [1, 1]
    @test pruner.ids == [1, 3]
end
