# Important that oldest is removed first !
# Example showing why this is important:
# For instance, if the lattice is:
#            :Unbounded - ...
#           /
# :Unbounded
#           \
#            :Optimal - ...
# The lower node will give optimality cuts and fill the cutpruner of the root node
# Meanwhile, the upper node will not give any cut so the root not will stay unbounded
# Since the root node is :Unbounded it will not have dual variables for the stats, so
# it may happen (for instance with AvgCutPruner) that all the cuts are at tie !
# When the upper node stops being :Unbounded and gives its cut that is necessary for
# the root node to stop being :Unbounded. However, if the lower node replaces it with
# another optimality cut it will be ignored !
for algo in [AvgCutPruningAlgo(3), DecayCutPruningAlgo(3)]
    @testset "Oldest removed first for $(typeof(algo))" begin
        pruner = CutPruner{2, Int}(algo)
        addcuts!(pruner, [1 0], [0], false, [true])
        addcuts!(pruner, [2 0], [0], false, [true])
        addcuts!(pruner, [3 0], [0], false, [true])
        addcuts!(pruner, [4 0], [0], false, [true])
        @test sort(pruner.cuts_DE[:,1]) == [2, 3, 4]
        addcuts!(pruner, [5 0], [0], false, [true])
        @test sort(pruner.cuts_DE[:,1]) == [3, 4, 5]
        # However, order should only be used to break ties
        addcuts!(pruner, [6 0], [0], false, [false])
        @test sort(pruner.cuts_DE[:,1]) == [3, 4, 5]
        addcuts!(pruner, [7 0; 8 0], [0, 0], false, [true, true])
        @test sort(pruner.cuts_DE[:,1]) == [5, 7, 8]
    end
end
