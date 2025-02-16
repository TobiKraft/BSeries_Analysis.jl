using RootedTrees
using BSeries
using Revise
using BSeries_Analysis.BSeries_Analysis
using LinearAlgebra: dot
using Counters
using RootedTrees_SubtreeStructures
using BenchmarkTools
tree_list=generateTrees_subtrees(6)[1]
tree_list2=generateTrees_butcher(6)
a=RootedTrees_SubtreeStructures.orderedSubtrees_and_Forests(tree_list[2],tree_list)
display(a)
series=TruncatedBSeries{Int,Int}()