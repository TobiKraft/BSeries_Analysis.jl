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
series=OrderedDict{Union{Int64,Vector{Int64}},Rational}()
series[Int64[]]=1

for t in tree_list
    series[t.index]=exact_value(t,tree_list)
    print(series[t.index],"\n")
end
series2=series
series3=BSeries_Analysis.compose(series,series,tree_list)
for x in series3
    print(x,"\n")
end