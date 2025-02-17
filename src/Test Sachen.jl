using RootedTrees
using BSeries
using Revise
using BSeries_Analysis.BSeries_Analysis
using LinearAlgebra: dot
using Counters
using RootedTrees_SubtreeStructures
using BenchmarkTools
tree_list=generateTrees_subtrees(6)[1]
data=generateTrees_butcher(6)
series_butcher=BSeries_Analysis.TruncatedBSeries{Union{Int64,Vector{Int64}},Rational}()
series_subtrees=BSeries_Analysis.TruncatedBSeries{Union{Int64,Vector{Int64}},Rational}()
series_butcher[Int64[]]=1
series_subtrees[Int64[]]=1
len=length(data.tree_list)

for t in range(1,len)
    series_butcher[t]=exact_value(t,data)
end
for tree in tree_list
    series_subtrees[tree.index]=exact_value(tree,tree_list)
end

for tree in tree_list
    subtree_value=series_subtrees[tree.index]
    butcher_index=Subtrees_to_ButcherProduct(tree,tree_list,data)
    butcher_value=series_butcher[butcher_index]
    bool=(subtree_value==butcher_value)
    if !bool
     print("Subtree:", subtree_value, "     Butcher:",butcher_value,"   ", bool,"\n")
    end
end
series_butcher2=BSeries_Analysis.compose(series_butcher,series_butcher,data)
series_subtrees2=BSeries_Analysis.compose(series_subtrees,series_subtrees,tree_list)
print("Index:Int64[]    Subtree:", series_subtrees2[Int64[]], "     Butcher:",series_butcher2[Int64[]],"\n")
for tree in tree_list
    subtree_value=series_subtrees2[tree.index]
    butcher_index=Subtrees_to_ButcherProduct(tree,tree_list,data)
    butcher_value=series_butcher2[butcher_index]
    bool=(subtree_value==butcher_value)
    print("Index:",tree.index,"   Subtree:", subtree_value, "     Butcher:",butcher_value,"   ", bool,"\n")
end