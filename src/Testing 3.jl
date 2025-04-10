using BSeries_Analysis
using RootedTrees_SubtreeStructures
using BenchmarkTools
i=4
data=set_order(i)
(tree_list,order_list)=generateTrees_subtrees(i)
subtree_dict=Dict{Vector{Int},Int}()
for x in tree_list
    subtree_dict[x.subtrees]=x.index
end
@time substitution_partitions_lower_memory2(tree_list,subtree_dict)
@time _lambda_sub(data)
print("Ende")