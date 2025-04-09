using RootedTrees
using BSeries
using StaticArrays
using Revise
using BSeries_Analysis.BSeries_Analysis
using LinearAlgebra: dot
using Counters
using RootedTrees_SubtreeStructures
using BenchmarkTools
using OrderedCollections: OrderedDict
#add https://github.com/TobiKraft/RootedTrees_SubtreeStructures.jl
#readablesize(x) =  Base.format_bytes(Base.summarysize(x))

n=9
(tree_list,order_list)=generateTrees_subtrees(n)
subtree_dict=Dict{Vector{Int},Int}()
a=OrderedDict{Int,Rational}()
a[0]=1

memory_dict=Dict{Int,Vector{Tuple{Vector{Int64},Vector{Int64}}}}()
for x in tree_list
    subtree_dict[x.subtrees]=x.index
    a[x.index]=1//RootedTrees_SubtreeStructures.density(x,tree_list)
    #print(a[x.index],"\n")
end
partition_list=RootedTrees_SubtreeStructures.substitution_partitions(tree_list,subtree_dict)   
print(typeof(partition_list))
function modified_equation(tree_list,partition_list)
    A=[0 0;1//2 0]
    b=[0,1//1]
    c=[0,1//2]
    rk=BSeries_Analysis.RungeKuttaMethod(A,b,c)
    series=OrderedDict{Int,Rational}()
    series[0]=1
    exact_series=OrderedDict{Int,Rational}()
    exact_series[0]=1
    result=OrderedDict{Int,Rational}()
    series=OrderedDict{Int,Rational}()
    series[0]=1
    exact_series=OrderedDict{Int,Rational}()
    exact_series[0]=1//1
    result[0]=0
    for tree in tree_list
        series[tree.index]=BSeries_Analysis.elementary_weight(tree,tree_list,rk)
        exact_series[tree.index]=1//RootedTrees_SubtreeStructures.density(tree.index,tree_list)
        result[tree.index]=0//1
    end
    result[1]=series[1]
    for tree in tree_list            
        result[tree.index]+=series[tree.index]-BSeries_Analysis.substitute(result,exact_series,partition_list[tree.index])
    end
    return result
end    

function BSeries_Modified_eq(A,b,c,up_to_order)
    series = BSeries.modified_equation(A, b, c, up_to_order)
end
A = @SArray [0 0; 1//2 0]
b = @SArray [0, 1//1]
c = @SArray [0, 1//2]
up_to_order = 9
result=BSeries_Modified_eq(A,b,c,up_to_order)
println("BSeries:Modified equation     Result:",sum(values(result)),"\n")
io=IOContext(stdout,:histmin=>40000000,:logbins=>true)
b=@benchmark BSeries_Modified_eq(A,b,c,up_to_order)
BenchmarkTools.save("BSeries_result.json",b)
show(io,MIME("text/plain"),b)
println("\n\n\n")
result=modified_equation(tree_list,partition_list)
println("BSeries_Analysis:Modified equation,     Result:",sum(values(result)),"\n")
b=@benchmark modified_equation(tree_list,partition_list)
show(io,MIME("text/plain"),b)