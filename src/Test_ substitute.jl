using RootedTrees
using BSeries
using Revise
using BSeries_Analysis.BSeries_Analysis
using LinearAlgebra: dot
using Counters
using RootedTrees_SubtreeStructures
using BenchmarkTools
using OrderedCollections: OrderedDict
using ProfileView
using PProf


n=9
(tree_list,order_list)=generateTrees_subtrees(n)
subtree_dict=Dict{Vector{Int},Int}()
exact_series=OrderedDict{Int,Rational{Int64}}()
exact_series[0]=1

memory_dict=Dict{Int,Vector{Tuple{Vector{Int64},Vector{Int64}}}}()
for x in tree_list
    subtree_dict[x.subtrees]=x.index
    exact_series[x.index]=1//RootedTrees_SubtreeStructures.density(x,tree_list)
    #print(a[x.index],"\n")
end
partition_list=RootedTrees_SubtreeStructures.substitution_partitions(tree_list,subtree_dict)
A=[0 0;1//2 0]
b=[0,1//1]
c=[0,1//2]
rk=BSeries_Analysis.RungeKuttaMethod(A,b,c)
series=OrderedDict{Int,Rational{Int64}}()
series[0]=1
for tree in tree_list
    series[tree.index]=BSeries_Analysis.elementary_weight(tree,tree_list,rk)
end
f(t,series)=1//RootedTrees.density(t)
io=IOContext(stdout,:histmin=>15000000,:logbins=>true)
bseries1=BSeries.bseries(A,b,c,n)
bseries2=BSeries.bseries(f,n)
result=BSeries.substitute(bseries1,bseries2)
copy1=copy(series)
copy2=copy(exact_series)
copy3=copy(series)
copy4=copy(exact_series)
copy5=copy(series)
copy6=copy(exact_series)
println("BSeries    Result:",sum(values(result)),"\n")
#b=@benchmark BSeries.substitute(BSeries.bseries(A,b,c,n),BSeries.bseries(f,n))
@time BSeries.substitute(bseries1,bseries2)
#show(io,MIME("text/plain"),b)
println("\n")
result=BSeries_Analysis.substitute(copy1,copy2,partition_list)
println("Substitute     Result:",sum(values(result)),"\n")
#b=@benchmark BSeries_Analysis.substitute(copy(series),copy(exact_series),partition_list)
@time BSeries_Analysis.substitute(copy1,copy2,partition_list)
#show(io,MIME("text/plain"),b)
println("\n")
result=BSeries_Analysis.substitute_threads(copy3,copy4,partition_list)
println("Substitute_threads     Result:",sum(values(result)),"\n")
@time BSeries_Analysis.substitute_threads(copy3,copy4,partition_list)
#b=@benchmark BSeries_Analysis.substitute_threads(copy(series),copy(exact_series),partition_list)
#show(io,MIME("text/plain"),b)
result=BSeries_Analysis.substitute_threads_spawn(copy5,copy6,partition_list)
println("Substitute_spawn    Result:",sum(values(result)),"\n")
@time BSeries_Analysis.substitute_threads_spawn(copy5,copy6,partition_list)
print("\n\n\n\n\n")
