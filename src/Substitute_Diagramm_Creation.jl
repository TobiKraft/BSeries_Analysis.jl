using RootedTrees
using BSeries
using Revise
using BSeries_Analysis.BSeries_Analysis
using LinearAlgebra: dot
using Counters
using RootedTrees_SubtreeStructures
using BenchmarkTools
using OrderedCollections: OrderedDict
s=11
l=11
BSeries_times=Vector{Tuple{Int,Float64}}()
Substitute_times=Vector{Tuple{Int,Float64}}()
Substitute_Threads_times=Vector{Tuple{Int,Float64}}()
Substitute_Threads_spawn_times=Vector{Tuple{Int,Float64}}()
for n in range(s,l)
    ( tree_list,order_list)=generateTrees_subtrees(n)
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
    bseries1=BSeries.bseries(A,b,c,n)
    bseries2=BSeries.bseries(f,n)
    a=copy(series)
    b=copy(exact_series)
    c=copy(series)
    d=copy(exact_series)
    e=copy(series)
    f=copy(exact_series)
    b1=@benchmark BSeries.substitute($bseries1,$bseries2)
    print(minimum(b1),"       ")
    push!(BSeries_times,(n,minimum(b1).time))
    b2=@benchmark BSeries_Analysis.substitute($a,$b,$partition_list)
    push!(Substitute_times,(n,minimum(b2).time))
    b3=@benchmark BSeries_Analysis.substitute_threads($c,$d,$partition_list)
    push!(Substitute_Threads_times,(n,minimum(b3).time))
    b4=@benchmark BSeries_Analysis.substitute_threads_spawn($e,$f,$partition_list)
    push!(Substitute_Threads_spawn_times,(n,minimum(b4).time))
    print("BSeries:",BSeries_times[n-s+1],"   Substitute:",Substitute_times[n-s+1],"   Threads:",Substitute_Threads_times[n-s+1],"   Spawn:",Substitute_Threads_spawn_times[n-s+1],"\n")
end