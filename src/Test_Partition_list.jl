using RootedTrees
using BSeries
using Revise
using BSeries_Analysis
using LinearAlgebra: dot
using Counters
using RootedTrees_SubtreeStructures
using BenchmarkTools
using OrderedCollections: OrderedDict
using ProfileView
using PProf
function readablesize(x)
    Base.format_bytes(Base.summarysize(x))
end
i=10
size_storage=Vector{Vector{String}}()
for n in range(1,i)
println(n,":")
tree_list=[]
order_list=[]
partition_list=[]
partition_list2=[]
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
print("\n BSeries:")
result= @time BSeries.substitute(bseries1,bseries2)
print("\n Partition_Normal:")
partition_list=@time RootedTrees_SubtreeStructures.substitution_partitions(tree_list,subtree_dict)
print("\n BSeries_Analysis_Partition_Normal ")
result2=@time BSeries_Analysis.substitute(copy1,copy2,partition_list)
print("\n Partition_Lower_Memory:")
size_1=readablesize(partition_list)
partition_list=@time RootedTrees_SubtreeStructures.substitution_partitions_lower_memory(tree_list,subtree_dict)
print("\n BSeries_Analysis_Partition_Memory_optimized: ")
result3=@time BSeries_Analysis.substitute(copy3,copy4,partition_list,tree_list)
print("\n Partition_Lower_Memory2:")
size_2=readablesize(partition_list)
partition_list=@time RootedTrees_SubtreeStructures.substitution_partitions_lower_memory2(tree_list,subtree_dict)
print("\n BSeries_Analysis_Partition_Memory_optimized2: ")
result4=@time BSeries_Analysis.substitute(copy5,copy6,partition_list,tree_list)
size_3=readablesize(partition_list)
push!(size_storage,[size_1,size_2,size_3])
println("Normal: ",size_1,"  Optimized:", size_2,"  Optimized2:", size_3)

println("-----------------------------------")



println("\n Result:"  ,sum(values(result)),"   ",sum(values(result2)),"   ",sum(values(result3)), "   ",sum(values(result4))) 
println("---------------------------------------------------------------------------------------------------------------------------------------")

end
display(size_storage)
open("storage_memory_data_for_diagramm6.txt","w") do h
    for (x,y) in size_storage
        println(h,"Normal:  ",x,"     Optimized:  ",y,"    Ratio:", "\n") #,round(x//y,digits=3)
    end
end