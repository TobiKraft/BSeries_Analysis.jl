using RootedTrees
using BSeries
using Revise
using BSeries_Analysis.BSeries_Analysis
using LinearAlgebra: dot
using Counters
using RootedTrees_SubtreeStructures
using BenchmarkTools
using OrderedCollections: OrderedDict
#add https://github.com/TobiKraft/RootedTrees_SubtreeStructures.jl

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



"""
println("Result:")
for y in range(1,10)
    println(y,":",result[y])
    
end
println("Series:")
for (x,y) in series
    if series[x]!=0
    println(x,":",y)
    end
    
end
println("Exact:")
for (x,y) in exact_series
    if exact_series[x]!=0
    println(x,":",y)
    end
    
end"""

function Partition_test(level_seq)
    for x in level_seq
        y=RootedTrees.PartitionIterator(x)
        for a in y
            b=a
        end
    end
end


print("C:_-----------------\n")
#for x in keys(c)
#    print("Key:",x,"    Wert:",c[x],"\n")
#end
#print("Compose:")
#@time c=BSeries_Analysis.compose(a,b,tree_list,subtree_dict)
#print("\n Memory:")
#@time d=BSeries_Analysis.compose_mem(a,a,tree_list,subtree_dict,memory_dict)
#data=generateTrees_butcher(n)

#print("\n BSeries:")
#@time c=BSeries.compose(b,b)
#for i in range(1,length(tree_list))
#     d[tree_list[i].index]!=c[level_seq[i]] ? error(d[tree_list[i].index],"  ", c[level_seq[i]]," ","false at index:",i) : continue
#end
#print(true,"\n")
function subtree_test(tree_list)
    
    for i in range(1,length(tree_list))
        #print(i,"\n")
        b=RootedTrees_SubtreeStructures.orderedSubtrees_and_Forests(tree_list[i],tree_list,subtree_dict)
        for x in b
            a=x
        end
    end
end
function butcher_test(data)
    for i in range(1,length(data.tree_list))
        #print(i,"\n")
        b=RootedTrees_SubtreeStructures.orderedSubtrees_and_Forests(i,data)
        for x in b
            a=x
        end
    end
end
function rootedTree_test(rootedTree_list)
    for x in rootedTree_list
        b=SplittingIterator(x)
        for x in b
            a=x
        end
    end
end
function substitute_test(tree_list,subtree_dict)
    list=RootedTrees_SubtreeStructures.substitution_partitions(tree_list,subtree_dict)
    for x in list
        a=x
        #print(x,"\n")
    end   
    
end

function Substitution_Testing(a,b,partition_list) #tree_list,subtree_dict
    c=BSeries_Analysis.substitute(a,b,partition_list)#tree_list,subtree_dict
    return c
end

function BSeries_substitute_test(a)
    c=BSeries.substitute(a,a)
    return c
end
function PartitionIterator_test(level_seq)
    result=RootedTrees.PartitionIterator(level_seq)
    for (forest,skeleton) in result
        print("\n Skeleton:",skeleton,"\n")
        for x in forest
            print(x,"   ")
        end
    end

end

    """
    partition_list=RootedTrees_SubtreeStructures.substitution_partitions(tree_list,subtree_dict)
    partition_list[0]=[]
    print("Meins")
    
    print("\nSeins:")
    @time BSeries_substitute_test(b)
    @time substitute_test(tree_list,subtree_dict)
    print("Done")
    print("RootedTree")
    @time rootedTree_test(level_seq)
    print("\n Selbst:\n")
    i=1
    for (skeleton,forest) in liste[m]
        read(stdin, 1); nothing
        print("\n i:",i,"   Skeleton:",skeleton,"     Forest:\n")
        for x in forest
            print(x,"     ")
                
        end
        global i+=1
    
    end
    io = IOContext(stdout, :histmin=>20000000, :histmax=>300000000)#, :logbins=>true)
    d=@benchmark Substitution_Testing(a,a,partition_liste)
    show(io,MIME("text/plain"),d)
    #readline()


    d=@benchmark BSeries.substitute(b,b)
    print("\nProf:")
    show(io,MIME("text/plain"),d)


    for (i,x) in enumerate(tree_list)
        print("\n",i,":\n")
        rooted_tree=RootedTrees.rootedtree(RootedTrees_SubtreeStructures.Subtrees_to_Levelsequence(tree_list[x.index],tree_list))
        iter=RootedTrees.PartitionIterator(rooted_tree)
        for (forest,skeleton) in iter
            print("Skeleton:",skeleton,"\n")
            print("Forest:")
            for y in forest
                print(y,"  ")
            end
            print("\n")
        end
    end


    for i in range(1,length(tree_list))
    print("\n",i,":\n")
    x=partitions_transfer(tree_list[i],tree_list,subtree_dict)
    for (a,b) in x
        print("Skeleton:",a,"\n Forest:")
        for t in b
            print(t,"   ")
        end
        print("\n")
    end
end
    """

m=10
#PartitionIterator_test(level_seq[m])


 #tree_list,subtree_dict

#print("\n")


#print("ende")
function partitions_transfer(tree,tree_list,subtree_dict)
    if tree.index==1
        return [(1,[1])]
    end
    t1=RootedTrees
    t2=RootedTrees
    answer=Vector{Tuple{Int,Vector{Int}}}()
    right=last(tree.subtrees)   
    left=subtree_dict[tree.subtrees[1:end-1]]
    right_transfered=rootedtree(Subtrees_to_Levelsequence(tree_list[right],tree_list))
    left_transfered=rootedtree(Subtrees_to_Levelsequence(tree_list[left],tree_list))
    iterator_left=collect(RootedTrees.PartitionIterator(left_transfered))
    iterator_right=collect(RootedTrees.PartitionIterator(right_transfered))
    for  (left_partition_forest,left_partition_skeleton) in iterator_left, (right_partition_forest,right_partition_skeleton) in iterator_right
        left_skeleton=RootedTrees_SubtreeStructures.Levelsequence_to_Subtrees(left_partition_skeleton.level_sequence,tree_list,subtree_dict)
        right_skeleton=RootedTrees_SubtreeStructures.Levelsequence_to_Subtrees(right_partition_skeleton.level_sequence,tree_list,subtree_dict)
        left_forest=Vector{Int}()
        right_forest=Vector{Int}()
        for x in left_partition_forest
            append!(left_forest,RootedTrees_SubtreeStructures.Levelsequence_to_Subtrees(x.level_sequence,tree_list,subtree_dict).index)
        end
        for x in right_partition_forest
            append!(right_forest,RootedTrees_SubtreeStructures.Levelsequence_to_Subtrees(x.level_sequence,tree_list,subtree_dict).index)
        end
        circ_subtrees=copy(left_skeleton.subtrees)
        append!(circ_subtrees,right_skeleton.index)
        circ_product_skeleton=subtree_dict[sort(circ_subtrees,rev=true)]
        circ_product_forest=[]
        append!(circ_product_forest,left_forest,right_forest)
        merge_subtrees=copy(left_skeleton.subtrees)
        append!(merge_subtrees,right_skeleton.subtrees)
        merge_product_skeleton=subtree_dict[sort(merge_subtrees,rev=true)]
        merge_forest_subtrees=copy(tree_list[left_forest[1]].subtrees)
        append!(merge_forest_subtrees,right_forest[1])
        merge_forest=[subtree_dict[sort(merge_forest_subtrees,rev=true)]]
        append!(merge_forest,left_forest[2:end],right_forest[2:end])
        append!(answer,[(circ_product_skeleton,circ_product_forest)])
        append!(answer,[(merge_product_skeleton,merge_forest)])
    end
    return answer
end
function transfer_test(tree_list,subtree_dict)
    for i in range(1,length(tree_list))
        x=partitions_transfer(tree_list[i],tree_list,subtree_dict)
        for (a,b) in x
            f=a
            g=b
        end
    end
end
function substitute_transfer_test(a,b,tree_list,subtree_dict)
    series_keys=keys(b)
    series=empty(b)
    t=first(series_keys)
    @assert iszero(t)
    series[t]=a[t]
    for t in Iterators.drop(series_keys,1)
        coefficient=BSeries_Analysis.substitute(a,b,partitions_transfer(tree_list[t],tree_list,subtree_dict))
        series[t]=coefficient
    end
    return series    
end


#c=BSeries_Analysis.compose(a,b,tree_list,subtree_dict)
level_seq=[RootedTrees.rootedtree(Subtrees_to_Levelsequence(x,tree_list)) for x in tree_list]
f(x,series)=1//RootedTrees.density(x)
#typeof(f)
b=BSeries.bseries(f,n)

print("BSeries:\n")
for x in range(1,5)
    print("BSeries:\n")
    @time BSeries.substitute(b,b)
    print("BSeries_Analysis:\n")
    #@time Substitution_Testing(a,a,partition_list)
    print("\n\n\n")
end 
#print("\n BSeries_Analysis Transfer")
#@time substitute_transfer_test(a,a,tree_list,subtree_dict)
#print("\n BSeries_Analysis Substitute")

#print("BSeries_Analysis:\n")
#for x in range(1,10)
#    @time Substitution_Testing(a,a,partition_list)
#    print("\n")
#end
#print("Ende")


partition_list=RootedTrees_SubtreeStructures.substitution_partitions(tree_list,subtree_dict)
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
    end
using BSeries, StaticArrays
function BSeries_Modified_eq(A,b,c,up_to_order)
    series = BSeries.modified_equation(A, b, c, up_to_order)
end
#A = @SArray [0 0; 1//2 0]
#b = @SArray [0, 1//1]
#c = @SArray [0, 1//2]
#up_to_order = 9
#println("BSeries:Modified equation")
#io=IOContext(stdout,:histmin=>40000000,:logbins=>true)
#b=@benchmark BSeries_Modified_eq(A,b,c,up_to_order)
#BenchmarkTools.save("BSeries_result.json",b)
#show(io,MIME("text/plain"),b)
#println()
#println("BSeries_Analysis:Modified equation")
#b=@benchmark modified_equation(tree_list,partition_list)
#show(io,MIME("text/plain"),b)

