using Revise
using RootedTrees_SubtreeStructures
using RootedTrees
using OrderedCollections
using Counters
using LinearAlgebra: dot

export exact_value
export bseries
export substitute, substitute_threads,substitute_threads_spawn
export modified_equation
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                           BSeries
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
    exact_value(tree::Union{Int,RootedTree_given_by_subtrees},tree_list::Array{RootedTree_given_by_subtrees}) -> Rational
``tree`` can be a ``RootedTree_given_by_subtrees`` or an Integer representing the index of ``tree`` in ``tree_list``.\\
For a given autonomous problem
* y'=f(y),        y(t_0)=y_0

and a tree tau, this function returns the coefficient alpha(tau) as a Rational of the BSeries, that is representing the exact solution of the given problem.
# Example
```julia
julia>(tree_List,order_list)=generateTrees_subtrees(5)
julia>exact_value(6,tree_list)
1//8
julia>exact_value(tree_list[6],tree_list)
1//8
```       
"""
function exact_value(tree::Union{Int,RootedTree_given_by_subtrees},tree_list::Array{RootedTree_given_by_subtrees})
    return 1//RootedTrees_SubtreeStructures.density(tree,tree_list)
end

"""
    bseries(rk::RungeKuttaMethod,order::Int,tree_list::Vector{RootedTree_given_by_subtrees}) -> OrderedDict{(Int,V)} where {V<:Number}
Returns the coefficients of the BSeries representing the given Runge-Kutta_method ``rk`` for every tree up to the given ``order``.\\
The coefficients are stored as an OrderedDictionary, using the index of a tree as key and its coefficient as its value.\\
The type of the value is either equal to the type of the Runge-Kutta-Metho or is Rational if the type of the Runge-Kutta-Method is a subtype of Integer.
# Example
```julia
julia>A=[0 0;1//2 0]; b=[0,1]; c=[0,1//2]
julia>rk=RungeKuttaMethod(A,b,c)
julia>bseries(rk,3,tree_list)
OrderedDict{Int64, Rational{Int64}} with 5 entries:
  0 => 1
  1 => 1
  2 => 1//2
  3 => 1//4
  4 => 0
```
"""
function bseries(rk::RungeKuttaMethod,order::Int,tree_list::Vector{RootedTree_given_by_subtrees}) 
    V_tmp=eltype(rk)
    if V_tmp<:Integer
        V=Rational{V_tmp}
    else
        V=V_tmp
    end
    series=OrderedDict{Int,V}()
    series[0]=one(V)
    iter=iterate(tree_list)
    (tree,state)=iter
    end_value=length(tree_list)
    #Iterating about tree_list until either
    #the end of the tree_list is reached
    #or the order of the current tree is higher than the given order.
    while state<=end_value&&(tree.order<=order)
        series[tree.index]=elementary_weight(tree,tree_list,rk)
        (tree,state)=iterate(tree_list,state)
    end
    if state>end_value
        series[tree.index]=elementary_weight(tree,tree_list,rk)
    end
    return series
end

function bseries(rk::RungeKuttaMethod,order::Int,tree_list::Vector{RootedTree_given_by_subtrees},order_list) 
    V_tmp=eltype(rk)
    if V_tmp<:Integer
        V=Rational{V_tmp}
    else
        V=V_tmp
    end
    series=OrderedDict{Int,V}()
    series[0]=one(V)
    iter=iterate(tree_list)
    (tree,state)=iter
    #Iterating about the tree_list, up to the last tree with the given order by looking up its index in order_list
    for tree in tree_list[1:(order_list[order+1]-1)]
        series[tree.index]=elementary_weight(tree,tree_list,rk)
    end
    return series
end
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                           Substitute
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function substitute(b::OrderedDict{Int,V},a::OrderedDict{Int,V},partitions::Vector{Tuple{Int,Vector{Int}}}) where {V<:Number}
    result=zero(first(values(a))*first(values(b)))
    for (skeleton,forest) in partitions
        update=a[skeleton]
        for tree in forest
            update*=b[tree]
        end
        result+=update
    end
    return result
end
function substitute(b::OrderedDict{Int,V},a::OrderedDict{Int,V},partitions::Vector{Tuple{Int,Tuple{Int,Int}}},tree_list) where {V<:Number}
    result=zero(first(values(a))*first(values(b)))
    for (skeleton,forest) in partitions
        update=a[skeleton]
        for tree in tree_list[forest[2]].subtrees
            update*=b[tree]
        end
        update*=b[forest[1]]
        result+=update
    end
    return result
end
function substitute(b::OrderedDict{Int,V},a::OrderedDict{Int,V},tree_list::Vector{RootedTree_given_by_subtrees}
                    ,subtree_dict::Dict{Vector{Int},Int}) where {V<:Number}
    series_keys=keys(b)
    series=empty(b)
    t=first(series_keys)
    @assert iszero(t)
    
    series[t]=a[t]
    partition_list=RootedTrees_SubtreeStructures.substitution_partitions(tree_list,subtree_dict)
    partition_list[0]=[]
    for t in Iterators.drop(series_keys,1)
        coefficient=substitute(b,a,partition_list[t])
        series[t]=coefficient
    end
    return series
end
function substitute(b::OrderedDict{Int,V},a::OrderedDict{Int,V},partition_list::OrderedDict{Int,Vector{Tuple{Int,Vector{Int}}}}) where {V<:Number}
    series_keys=keys(b)
    series=empty(b)
    t=first(series_keys)
    @assert iszero(t)
    series[t]=a[t]
    for t in Iterators.drop(series_keys,1)
        series[t]=substitute(b,a,partition_list[t])
        
    end
    return series
end
function substitute(b::OrderedDict{Int,V},a::OrderedDict{Int,V},partition_list::OrderedDict{Int,Vector{Tuple{Int,Tuple{Int,Int}}}},tree_list) where {V<:Number}
    series_keys=keys(b)
    series=empty(b)
    t=first(series_keys)
    @assert iszero(t)
    series[t]=a[t]
    for t in Iterators.drop(series_keys,1)
        series[t]=substitute(b,a,partition_list[t],tree_list)
        
    end
    return series
end
function substitute_threads(b::OrderedDict{Int,V},a::OrderedDict{Int,V},partition_list::OrderedDict{Int,Vector{Tuple{Int,Vector{Int}}}}) where {V<:Number}
    series_keys::Base.KeySet{Int64, OrderedDict{Int64, V}}=keys(b)
    series::OrderedDict{Int,V}=copy(b)
    t=first(series_keys)
    @assert iszero(t)
    series[t]=a[t]
    #series_keys=Iterators.drop(series_keys,1)
    #No data race inside differnt Threads, as every key t is only used once 
    Threads.@threads for t in collect(Iterators.drop(series_keys,1))
        #println("Key:",t,"   Thread_Id:",Threads.threadid())
        series[t]=substitute(b,a,partition_list[t])       
    end
    return series
end
function substitute_threads_spawn(b::OrderedDict{Int,V},a::OrderedDict{Int,V},partition_list::OrderedDict{Int,Vector{Tuple{Int,Vector{Int}}}}) where {V<:Number}
    series_keys::Base.KeySet{Int64, OrderedDict{Int64, V}}=keys(b)
    series::OrderedDict{Int,V}=copy(b)
    t=first(series_keys)
    @assert iszero(t)
    series[t]=a[t]
    #series_keys=Iterators.drop(series_keys,1)
    #No data race inside differnt Threads, as every key t is only used once 
    Threads.@sync for t in collect(Iterators.drop(series_keys,1))
        #println("Key:",t,"   Thread_Id:",Threads.threadid())
        Threads.@spawn series[t]=substitute(b,a,partition_list[t])       
    end
    return series
end
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                           Modified Equations
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function modified_equation(series::OrderedDict{Int,V},tree_list::Vector{RootedTree_given_by_subtrees},partition_list::OrderedDict{Int64, Vector{Tuple{Int64, Vector{Int64}}}}) where {V<:Number}
    exact_series=empty(series)
    result=empty(series)
    exact_series[0]=one(V)
    series_keys=keys(series)
    result[0]=0
    for index in Iterators.drop(series_keys,1)
        exact_series[index]=convert(V,1//RootedTrees_SubtreeStructures.density(tree_list[index],tree_list))
        result[index]=zero(V)
    end
    iter=iterate(series_keys)
    if iter !==nothing
        t, t_state = iter
        if iszero(t)
            iter = iterate(series_keys, t_state)
            if iter !== nothing
                t, t_state = iter
            end
        end
        result[t] = series[t]
        iter = iterate(series_keys, t_state)
    end
    while iter !==nothing
        t, t_state = iter
        result[t] += series[t] - substitute(result, exact_series, partition_list[t])
        iter = iterate(series_keys, t_state)
    end
    return result
end
function modified_equation(rk::RungeKuttaMethod,tree_list,partition_list)
    series=series=OrderedDict{Int,Rational}()
    series[0]=1
    exact_series=OrderedDict{Int,Rational}()
    exact_series[0]=1
end
