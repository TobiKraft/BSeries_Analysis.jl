module BSeries_Analysis
using OrderedCollections
using Revise
using Counters
using LinearAlgebra: dot
using RootedTrees_SubtreeStructures
export AbstractTimeIntegrationMethod,RungeKuttaMethod
export TruncatedBSeries
export exact_value
export compose,compose_mem
export bseries

#--------------------------------------------------------------------------------------
#           Export RootedTrees_given_by_ButcherProduct
#--------------------------------------------------------------------------------------
export exact_value_butcher, exact_value_butcher!
export exact_series_butcher, exact_series_butcher!
include("BSeries_RootedTrees_given_by_subtrees.jl")
include("BSeries_RootedTrees_given_by_ButcherProduct.jl")
include("BSeries_CircProduct.jl")
#Github RootedTrees_SubtreeStructures
#add https://github.com/TobiKraft/RootedTrees_SubtreeStructures.jl
"""
abstract type AbstractTimeIntegrationMethod end


struct RungeKuttaMethod{T, MatT <: AbstractMatrix{T}, VecT <: AbstractVector{T}} <:
    AbstractTimeIntegrationMethod
 A::MatT
 b::VecT
 c::VecT
end
Base.eltype(rk::RungeKuttaMethod{T}) where {T} = T
function RungeKuttaMethod(A::AbstractMatrix, b::AbstractVector,
                       c::AbstractVector = vec(sum(A, dims = 2)))
    T = promote_type(eltype(A), eltype(b), eltype(c))
    _A = T.(A)
    _b = T.(b)
    _c = T.(c)
    return RungeKuttaMethod(_A, _b, _c)
end

function elementary_weight(tree::Union{Int,RootedTree_given_by_subtrees},tree_list::Union{Vector{RootedTree_given_by_subtrees},Data_given_by_ButcherProduct},rk::RungeKuttaMethod)
    dot(rk.b, derivative_weight(tree,tree_list,rk))
end

function derivative_weight(tree::RootedTree_given_by_subtrees,tree_list::Vector{RootedTree_given_by_subtrees},rk::RungeKuttaMethod) 
    A = rk.A
    c = rk.c

    # Initialize `result` with the identity element of pointwise multiplication `.*`
    result = zero(c) .+ one(eltype(c))

    # Iterate over all subtrees and update the `result` using recursion
    count= counter(tree.subtrees)
    for (subtree,quantity) in zip(keys(count),values(count))
        tmp = A * derivative_weight(tree_list[subtree],tree_list,rk)
        result = result .* tmp.^(quantity)
    end

    return result
end
function derivative_weight(index::Int,tree_list::Data_given_by_ButcherProduct,rk::RungeKuttaMethod) 
    A = rk.A
    c = rk.c

    # Initialize `result` with the identity element of pointwise multiplication `.*`
    result = zero(c) .+ one(eltype(c))

    # Iterate over all subtrees and update the `result` using recursion
    counter=children(index,tree_list)
    for (subtree,quantity) in zip(keys(counter),values(counter))
        tmp = A * derivative_weight(subtree,tree_list,rk)
        result = result .* tmp.^(quantity)
    end

    return result
end
"""
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                           BSeries
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
function bseries_ohne_order(rk::RungeKuttaMethod,order::Int,tree_list::Vector{RootedTree_given_by_subtrees}) 
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
    while state<=end_value&&(tree.order<=order)
        series[tree.index]=elementary_weight(tree,tree_list,rk)
        (tree,state)=iterate(tree_list,state)
    end
    #for tree in tree_list[1:(order_list[order+1]-1)]
    #   series[tree.index]=elementary_weight(tree,tree_list,rk)
    #end
    if state>end_value
        series[tree.index]=elementary_weight(tree,tree_list,rk)
    end
    return series
end
function bseries_mit_order(rk::RungeKuttaMethod,order::Int,tree_list::Vector{RootedTree_given_by_subtrees},order_list) 
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
    #while tree.order<=order
    #    series[tree.index]=elementary_weight(tree,tree_list,rk)
    #    (tree,state)=iterate(tree_list,state)
    #end
    for tree in tree_list[1:(order_list[order+1]-1)]
        series[tree.index]=elementary_weight(tree,tree_list,rk)
    end
    return series
end
function bseries(A::AbstractMatrix,b::AbstractVector,c::AbstractVector,order,tree_list)
    rk=RungeKuttaMethod(A,b,c)
    return bseries(rk,order,tree_list)
end
function bseries(f::Function,order,tree_list,order_list)
    series=OrderedDict{Int,V}()
    t=0
    v=f(t,nothing)
    V_tmp=typeof(v)
    if V_tmp<:Integer
        V = Rational{V_tmp}
    else
        V = V_tmp
    end
    series[t]=v
    (tree,state)=iterate(tree_list)
    while tree.order<=order
        series[tree.index]=elementary_weight(tree,tree_list,rk)
        (tree,state)=iterate(tree_list,state)
    end
end




#function exact_value(tree::RootedTree_given_by_subtrees,tree_list::Array{RootedTree_given_by_subtrees})
#    return 1//density(tree,tree_list)
#end

#function exact_value(index::Int,data::Data_given_by_ButcherProduct)
#    return 1//density(index,data)
#end

function compose(a,b,t::Vector{},tree_list::Array{RootedTree_given_by_subtrees})
    return a[Int64[]]*b[Int64[]]
end
function compose(a,b,t::Vector{},data::Data_given_by_ButcherProduct)
    return a[Int64[]]*b[Int64[]]
end

function compose(a,b,t::RootedTree_given_by_subtrees,tree_list::Array{RootedTree_given_by_subtrees},subtree_dict)
    result = zero(first(values(a)) * first(values(b)))
    vector_subtrees_forests=orderedSubtrees_and_Forests(t,tree_list,subtree_dict)
    #catch the empty tree 
    #result+=b[Int64[]]*a[t.index]
    for (subtree,forest) in vector_subtrees_forests
        #subtree is a vector containing one value
        tmp=b[subtree[1]]
        for tree in forest 
            tmp*=a[tree]
        end
        result+=tmp
    end
    return result
end


function compose(a,b,tree_list::Array{RootedTree_given_by_subtrees},subtree_dict; normalize_stepsize=false)
    c=copy(a)
    series_keys= keys(delete!(c,0))
    series = Dict{Int,Rational}()
    series[0]=a[0]*b[0]
    for t in series_keys      
            coefficient = compose(a,b,tree_list[t], tree_list,subtree_dict)
        if normalize_stepsize
            coefficient /= 2^order(t)
        end
        series[t] = coefficient
    end

    return series
end

function compose_mem(a,b,t::RootedTree_given_by_subtrees,tree_list::Array{RootedTree_given_by_subtrees},subtree_dict,memory_dict)
    result = zero(first(values(a)) * first(values(b)))
    (vector_subtrees_forests,memory_dict_new)=mem_orderedSubtrees_and_Forests(t,tree_list,subtree_dict,memory_dict)
    #catch the empty tree 
    #result+=b[Int64[]]*a[t.index]
    for (subtree,forest) in vector_subtrees_forests
        #subtree is a vector containing one value
        tmp=b[subtree[1]]
        for tree in forest 
            tmp*=a[tree]
        end
        result+=tmp
    end
    return (result,memory_dict_new)
end

function compose_mem(a,b,tree_list::Array{RootedTree_given_by_subtrees},subtree_dict,memory_dict; normalize_stepsize=false)
    c=copy(a)
    series_keys= keys(delete!(c,0))
    series = Dict{Int,Rational}()
    series[0]=a[0]*b[0]
    memory_dict_new=memory_dict
    for t in series_keys      
            (coefficient,memory_dict_new) = compose_mem(a,b,tree_list[t], tree_list,subtree_dict,memory_dict_new)
        if normalize_stepsize
            coefficient /= 2^order(t)
        end
        series[t] = coefficient
    end

    return series
end

function compose(a,b,index::Int,data::Data_given_by_CircProduct)
    result = zero(first(values(a)) * first(values(b)))
    vector_subtrees_forests=orderedSubtrees_and_Forests(index,data)
    #catch the empty tree 
    result+=b[Int64[]]*a[index]
    for (subtree,forest) in vector_subtrees_forests[2:end]
        tmp=b[subtree[1]]
        for tree in forest 
            tmp*=a[tree]
        end
        result+=tmp
    end
    return result
end

function compose(a,b,data::Data_given_by_CircProduct; normalize_stepsize=false)
   
    series_keys = keys(a)
    series = empty(a)
    for t in series_keys     
            coefficient = compose(a,b,t,data)
        if normalize_stepsize
            coefficient /= 2^order(t)
        end
        series[t] = coefficient
    end

    return series
end
"""
function substitute(b,a,partitions::Vector{Tuple{Int,Vector{Int}}})
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
function substitute(b,a,tree_list,subtree_dict)
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
function substitute(b,a,partition_list::OrderedDict{Int,Vector{Tuple{Int,Vector{Int}}}})
    series_keys=keys(b)
    series=empty(b)
    t=first(series_keys)
    @assert iszero(t)
    series[t]=a[t]
    for t in Iterators.drop(series_keys,1)
        coefficient=substitute(b,a,partition_list[t])
        series[t]=coefficient
    end
    return series
end
"""
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                           Modified Equations
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
end