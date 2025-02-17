module BSeries_Analysis
using OrderedCollections
using Revise
using Counters
using LinearAlgebra: dot
using RootedTrees_SubtreeStructures
export TruncatedBSeries
export exact_value
export compose

abstract type AbstractTimeIntegrationMethod end


struct RungeKuttaMethod{T, MatT <: AbstractMatrix{T}, VecT <: AbstractVector{T}} <:
    AbstractTimeIntegrationMethod
 A::MatT
 b::VecT
 c::VecT
end

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


struct TruncatedBSeries{T,V}
    coef::OrderedDict{T,V}
    
end
TruncatedBSeries{T, V}() where {T, V} = TruncatedBSeries{T, V}(OrderedDict{T, V}())

# general interface methods of `AbstractDict` for `TruncatedBSeries`

@inline function Base.getindex(series::TruncatedBSeries, t::Int)
    getindex(series.coef, t)    
end
@inline function Base.getindex(series::TruncatedBSeries, t::Int64)
    getindex(series.coef, t)
    
end
@inline function Base.setindex!(series::TruncatedBSeries, val, t::Int)
    setindex!(series.coef, val, t)
end
@inline function Base.setindex!(series::TruncatedBSeries, val, t::Int64)
    setindex!(series.coef, val, t)
end



function exact_value(tree::RootedTree_given_by_subtrees,tree_list::Array{RootedTree_given_by_subtrees})
    return 1//density(tree,tree_list)
end

function exact_value(index::Int,data::Data_given_by_ButcherProduct)
    return 1//density(index,data)
end

function compose(a::TruncatedBSeries,b::TruncatedBSeries,t::Vector{},tree_list::Array{RootedTree_given_by_subtrees})
    return a[Int64[]]*b[Int64[]]
end
function compose(a::TruncatedBSeries,b::TruncatedBSeries,t::Vector{},data::Data_given_by_ButcherProduct)
    return a[Int64[]]*b[Int64[]]
end

function compose(a::TruncatedBSeries,b::TruncatedBSeries,t::RootedTree_given_by_subtrees,tree_list::Array{RootedTree_given_by_subtrees})
    result = zero(first(values(a)) * first(values(b)))
    vector_subtrees_forests=orderedSubtrees_and_Forests(t,tree_list)
    #catch the empty tree 
    result+=b[Int64[]]*a[t.index]
    for (subtree,forest) in vector_subtrees_forests[2:end]
        #subtree is a vector containing one value
        tmp=b[subtree[1]]
        for tree in forest 
            tmp*=a[tree]
        end
        result+=tmp
    end
    return result
end


function compose(a::TruncatedBSeries,b::TruncatedBSeries,tree_list::Array{RootedTree_given_by_subtrees}; normalize_stepsize=false)
    series_keys = keys(a)
    series = empty(a)
    for t in series_keys        
            coefficient = compose(a,b,tree_list[t], tree_list)
        if normalize_stepsize
            coefficient /= 2^order(t)
        end
        series[t] = coefficient
    end

    return series
end

function compose(a::TruncatedBSeries,b::TruncatedBSeries,index::Int,data::Data_given_by_ButcherProduct)
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

function compose(a::TruncatedBSeries,b::TruncatedBSeries,data::Data_given_by_ButcherProduct; normalize_stepsize=false)
   
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

end