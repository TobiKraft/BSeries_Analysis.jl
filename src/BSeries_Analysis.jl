module BSeries_Analysis
using OrderedCollections
using Revise
using Counters
using LinearAlgebra: dot
using RootedTrees_SubtreeStructures


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


struct TruncatedBSeries{T<:Union{Int,Vector{Any}},V} <:AbstractDict{T,V}
    coef::OrderedDict{T,V}
    TruncatedBSeries{T, V}() where {T, V} = TruncatedBSeries{T, V}(OrderedDict{T, V}())
end

function compose(a::AbstractDict,b::AbstractDict,t::RootedTree_given_by_subtrees,tree_list::Array{RootedTree_given_by_subtrees})
    result = zero(first(values(a)) * first(values(b)))
    for (subtree,forest) in orderedSubtrees_and_Forests(t,tree_list)
        tmp=b[subtree]
        for tree in forest 
            tmp*=a[tree]
        end
        result+=tmp
    end
    return result
end

function exact_value(tree::RootedTree_given_by_subtrees,tree_list::Array{RootedTree_given_by_subtrees})
    return 1/density(tree,tree_list)
end

function compose(a::TruncatedBSeries,b::TruncatedBSeries,tree_list::Array{RootedTree_given_by_subtrees}; normalize_stepsize=false)
    series_keys = keys(b)
    series = empty(b)

    for t in series_keys
        coefficient = compose(b, a, t, tree_list)
        if normalize_stepsize
            coefficient /= 2^order(t)
        end
        series[t] = coefficient
    end

    return series
end

function compose(a::AbstractDict,b::AbstractDict,t::Int,data::Data_given_by_ButcherProduct)
    result = zero(first(values(a)) * first(values(b)))
    for (subtree,forest) in orderedSubtrees_and_Forests(t,data)
        tmp=b[subtree]
        for tree in forest 
            tmp*=a[tree]
        end
        result+=tmp
    end
    return result
end
function compose(a::TruncatedBSeries,b::TruncatedBSeries,data::Data_given_by_ButcherProduct ; normalize_stepsize=false)
    series_keys = keys(b)
    series = empty(b)   

    for t in series_keys
        coefficient = compose(b, a, t, data)
        if normalize_stepsize
            coefficient /= 2^order(t)
        end
        series[t] = coefficient
    end

    return series
end
end