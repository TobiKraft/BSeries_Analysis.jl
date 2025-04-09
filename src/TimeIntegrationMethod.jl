using LinearAlgebra: dot
using Revise
using RootedTrees_SubtreeStructures
using OrderedCollections
using Counters
export AbstractTimeIntegrationMethod
export RungeKuttaMethod
export elementary_weight,derivative_weight

abstract type AbstractTimeIntegrationMethod end

"""
    RungeKuttaMethod(A, b, c=vec(sum(A, dims=2)))

Represent a Runge-Kutta method with Butcher coefficients `A`, `b`, and `c`.
If `c` is not provided, the usual "row sum" requirement of consistency with
autonomous problems is applied.
"""
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

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                           RootedTrees_given_by_subtrees
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function elementary_weight(tree::RootedTree_given_by_subtrees,tree_list::Vector{RootedTree_given_by_subtrees},rk::RungeKuttaMethod)
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


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                           ButcherProduct
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

function elementary_weight(index::Int,data::Data_given_by_ButcherProduct,rk::RungeKuttaMethod)
    dot(rk.b, derivative_weight(index,data,rk))
end

function derivative_weight(index::Int,data::Data_given_by_ButcherProduct,rk::RungeKuttaMethod) 
    A = rk.A
    c = rk.c

    # Initialize `result` with the identity element of pointwise multiplication `.*`
    result = zero(c) .+ one(eltype(c))

    # Iterate over all subtrees and update the `result` using recursion
    counter=children(index,data)
    for (subtree,quantity) in zip(keys(counter),values(counter))
        tmp = A * derivative_weight(subtree,data,rk)
        result = result .* tmp.^(quantity)
    end

    return result
end


function Base.show(io::IO, rk::RungeKuttaMethod)
    print(io, "RungeKuttaMethod{", eltype(rk), "}")
    if get(io, :compact, false)
        print(io, "(")
        show(io, rk.A)
        print(io, ", ")
        show(io, rk.b)
        print(io, ", ")
        show(io, rk.c)
        print(io, ")")
    else
        print(io, " with\nA: ")
        show(io, MIME"text/plain"(), rk.A)
        print(io, "\nb: ")
        show(io, MIME"text/plain"(), rk.b)
        print(io, "\nc: ")
        show(io, MIME"text/plain"(), rk.c)
        print(io, "\n")
    end
end
