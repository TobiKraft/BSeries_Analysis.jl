
#=
#= Reference:
    This file is a Julia implementation of the file "BSeries.py" [2] of the already existing Python package "orderConditions" of Valentin Dallerit [1].
    It is nearly a one to one implentation of the Python code in Julia. Some changes were made due to the differences in Julia and Python.
    The displaying part of the code is left out.
=#
=#




using Revise
using RootedTrees_SubtreeStructures
using OrderedCollections: OrderedDict
using Counters

export Data_BSeries
export set_order
export exact, y, hf
export compo_hf
export _lambda_sub,substitution
export modified_equation,modifying_integrator
mutable struct Data_BSeries
    order_max::Int
    len::Int
    lambda_cache::Dict{Int,Number}
    circProduct_data::Data_given_by_CircProduct
    function Data_BSeries()
        new(0,1,Dict{Int,Number}(),Data_given_by_CircProduct())
    end
end

"""
    set_order(order::Int)   ->  Data_BSeries
Sets the order to be used for the B-Series.\\
This function must be run first to set the size of the B-Series and generate associated trees, i. e. Data_given_by_CircProduct.
# Example
```julia
julia>data=set_order(4)
Data_BSeries(4,9,Dict{Int64,Number}(),Data_given_by_CircProduct([(0,0),(2,0),(2,2),(3,2),(2,3),(4,2),(5,2),(2,4),(2,5)],4,[2,3,4,6,10]))
```
"""
function set_order(order::Int)
    data_new=Data_BSeries()
    data_new.circProduct_data=generateTrees_CircProduct(order)
    data_new.order_max=order
    data_new.len=length(data_new.circProduct_data.tree_list)
    return data_new
end

#--------------------------------------------------------------------------------
#                    Special B-Series
#--------------------------------------------------------------------------------
"""
    exact(data::Data_BSeries)   ->  Array{Rational{Int}}
Returns the B-Series of the exact solution at time ``t_{n+1}:y_{n+1}``
# Example
```julia
julia>data=set_order(4)
julia>exact(data)
8-element Vector{Rational{Int64}}:
  1
 1//2
 1//3
 1//6
 1//4
 1//8
 1//12
 1//24
```
"""
function exact(data::Data_BSeries)
    return [1//RootedTrees_SubtreeStructures.density(i,data.circProduct_data) for i in range(1,data.len)]
end

"""
    y(data::Data_BSeries)   ->  Vector{Int} representing a B-Series
Returns the B-Series of the exact solution
"""
function y(data::Data_BSeries,::Type{T}=Int) where{T<:Number}
    yn=zeros(T,data.len)
    yn[1]=T(1)
    return yn
end

"""
    hf(data::Data_BSeries)  ->  Vector{Int} representing a B-Series
Returns the B-Series of ``hf(y)``.
"""
function hf(data::Data_BSeries,::Type{T}=Int) where{T<:Number}
    f=zeros(T,data.len)
    f[2]=T(1)
    return f
end
#--------------------------------------------------------------------------------
#                    Composition Rule
#--------------------------------------------------------------------------------
function _lambda_compo(alpha::Dict{Int,Number},tree::Int,data::Data_BSeries)
    if tree==1
        return Dict{Int,Number}()
    end
    if tree==2
        return Dict(1=>Integer(1))
    end
    (u,v)=decompose(tree,data.circProduct_data)
    lambda_u=data.lambda_cache[u]
    lambda_v=data.lambda_cache[v]
    ans=Dict{Int,Number}
end

function compo_hf(alpha,data::Data_BSeries)
    ans=zeros(Rational{Int},data.len)
    ans[2]=1
    for tree in range(3,data.len)
        (u,v)=decompose(tree,data.circProduct_data)
        ans[tree]=alpha[v]*ans[u]
    end
    return ans
end
#--------------------------------------------------------------------------------
#                    Composition Rule
#--------------------------------------------------------------------------------
function _lambda_sub(data::Data_BSeries)
    l=[Dict{Int,Vector{Vector{Int}}}() for x in range(1,data.len)]
    l[1][1]=[[]]
    l[2][2]=[[2]]
    for tree in range(3,data.len)
        (u,v)=decompose(tree,data.circProduct_data)
        for (skeleton_u,forestlist_u) in pairs(l[u]),(skeleton_v,forestlist_v) in pairs(l[v])
            circ_uv::Int=circ(skeleton_u,skeleton_v,data.circProduct_data)
            get!(l[tree],circ_uv,Vector{Int}())
            #l[tree][circ_uv]=Vector{Int}()
            merge_uv=merge_root(skeleton_u,skeleton_v,data.circProduct_data)
            #l[tree][merge_uv]=Vector{Int}()
            get!(l[tree],merge_uv,Vector{Int}())
            for forest_u in forestlist_u, forest_v in forestlist_v
                tmp_circ=copy(forest_u)
                append!(tmp_circ,forest_v)
                push!(l[tree][circ_uv],tmp_circ)
                tmp_merge=[circ(forest_u[1],forest_v[1],data.circProduct_data)]
                append!(tmp_merge,forest_u[2:end],forest_v[2:end])
                push!(l[tree][merge_uv],tmp_merge)                
            end
        end
    end
    return l

end
"""
    substitution(a,b,data::Data_BSeries)    ->  Vector{Int} representing a B-Series

"""
function substitution(a::Union{Vector{Number},Dict{Int,Number}},b::Union{Vector{Number},Dict{Int,Number}},data::Data_BSeries)
    l=_lambda_sub(data)
    ans=zeros(Rational{Int},data.len)
    for i in range(1,data.len)
        for (b_tree,forest_list) in pairs(l[i])
            for forest in forest_list
                k=1
                for a_tree in forest
                    k*=a[a_tree]
                end
                ans[i]+=k*b[b_tree]
            end
        end
    end
    return ans
end

#--------------------------------------------------------------------------------
#                    Modifying Equation/Integrator
#--------------------------------------------------------------------------------
function modified_equation(a::Vector{V},data::Data_BSeries) where {V<:Number}
    ex_series=exact(data)
    series=zeros(V,data.len)
    series[2]=a[2]
    l=_lambda_sub(data)
    for t in range(3,data.len)
        tmp=0
        for (b_tree,list_forest) in pairs(l[t])
            for forest in list_forest
                k=1
                for a_tree in forest
                    k*=series[a_tree]
                end
                tmp+=k*ex_series[b_tree]
            end
        end
        series[t]=a[t]-tmp
    end
    return series
end

function modifying_integrator(a,data::Data_BSeries)
    ex_series=exact(data)
    series=zeros(Rational{Int},data.len)
    series[2]=a[2]
    l=lambda_sub()
    for t in range(3,data.len)
        tmp=0
        for (b_tree,list_forest) in pairs(l[t])
            for forest in list_forest
                k=1
                for a_tree in list_forest
                    k*=series[a_tree]
                end
                tmp+=k*a[b_tree]
            end
        end
        series[t]=ex_series[t]-tmp
    end
end
#
#function compo(alpha,beta,data)
#    ans=zeros(data.len)
#    _data.lambda_cache = [{}] * _data.len
#end
#=
References

[1] Valentin Dallerit: Gitlab package orderConditions (Programming Language Python)
    https://gitlab.com/v_dallerit/orderconditions/, accessed 09.04.2025

[2] Valentin Dallerit: Gitlab file "Bseries.py"
    https://gitlab.com/v_dallerit/orderconditions/-/blob/master/orderConditions/BSeries.py, accessed 09.04.2025
=#