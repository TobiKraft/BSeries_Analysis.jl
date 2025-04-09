
#=
#= Reference:
    This file is a Julia implementation of the file "BSeries.py" [2] of the already existing Python package "orderConditions" of Valentin Dallerit [1].
    It is nearly a one to one implentation of the Python code in Julia. Some changes were made due to the differences in Julia and Python.
    The plotting part of the code is left out.
=#
=#




using Revise
using RootedTrees_SubtreeStructures
using OrderedCollections: OrderedDict
using Counters

export Data_BSeries
export set_order
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
Data_BSeries(4,8,Dict{Int64,Number}(),Data_given_by_CircProduct([(1,0),(1,1),(2,1),(1,2),(3,1),(4,1),(1,3),(1,4)],4,[1,2,3,5,9]))
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

#=
References

[1] Valentin Dallerit: Gitlab package orderConditions (Programming Language Python)
    https://gitlab.com/v_dallerit/orderconditions/, accessed 09.04.2025

[2] Valentin Dallerit: Gitlab file "Bseries.py"
    https://gitlab.com/v_dallerit/orderconditions/-/blob/master/orderConditions/BSeries.py, accessed 09.04.2025
=#