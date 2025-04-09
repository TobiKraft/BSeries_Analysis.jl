using Revise
using RootedTrees_SubtreeStructures
using OrderedCollections: OrderedDict
using Counters
include("TimeIntegrationMethod.jl")

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                           BSeries
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
    exact_value_butcher(index::Int,data::Data_given_by_CircProduct) -> Rational{Integer}
``index`` needs to an Integer representing the index of the wanted tree  in ``data.tree_list``.\\
For a given autonomous problem
* y'=f(y),        y(t_0)=y_0

and a tree tau, this function returns the coefficient alpha(tau) as a Rational of the BSeries, that is representing the exact solution of the given problem.\\
For repeatedly calls use exact_value_butcher! instead, which changes ``data.density_Dict``.
# Example
```julia
julia>data=generateTrees_butcher(5)
julia>exact_value_butcher(6,data)
1//8
```       
"""
function exact_value_butcher(index::Int,data::Data_given_by_CircProduct)
    return 1//density(index,data)
end


"""
    exact_value_butcher!(index::Int,data::Data_given_by_CircProduct) -> Rational{Integer}
See also [`exact_value_butcher`] for more details.\\
Call this instead of ``exact_value_butcher`` to get the coefficient of ``data.tree_list[i]`` by accessing ``data.density_Dict``.\\
Needs some more memory and time than ``exact_value_butcher`` during the first call, but is much faster in repeated use.
# Example
```julia
julia>data=generateTrees_butcher(5)
julia>exact_value_butcher(6,data)
1//8
```    
"""
function exact_value_butcher!(index::Int,data::Data_given_by_CircProduct)
    return 1//density!(index,data)
end

"""
    exact_series_butcher(order::Int,data::Data_given_by_CircProduct,V::DataType=Rational{Integer}) -> OrderedDict{Int,V}
Returns an OrderedDict{Int,``V``} where ``V`` is the DataType Rational{Integer} per default. ``V`` needs to be a DataType a Rational can be converted into.\\
Returns the coefficients alpha(tau) of the BSeries representing the exact solution of the autonomous problem
* y'=f(y),        y(t_0)=y_0

for every tree tau in ``data.tree_list`` up to the given ``order``.\\
This version is faster if ``data.density_Dict`` is already computed for every tree.
"""
function exact_series_butcher(order::Int,data::Data_given_by_CircProduct,V::DataType=Rational{Integer})
    series=OrderedDict{Int,V}()
    for x in range(1,(data.index_first[order+1]-1))
        series[x]=convert(V,exact_value_butcher(x,data))
    end
    return series
end

"""
    exact_series_butcher!(order::Int,data::Data_given_by_CircProduct,V::DataType=Rational{Integer}) -> OrderedDict{Int,V}
Returns an OrderedDict{Int,``V``} where ``V`` is the DataType Rational{Integer} per default. ``V`` needs to be a DataType a Rational can be converted into.\\
See [`exact_series_butcher`] for more details.\\
Call this instead of ``exact_value_butcher`` to get the coefficient of ``data.tree_list[i]`` by accessing ``data.density_Dict``.\\
This version is much faster if ``data.density_Dict`` is already computed for every tree.
"""
function exact_series_butcher!(order::Int,data::Data_given_by_CircProduct,V::DataType=Rational{Integer})
    series=OrderedDict{Int,V}()
    for x in range(1,(data.index_first[order+1]-1))
        series[x]=convert(V,exact_value_butcher!(x,data))
    end
    return series
end