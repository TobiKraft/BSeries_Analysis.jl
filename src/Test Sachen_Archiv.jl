using BSeries
"""
A=[0 0 0 0;1/2 0 0 0;0 1/2 0 0;0 0 1 0]
b=[1/6,2/6,2/6,1/6]
c=[0,1/2,1/2,1]
rk=BSeries.RungeKuttaMethod(A,b,c)
rk2=BSeries_Analysis.RungeKuttaMethod(A,b,c)
print("------------------------------------------------------------")
function butcher_test(rk2,tree_list2)
    for index in range(1,37)
        butcher=BSeries_Analysis.elementary_weight(index,tree_list2,rk2)
    end
end
function subtrees_test(rk2,tree_list)
    for index in range(1,37)
        subtree=BSeries_Analysis.elementary_weight(tree_list[index],tree_list,rk2)
    end
end
function level_seq_test(liste,rk)
    for x in liste
        level_seq=BSeries.elementary_weight(RootedTree(x),rk)
    end
end
liste=[]
for index in range(1,37)
    push!(liste,Subtrees_to_Levelsequence(tree_list[index],tree_list))
end
io = IOContext(stdout, :histmin=>30000, :histmax=>400000, :logbins=>true)

print("\n Level-seq: \n")
b=@benchmark level_seq_test(liste,rk)
show(io, MIME("text/plain"), b)
print("\n------------------------------------------------------------")
print("\n Butcher:\n")
b=@benchmark butcher_test(rk2,tree_list2)
show(io, MIME("text/plain"), b)
print("\n------------------------------------------------------------")
print("\n Subtree:\n")
b=@benchmark subtrees_test(rk2,tree_list)
show(io, MIME("text/plain"), b)
print("\n------------------------------------------------------------")
print("\n------------------------------------------------------------")
print("\n------------------------------------------------------------\n")
print(b)

Test2: 
series_butcher=BSeries_Analysis.TruncatedBSeries{Union{Int64,Vector{Int64}},Rational}()
series_subtrees=BSeries_Analysis.TruncatedBSeries{Union{Int64,Vector{Int64}},Rational}()
series_butcher[Int64[]]=1
series_subtrees[Int64[]]=1
len=length(data.tree_list)

for t in range(1,len)
    series_butcher[t]=exact_value(t,data)
end
for tree in tree_list
    series_subtrees[tree.index]=exact_value(tree,tree_list)
end

for tree in tree_list
    subtree_value=series_subtrees[tree.index]
    butcher_index=Subtrees_to_ButcherProduct(tree,tree_list,data)
    butcher_value=series_butcher[butcher_index]
    bool=(subtree_value==butcher_value)
    if !bool
     print("Subtree:", subtree_value, "     Butcher:",butcher_value,"   ", bool,"\n")
    end
end
series_butcher2=BSeries_Analysis.compose(series_butcher,series_butcher,data)
series_subtrees2=BSeries_Analysis.compose(series_subtrees,series_subtrees,tree_list)
print("Index:Int64[]    Subtree:", series_subtrees2[Int64[]], "     Butcher:",series_butcher2[Int64[]],"\n")
for tree in tree_list
    subtree_value=series_subtrees2[tree.index]
    butcher_index=Subtrees_to_ButcherProduct(tree,tree_list,data)
    butcher_value=series_butcher2[butcher_index]
    bool=(subtree_value==butcher_value)
    print("Index:",tree.index,"   Subtree:", subtree_value, "     Butcher:",butcher_value,"   ", bool,"\n")
end"""