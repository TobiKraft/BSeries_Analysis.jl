
using RootedTrees_SubtreeStructures
using BenchmarkTools
using BSeries_Analysis  
A=[0 0;1//2 0]
b=[0,1]
c=[0,1//2]
i=15
rk=BSeries_Analysis.RungeKuttaMethod(A,b,c)
(tree_list,order_list)=generateTrees_subtrees(i)
io=IOContext(stdout,:logbins=>true)
for index_i in range(1,i)
    a=index_i
    println("___________________________________________________________________")
    println("Max_Order=",index_i)
    #print(a)
    bm1=@benchmark BSeries_Analysis.bseries_ohne_order(rk,$a,tree_list) 
    #println("\nohne_order")
    #show(io,MIME("text/plain"),bm1)
    #println("\nmit_order")
    bm2=@benchmark BSeries_Analysis.bseries_mit_order(rk,$a,tree_list,order_list)
    #show(io,MIME("text/plain"),bm2)
    show(judge(median(bm1),median(bm2)))
    println("\n___________________________________________________________________")
end