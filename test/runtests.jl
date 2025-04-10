using BSeries_Analysis
using RootedTrees_SubtreeStructures
using Test
using OrderedCollections: OrderedDict
@testset "BSeries_Analysis.jl" begin
    @testset "TimeIntegrationMethod.jl" begin
        A=[0 0;1//2 0]
        b=[0,1]
        c=[0,1//2]
        rk=RungeKuttaMethod(A,b,c)
        expected=[1//1,1//2,1//4,0//1,1//8,0//1,0//1,0//1,1//16,0//1,0//1,0//1,0//1,0//1,0//1,0//1,0//1]
        expected_circProduct=[1//1,1//1,1//2,1//4,0//1,1//8,0//1,0//1,0//1,1//16,0//1,0//1,0//1,0//1,0//1,0//1,0//1,0//1]
        (tree_list,order_list)=generateTrees_subtrees(5)
        data=generateTrees_CircProduct(5)
        @testset "elementary_weight" begin
            for (i,x) in enumerate(tree_list)
                ans_sub=elementary_weight(x,tree_list,rk)
                ans_circProduct=elementary_weight(i+1,data,rk)
                @test expected[i]==ans_sub
                @test expected[i]==ans_circProduct
            end
        end
    end
    @testset "BSeries_RootedTrees_given_by_subtrees.jl" begin
        
        @testset "RungeKuttaMethod_elementary_weights" begin
            A=[0 0;1//2 0]
            b=[0,1]
            c=[0,1//2]
            rk=RungeKuttaMethod(A,b,c)

        end
        (tree_list,order_list)=generateTrees_subtrees(6)
        @testset "exact_value" begin
            expected=[1,1//2,1//3,1//6,1//4,1//8,1//12,1//24]
            for i in range(1,8)
                tmp_ans=exact_value(tree_list[i],tree_list)
                tmp_ans2=exact_value(i,tree_list)
                @test expected[i]==tmp_ans2
                @test expected[i]==tmp_ans
            end
        end
        @testset "bseries" begin
            A=[0 0;1//2 0]
            b=[0,1]
            c=[0,1//2]
            rk=RungeKuttaMethod(A,b,c)
            ans=bseries(rk,3,tree_list)
            expected=OrderedDict{Int64, Rational{Int64}}()
            expected[0]=1
            expected[1]=1
            expected[2]=1//2
            expected[3]=1//4
            expected[4]=0
            @test expected==ans        
        end
    end
    @testset "BSeries_RootedTrees_given_by_ButcherProduct.jl" begin
        data=generateTrees_CircProduct(5)
        @testset "exact_series_butcher" begin
            expected=Dict{Int,Number}()
            expected_series=OrderedDict{Int,Rational{Integer}}()
            expected_series_float=OrderedDict{Int,Float64}()
            expected_values=[1,1//2,1//3,1//6,1//4,1//8,1//12,1//24]
            for (i,x) in enumerate(expected_values)
                expected[i]=convert(Int,1/x)
                expected_series[i]=x
                expected_series_float[i]=convert(Float64,x)
            end
            #series=exact_series_butcher!(4,data)
            #series_float=exact_series_butcher!(4,data,Float64)
            #@test expected==data.density_Dict
            #@test expected_series==series
            #@test expected_series_float==series_float
            #@test_throws InexactError exact_series_butcher(5,data,Integer)
        end
    end
    @testset "BSeries_CircProduct.jl" begin
        data=set_order(4)
        @test data.order_max==4
        @test data.len==9
        @test data.lambda_cache==Dict{Int,Number}()
        @test data.circProduct_data==Data_given_by_CircProduct([(0,0),(2,0),(2,2),(3,2),(2,3),(4,2),(5,2),(2,4),(2,5)],4,[2,3,4,6,10])
        exact_series=exact(data)
        @test exact_series==[1//1,1//1,1//2,1//3,1//6,1//4,1//8,1//12,1//24]
        @testset "Composition-Rule" begin
            data=set_order(5)
            alpha=exact(data)
            beta=compo_hf(alpha,data)
            @test beta==[0,1,1,1,1//2,1,1//2,1//3,1//6,1,1//2,1//3,1//6,1//4,1//4,1//8,1//12,1//24]
        end
    end
end
