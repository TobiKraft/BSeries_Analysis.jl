using BSeries
using RootedTrees
using RootedTrees_SubtreeStructures
using Base.Threads

a=RootedTrees.rootedtree([1,2,3])
b=RootedTrees.rootedtree([1,2])
const BUFFER_LENGTH = 128
const CANONICAL_REPRESENTATION_BUFFER = Vector{Vector{Int}}()
const PARTITION_ITERATOR_BUFFER_FOREST_T = Vector{Vector{Int}}()
const PARTITION_ITERATOR_BUFFER_FOREST_T_COLORS = Vector{Vector{Bool}}()
const PARTITION_ITERATOR_BUFFER_FOREST_LEVEL_SEQUENCE = Vector{Vector{Int}}()
const PARTITION_ITERATOR_BUFFER_FOREST_COLOR_SEQUENCE = Vector{Vector{Bool}}()
const PARTITION_ITERATOR_BUFFER_SKELETON = Vector{Vector{Int}}()
const PARTITION_ITERATOR_BUFFER_SKELETON_COLORS = Vector{Vector{Bool}}()
const PARTITION_ITERATOR_BUFFER_EDGE_SET = Vector{Vector{Bool}}()
const PARTITION_ITERATOR_BUFFER_EDGE_SET_TMP = Vector{Vector{Bool}}()
function PartitionIterator_neu(t::RootedTree{Int, Vector{Int}},id ) # = Threads.threadid()
    order_t = order(t)
    println("IDH:",Threads.threadid())
    if order_t <= BUFFER_LENGTH
        #id= Threads.threadid()
        buffer_forest_t = PARTITION_ITERATOR_BUFFER_FOREST_T[id]
        resize!(buffer_forest_t, order_t)
        level_sequence = PARTITION_ITERATOR_BUFFER_FOREST_LEVEL_SEQUENCE[id]
        resize!(level_sequence, order_t)
        buffer_skeleton = PARTITION_ITERATOR_BUFFER_SKELETON[id]
        resize!(buffer_skeleton, order_t)
        edge_set = PARTITION_ITERATOR_BUFFER_EDGE_SET[id]
        resize!(edge_set, order_t - 1)
        edge_set_tmp = PARTITION_ITERATOR_BUFFER_EDGE_SET_TMP[id]
        resize!(edge_set_tmp, order_t - 1)
    else
        buffer_forest_t = Vector{Int}(undef, order_t)
        level_sequence = similar(buffer_forest_t)
        buffer_skeleton = similar(buffer_forest_t)
        edge_set = Vector{Bool}(undef, order_t - 1)
        edge_set_tmp = similar(edge_set)
    end

    skeleton = RootedTree(buffer_skeleton, true)
    t_forest = RootedTree(buffer_forest_t, true)
    t_temp_forest = RootedTree(level_sequence, true)
    forest = RootedTrees.PartitionForestIterator(t_forest, t_temp_forest, edge_set_tmp)
    RootedTrees.PartitionIterator{typeof(t), RootedTree{Int, Vector{Int}}}(t, forest, skeleton,
                                                               edge_set, edge_set_tmp)
end
function Base.copy(iterator::PartitionForestIterator)
    return PartitionForestIterator(copy(iterator.t_iter),copy(iterator.t_temp),copy(iterator.edge_set))
end
Threads.resize_nthreads!(CANONICAL_REPRESENTATION_BUFFER,
                             Vector{Int}(undef, BUFFER_LENGTH))

    # PartitionIterator
    Threads.resize_nthreads!(PARTITION_ITERATOR_BUFFER_FOREST_T,
                             Vector{Int}(undef, BUFFER_LENGTH))
    Threads.resize_nthreads!(PARTITION_ITERATOR_BUFFER_FOREST_T_COLORS,
                             Vector{Bool}(undef, BUFFER_LENGTH))
    Threads.resize_nthreads!(PARTITION_ITERATOR_BUFFER_FOREST_LEVEL_SEQUENCE,
                             Vector{Int}(undef, BUFFER_LENGTH))
    Threads.resize_nthreads!(PARTITION_ITERATOR_BUFFER_FOREST_COLOR_SEQUENCE,
                             Vector{Bool}(undef, BUFFER_LENGTH))
    Threads.resize_nthreads!(PARTITION_ITERATOR_BUFFER_SKELETON,
                             Vector{Int}(undef, BUFFER_LENGTH))
    Threads.resize_nthreads!(PARTITION_ITERATOR_BUFFER_SKELETON_COLORS,
                             Vector{Bool}(undef, BUFFER_LENGTH))
    Threads.resize_nthreads!(PARTITION_ITERATOR_BUFFER_EDGE_SET,
                             Vector{Bool}(undef, BUFFER_LENGTH))
    Threads.resize_nthreads!(PARTITION_ITERATOR_BUFFER_EDGE_SET_TMP,
                             Vector{Bool}(undef, BUFFER_LENGTH))

RootedTrees.PartitionIterator(t::RootedTree{Int, Vector{Int}},id )=PartitionIterator_neu(t,id) 
#  =Threads.threadid()
iter_a=RootedTrees.PartitionIterator(a,1)
iter_b=RootedTrees.PartitionIterator(b,3)
for (forest,skeleton) in iter_a,(forest_y,skeleton_y) in iter_b    
    println("Skeleton_X:",skeleton,"   Sekeleton_Y:",skeleton_y)
    print("Forest_X:")
    for x in copy(forest)
        print(x,"  ")
    end
    print("\n Forest_Y:")
    for y in forest_y   
        print( y, "  ")
    end
    println("\n________________________________________________________")
    
end


tree_list=generateTrees_subtrees(3)[1]
subtree_dict=OrderedDict{Vector{Int},Int}()
for x in tree_list
    subtree_dict[x.subtrees]=x.index
end
#=
for (forest,skeleton) in RootedTrees.PartitionIterator(a)
    println("Skeleton_X:",skeleton)
    print("Forest_X:")
    for x in copy(forest)
        print( x, "  ")
    end
    println()
end
substitution_partitions(tree_list,subtree_dict)=#