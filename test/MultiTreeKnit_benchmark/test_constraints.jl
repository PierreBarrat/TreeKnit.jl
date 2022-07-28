using TreeKnit
using TreeTools
using Test

println("##### testing constraints #####")

function join_sets_slow(input_sets::Vector{Vector{Vector{String}}})
    start_set = input_sets[1]
    for i in 2:length(input_sets)
        joint_sets = Vector{String}[]
        to_be_joint_set = [Set{String}(s2) for s2 in input_sets[i]]
        for s1 in start_set
            nodes = length(s1)
            while nodes>0
                for s2 in to_be_joint_set
                    joint_set = intersect(Set{String}(s1), s2)
                    if !isempty(joint_set)
                        s2 = setdiff(s2,joint_set)
                        nodes -= length(joint_set)
                        append!(joint_sets, [sort(collect(joint_set))])
                    end
                end
            end
        end
        start_set = joint_sets
    end
    return sort(start_set; lt=clt)
end

@testset "MCC join is accurate" begin
	trees, arg = get_trees(3, 10000);
    rMCC_list = TreeKnit.get_real_MCCs(3, arg)
    slow = join_sets_fast([rMCC_list[1], rMCC_list[2]])
    fast = TreeKnit.join_sets([rMCC_list[1], rMCC_list[2]])
    @test slow == fast
end