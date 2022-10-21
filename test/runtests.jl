using Test

@testset verbose=true "TreeKnit" begin
	@testset "Resolving" begin
		println("# Resolving")
		include("resolving/test.jl")
	end

	@testset "SplitGraph" begin
		println("# SplitGraph")
		include("splitgraph/test.jl")
	end

	@testset "src/main.jl" begin
		println("# Main")
		include("main/test.jl")
	end
  
	@testset "src/mcc_tools.jl" begin
		println("# MCC tools")
		include("mcc_tools/test.jl")
	end

	@testset "SimpleReassortmentGraph" begin
		println("# SimpleReassortmentGraph")
		include("SRG/test.jl")
	end

	@testset "NY data" begin
		println("# NY data")
		include("NYdata/test.jl")
	end

	@testset "MultiTreeKnit" begin
		println("# MultiTreeKnit")
		include("MultiTreeKnit/test.jl")
	end
end
