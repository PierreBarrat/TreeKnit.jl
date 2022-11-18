using TreeKnit
using Documenter

Documenter.makedocs(;
	source = "src",
	clean = true,
	doctest = true,
	modules = Module[TreeKnit],
	repo = "",
	highlightsig = true,
	sitename = "TreeKnit documentation",
	expandfirst = [],
	pages = [
		"Index" => "index.md",
		"Usage" => [
			"Overview" => "overview.md",
			"MCCs" => "mccs.md",
			"Important options" => "options.md",
			"MultiTreeKnit" => "multitreeknit.md",
			"Visualizing MCCs" => "visualization.md",
		],
		"Under the hood" => [
			"`opttrees`" => "opttrees.md",
			"`runopt`" => "runopt.md",
			"Resolving" => "resolving.md",
		],
		"Library" => [
			"OptArgs" => "types.md",
			"Functions" => "functions.md",
		],
	]
)

deploydocs(
    repo = "github.com/PierreBarrat/TreeKnit.jl.git";
    versions = nothing,
)
