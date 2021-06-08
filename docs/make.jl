using RecombTools
using Documenter

Documenter.makedocs(;
	source = "src",
	clean = true,
	doctest = true,
	modules = Module[RecombTools],
	repo = "",
	highlightsig = true,
	sitename = "RecombTools documentation",
	expandfirst = [],
	pages = [
		"Index" => "index.md",
		"Usage" => "usage.md",
		"Resolving" => "resolving.md",
		"Library" => [
			"Types" => "types.md",
			"Functions" => "functions.md",
		]
	]
)
