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
		"Library" => [
			"Functions" => "functions.md",
		]
	]
)
