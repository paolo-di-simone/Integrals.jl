push!(LOAD_PATH,"../src/")

using Integrals, Documenter
makedocs(
	format = Documenter.HTML(
		prettyurls = get(ENV, "CI", nothing) == "true"
	),
	sitename = "Integrals.jl",
	modules = [Integrals],
	pages=[
    	"Home" => "index.md"
    ]
)

deploydocs(;
    repo="github.com/paolo-di-simone/Integrals.jl",
)
