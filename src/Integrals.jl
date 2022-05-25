module Integrals

	using LinearAlgebraicRepresentation

	const Points = Matrix
	const Cells = Array{Array{Int,1},1}
	const LAR = Union{ Tuple{Points, Cells},Tuple{Points, Cells, Cells} }

	include("./integrals.jl")

end
