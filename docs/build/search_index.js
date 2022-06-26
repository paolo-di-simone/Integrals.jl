var documenterSearchIndex = {"docs":
[{"location":"index.html#Integrals.jl","page":"Home","title":"Integrals.jl","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Questo modulo implementa un metodo di integrazione finita di polinomi. ","category":"page"},{"location":"index.html#Interfacce-principali","page":"Home","title":"Interfacce principali","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"II(P::LAR, alpha::Int, beta::Int, gamma::Int, signedInt=false)::Float64\nIII(P::LAR, alpha::Int, beta::Int, gamma::Int, signedInt::Bool=false)::Float64\nsurface(P::LAR, signedInt::Bool=false)::Float64\nvolume(P::LAR, signedInt::Bool=false)::Float64\nfirstMoment(P::LAR)::Array{Float64,1}\nsecondMoment(P::Lar.LAR)::Array{Float64,1}\ninertiaProduct(P::LAR)::Array{Float64,1}\ncentroid(P::LAR)::Array{Float64,1}\ninertiaMoment(P::LAR)::Array{Float64,1}","category":"page"},{"location":"index.html#Documentazione-funzioni","page":"Home","title":"Documentazione funzioni","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Integrals.M","category":"page"},{"location":"index.html#Integrals.M","page":"Home","title":"Integrals.M","text":"M(alpha::Int, beta::Int)::Float64\n\nModule for integration of polynomials over 3D volumes and surfaces\n\n\n\n\n\n","category":"function"},{"location":"index.html","page":"Home","title":"Home","text":"Integrals.TT","category":"page"},{"location":"index.html#Integrals.TT","page":"Home","title":"Integrals.TT","text":"TT(tau::Array{Float64,2}, alpha::Int, beta::Int, gamma::Int, signedInt::Bool=false)\n\nThe main integration routine \n\n\n\n\n\n","category":"function"},{"location":"index.html","page":"Home","title":"Home","text":"Integrals.II","category":"page"},{"location":"index.html#Integrals.II","page":"Home","title":"Integrals.II","text":"II(P::Lar.LAR, alpha::Int, beta::Int, gamma::Int, signedInt=false)\n\nBasic integration function on 2D plane.\n\nExample  unit 3D triangle\n\njulia> V = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]\n3×3 Array{Float64,2}:\n 0.0  1.0  0.0\n 0.0  0.0  1.0\n 0.0  0.0  0.0\n \njulia> FV = [[1,2,3]]\n1-element Array{Array{Int64,1},1}:\n [1, 2, 3]\n \njulia> P = V,FV\n([0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0], Array{Int64,1}[[1, 2, 3]])\n\njulia> Integrals.II(P, 0,0,0)\n0.5\n\n\n\n\n\n","category":"function"},{"location":"index.html","page":"Home","title":"Home","text":"Integrals.III","category":"page"},{"location":"index.html#Integrals.III","page":"Home","title":"Integrals.III","text":"III(P::LAR, alpha::Int, beta::Int, gamma::Int, signedInt::Bool=false)::Float64\n\nBasic integration function on 3D space.\n\nExample # unit 3D tetrahedron\n\njulia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]\n3×4 Array{Float64,2}:\n 0.0  1.0  0.0  0.0\n 0.0  0.0  1.0  0.0\n 0.0  0.0  0.0  1.0\n\njulia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]]\n4-element Array{Array{Int64,1},1}:\n [1, 2, 4]\n [1, 3, 2]\n [4, 3, 1]\n [2, 3, 4]\n\njulia> P = V,FV\n([0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0], \nArray{Int64,1}[[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]])\n\njulia> Integrals.III(P, 0,0,0)\n0.16666666666666674\n\n\n\n\n\n","category":"function"},{"location":"index.html","page":"Home","title":"Home","text":"Integrals.surface","category":"page"},{"location":"index.html#Integrals.surface","page":"Home","title":"Integrals.surface","text":"surface(P::Lar.LAR, signedInt::Bool=false)::Float64\n\nsurface integral on polyhedron P.\n\nExample unit 3D triangle\n\njulia> V = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]\n3×3 Array{Float64,2}:\n 0.0  1.0  0.0\n 0.0  0.0  1.0\n 0.0  0.0  0.0\n \njulia> FV = [[1,2,3]]\n1-element Array{Array{Int64,1},1}:\n [1, 2, 3]\n \njulia> P = V,FV\n([0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0], Array{Int64,1}[[1, 2, 3]])\n\njulia> Integrals.surface(P)\n0.5\n\n\n\n\n\n","category":"function"},{"location":"index.html","page":"Home","title":"Home","text":"Integrals.volume","category":"page"},{"location":"index.html#Integrals.volume","page":"Home","title":"Integrals.volume","text":"volume(P::LAR, signedInt::Bool=false)::Float64\n\nvolume integral on polyhedron P.\n\nExample # unit 3D tetrahedron\n\njulia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]\n3×4 Array{Float64,2}:\n 0.0  1.0  0.0  0.0\n 0.0  0.0  1.0  0.0\n 0.0  0.0  0.0  1.0\n\njulia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]]\n4-element Array{Array{Int64,1},1}:\n [1, 2, 4]\n [1, 3, 2]\n [4, 3, 1]\n [2, 3, 4]\n\njulia> P = V,FV\n([0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0], \nArray{Int64,1}[[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]])\n\njulia> Integrals.volume(P)\n0.16666666666666674\n\n\n\n\n\n","category":"function"},{"location":"index.html","page":"Home","title":"Home","text":"Integrals.firstMoment","category":"page"},{"location":"index.html#Integrals.firstMoment","page":"Home","title":"Integrals.firstMoment","text":"firstMoment(P::Lar.LAR)::Array{Float64,1}\n\nFirst moments as terms of the Euler tensor. Remember that the integration algorithm is a boundary integration. Hence the model must be a boundary model. In this case, a 2-complex of triangles. \n\nExample # unit 3D tetrahedron\n\njulia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];\n\njulia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];\n\njulia> P = V,FV;\n\njulia> Integrals.firstMoment(P)\n3-element Array{Float64,1}:\n 0.0416667\n 0.0416667\n 0.0416667\n\n\n\n\n\n","category":"function"},{"location":"index.html","page":"Home","title":"Home","text":"Integrals.secondMoment","category":"page"},{"location":"index.html#Integrals.secondMoment","page":"Home","title":"Integrals.secondMoment","text":"secondMoment(P::Lar.LAR)::Array{Float64,1}\n\nSecond moments as terms of the Euler tensor.\n\nExample # unit 3D tetrahedron\n\njulia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];\n\njulia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];\n\njulia> P = V,FV;\n\njulia> Integrals.secondMoment(P)\n3-element Array{Float64,1}:\n 0.0166667\n 0.0166667\n 0.0166667\n\n\n\n\n\n","category":"function"},{"location":"index.html","page":"Home","title":"Home","text":"Integrals.inertiaProduct","category":"page"},{"location":"index.html#Integrals.inertiaProduct","page":"Home","title":"Integrals.inertiaProduct","text":"inertiaProduct(P::Lar.LAR)::Array{Float64,1}\n\nInertia products as terms of the Euler tensor.\n\nExample # unit 3D tetrahedron\n\njulia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];\n\njulia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];\n\njulia> P = V,FV;\n\njulia> Integrals.inertiaProduct(P)\n3-element Array{Float64,1}:\n 0.00833333\n 0.00833333\n 0.00833333\n\n\n\n\n\n","category":"function"},{"location":"index.html","page":"Home","title":"Home","text":"Integrals.centroid","category":"page"},{"location":"index.html#Integrals.centroid","page":"Home","title":"Integrals.centroid","text":"centroid(P::Lar.LAR)::Array{Float64,1}\n\nBarycenter or centroid of polyhedron P.\n\nExample # unit 3D tetrahedron\n\njulia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];\n\njulia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];\n\njulia> P = V,FV;\n\njulia> Integrals.centroid(P)\n3-element Array{Float64,1}:\n 0.25\n 0.25\n 0.25\n\n\n\n\n\n","category":"function"},{"location":"index.html","page":"Home","title":"Home","text":"Integrals.inertiaMoment","category":"page"},{"location":"index.html#Integrals.inertiaMoment","page":"Home","title":"Integrals.inertiaMoment","text":"inertiaMoment(P::Lar.LAR)::Array{Float64,1}\n\nInertia moments  of polyhedron P.\n\nExample # unit 3D tetrahedron\n\njulia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];\n\njulia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];\n\njulia> P = V,FV;\n\njulia> Integrals.inertiaMoment(P)\n3-element Array{Float64,1}:\n 0.0333333\n 0.0333333\n 0.0333333\n\n\n\n\n\n","category":"function"}]
}