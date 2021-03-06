"""
    M(alpha::Int, beta::Int)::Float64
Module for integration of polynomials over 3D volumes and surfaces
"""
function M(alpha::Int, beta::Int)::Float64
    a = 0
    for l=0:(alpha + 1)
        a += binomial(alpha+1,l) * (-1)^l/(l+beta+1)
    end
    return a/(alpha + 1)
end

function s4(a, b, h, k, m, i, j)
    ss4 = 0.0
    for l=0:m
		ss4 += binomial(m,l) * a[3]^(m-l) * b[3]^l * M(h+k+m-i-j-l, i+j+l)
    end
    return ss4
end

function s3(a, b, h, k, m, i)
    ss3 = 0.0
    for j=0:k
		ss3 += binomial(k,j) * a[2]^(k-j) * b[2]^j * s4(a, b, h, k, m, i, j)
    end
    return ss3
end

function s2(a, b, h, k, m)
    ss2 = 0.0
    for i=0:h 
		ss2 += binomial(h,i) * a[1]^(h-i) * b[1]^i * s3(a, b, h, k, m, i);
    end
    return ss2
end

function s1(a, b, alpha, beta, gamma, vo)
	ss1 = 0.0
	for x=0:((alpha+1) * (beta+1) * (gamma+1))
		h = x ÷ ((beta+1) * (gamma+1))
		k = (x - h * (beta+1) * (gamma+1)) ÷ (gamma + 1) 
		m = (x - h * (beta+1) * (gamma+1)) % (gamma + 1) 
		ss1 += binomial(alpha,h) * binomial(beta,k) * binomial(gamma,m) * vo[1]^(alpha-h) * vo[2]^(beta-k) * vo[3]^(gamma-m) * s2(a, b, h, k, m)
	end
	return ss1
end

""" 
	TT(tau::Array{Float64,2}, alpha::Int, beta::Int, gamma::Int, signedInt::Bool=false)
The main integration routine 
"""
function TT(tau::Array{Float64,2}, alpha::Int, beta::Int, gamma::Int, signedInt::Bool=false)
    vo,va,vb = tau[:,1],tau[:,2],tau[:,3]
    a = va - vo
    b = vb - vo
    c = cross(a,b)
    if alpha == 0 && beta == 0 && gamma == 0 #surface
        area_tt = M(0,0)
    elseif alpha == 1 && beta == 0 && gamma == 0 #volume
        area_tt = vo[1]*M(0,0) + a[1]*M(1,0) + b[1]*M(0,1)
    elseif alpha == 2 && beta == 0 && gamma == 0 #firstMoment
        area_tt = vo[1]^2*M(0,0) + 
                  2*vo[1]*(a[1]*M(1,0) + b[1]*M(0,1)) + 
                  a[1]^2*M(2,0) + 2*a[1]*b[1]*M(1,1) + b[1]^2*M(0,2)
    elseif alpha == 1 && beta == 1 && gamma == 0 #firstMoment
        area_tt = vo[1]*vo[2]*M(0,0) + 
                  vo[1]*(a[2]*M(1,0) + b[2]*M(0,1)) + 
                  vo[2]*(a[1]*M(1,0) + b[1]*M(0,1)) + 
                  a[1]*(a[2]*M(2,0) + b[2]*M(1,1)) + b[1]*(a[2]*M(1,1) + b[2]*M(0,2))
    elseif alpha == 1 && beta == 0 && gamma == 1 #firstMoment
        area_tt = vo[1]*vo[3]*M(0,0) + 
                  vo[1]*(a[3]*M(1,0) + b[3]*M(0,1)) + 
                  vo[3]*(a[1]*M(1,0) + b[1]*M(0,1)) + 
                  a[1]*(a[3]*M(2,0) + b[3]*M(1,1)) + b[1]*(a[3]*M(1,1) + b[3]*M(0,2))
    elseif alpha == 3 && beta == 0 && gamma == 0 #secondMoment
        area_tt = vo[1]^3*M(0,0) + 
                  3*vo[1]^2*(a[1]*M(1,0) + b[1]*M(0,1)) + 
                  3*vo[1]*(a[1]^2*M(2,0) + 2*a[1]*b[1]*M(1,1) + b[1]^2*M(0,2)) + 
                  a[1]^3*M(3,0) + 3*a[1]^2*b[1]*M(2,1) + 3*a[1]*b[1]^2*M(1,2) + b[1]^3*M(0,3)
    elseif alpha == 1 && beta == 2 && gamma == 0 #secondMoment
        area_tt = vo[1]*vo[2]^2*M(0,0) + 
                  2*vo[1]*vo[2]*(a[2]*M(1,0) + b[2]*M(0,1)) + 
                  vo[1]*(a[2]^2*M(2,0) + 2*a[2]*b[2]*M(1,1) + b[2]^2*M(0,2)) + 
                  vo[2]^2*(a[1]*M(1,0) + b[1]*M(0,1)) + 
                  2*vo[2]*(a[1]*(a[2]*M(2,0) + b[2]*M(1,1)) + b[1]*(a[2]*M(1,1)+b[2]*M(0,2))) + 
                  a[1]*(a[2]^2*M(3,0) + 2*a[2]*b[2]*M(2,1) + b[2]^2*M(1,2)) + 
                  b[1]*(a[2]^2*M(2,1) + 2*a[2]*b[2]*M(1,2) + b[2]^2*M(0,3))
    elseif alpha == 1 && beta == 0 && gamma == 2 #secondMoment
        area_tt = vo[1]*vo[3]^2*M(0,0) + 
                  2*vo[1]*vo[3]*(a[3]*M(1,0) + b[3]*M(0,1)) + 
                  vo[1]*(a[3]^2*M(2,0) + 2*a[3]*b[3]*M(1,1) + b[3]^2*M(0,2)) + 
                  vo[3]^2*(a[1]*M(1,0) + b[1]*M(0,1)) + 
                  2*vo[3]*(a[1]*(a[3]*M(2,0) + b[3]*M(1,1)) + b[1]*(a[3]*M(1,1) + b[3]*M(0,2))) + 
                  a[1]*(a[3]^2*M(3,0) + 2*a[3]*b[3]*M(2,1) + b[3]^2*M(1,2)) + 
                  b[1]*(a[3]^2*M(2,1) + 2*a[3]*b[3]*M(1,2) + b[3]^2*M(0,3))
    elseif alpha == 1 && beta == 1 && gamma == 1 #inertiaProduct
        area_tt = vo[1]*vo[2]*vo[3]*M(0,0) + 
                  vo[1]*vo[2]*(a[3]*M(1,0) + b[3]*M(0,1)) + 
                  vo[1]*vo[3]*(a[2]*M(1,0) + b[2]*M(0,1)) + 
                  vo[1]*(a[2]*(a[3]*M(2,0) + b[3]*M(1,1)) + b[2]*(a[3]*M(1,1) + b[3]*M(0,2))) + 
                  vo[2]*vo[3]*(a[1]*M(1,0) + b[1]*M(0,1)) + 
                  vo[2]*(a[1]*(a[3]*M(2,0) + b[3]*M(1,1)) + b[1]*(a[3]*M(1,1) + b[3]*M(0,2))) + 
                  vo[3]*(a[1]*(a[2]*M(2,0) + b[2]*M(1,1)) + b[1]*(a[2]*M(1,1) + b[2]*M(0,2))) +
                  a[1]*(a[2]*(a[3]*M(3,0) + b[3]*M(2,1)) + b[2]*(a[3]*M(2,1) + b[3]*M(1,2))) + 
                  b[1]*(a[2]*(a[3]*M(2,1) + b[3]*M(1,2)) + b[2]*(a[3]*M(1,2) + b[3]*M(3,0)))
    elseif alpha == 2 && beta == 0 && gamma == 1 #inertiaProduct
        area_tt = vo[1]^2*vo[3]*M(0,0) +
                  vo[1]^2*(a[3]*M(1,0) + b[3]*M(0,1)) +
                  2*vo[1]*vo[3]*(a[1]*M(1,0) + b[1]*M(0,1)) +
                  2*vo[1]*(a[1]*(a[3]*M(2,0) + b[3]*M(1,1)) + b[1]*(a[3]*M(1,1) + b[3]*M(0,2))) +
                  vo[3]*(a[1]^2*M(2,0) + 2*a[1]*b[1]*M(1,1) + b[1]^2*M(0,2)) +
                  a[1]^2*(a[3]*M(3,0) + b[3]*M(2,1)) +
                  2*a[1]*b[1]*(a[3]*M(2,1) + b[3]*M(1,2)) +
                  b[1]^2*(a[3]*M(1,2) + b[3]*M(0,3))
    elseif alpha == 2 && beta == 1 && gamma == 0 #inertiaProduct
        area_tt = vo[1]^2*vo[2]*M(0,0) +
                  vo[1]^2*(a[2]*M(1,0) + b[2]*M(0,1)) +
                  2*vo[1]*vo[2]*(a[1]*M(1,0) + b[1]*M(0,1)) +
                  2*vo[1]*(a[1]*(a[2]*M(2,0) + b[2]*M(1,1)) + b[1]*(a[2]*M(1,1) + b[2]*M(0,2))) +
                  vo[2]*(a[1]^2*M(2,0) + 2*a[1]*b[1]*M(1,1) + b[1]^2*M(0,2)) +
                  a[1]^2*(a[2]*M(3,0) + b[2]*M(2,1)) +
                  2*a[1]*b[1]*(a[2]*M(2,1) + b[2]*M(1,2)) +
                  b[1]^2*(a[2]*M(1,2) + b[2]*M(0,3))
    else
        area_tt = s1(a, b, alpha, beta, gamma, vo)
    end
    if signedInt == true
        return area_tt * norm(c) * sign(c[3])
    else
        return area_tt * norm(c)
    end
end

""" 
	II(P::Lar.LAR, alpha::Int, beta::Int, gamma::Int, signedInt=false)
	
Basic integration function on 2D plane.
# Example  unit 3D triangle
```julia
julia> V = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
3×3 Array{Float64,2}:
 0.0  1.0  0.0
 0.0  0.0  1.0
 0.0  0.0  0.0
 
julia> FV = [[1,2,3]]
1-element Array{Array{Int64,1},1}:
 [1, 2, 3]
 
julia> P = V,FV
([0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0], Array{Int64,1}[[1, 2, 3]])
julia> Integrals.II(P, 0,0,0)
0.5
```
"""
function II(P::LAR, alpha::Int, beta::Int, gamma::Int, signedInt=false)::Float64
    V, FV = P
    partialSum = zeros(nthreads())
    @threads for i=1:length(FV)
        tau = hcat([V[:,v] for v in FV[i]]...)
        if size(tau,2) == 3
            term = TT(tau, alpha, beta, gamma, signedInt)
            if signedInt
		        @inbounds partialSum[threadid()] += term
		    else
		        @inbounds partialSum[threadid()] += abs(term)
		    end
        elseif size(tau,2) > 3
            println("ERROR: FV[$(i)] is not a triangle")
        else
            println("ERROR: FV[$(i)] is degenerate")
        end
    end    
    return sum(partialSum)
end

""" 
	III(P::LAR, alpha::Int, beta::Int, gamma::Int, signedInt::Bool=false)::Float64
Basic integration function on 3D space.
# Example # unit 3D tetrahedron
```julia
julia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
3×4 Array{Float64,2}:
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0
julia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]]
4-element Array{Array{Int64,1},1}:
 [1, 2, 4]
 [1, 3, 2]
 [4, 3, 1]
 [2, 3, 4]
julia> P = V,FV
([0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0], 
Array{Int64,1}[[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]])
julia> Integrals.III(P, 0,0,0)
0.16666666666666674
```
"""
function III(P::LAR, alpha::Int, beta::Int, gamma::Int, signedInt::Bool=false)::Float64
    V, FV = P
    partialSum = zeros(nthreads())
    @threads for i=1:length(FV)
        tau = hcat([V[:,v] for v in FV[i]]...)
        vo,va,vb = tau[:,1],tau[:,2],tau[:,3]
        a = va - vo
        b = vb - vo
        c = cross(a,b)
        term = c[1]/norm(c) * TT(tau, alpha+1, beta, gamma, signedInt)
        if signedInt
            @inbounds partialSum[threadid()] += term
        else
            @inbounds partialSum[threadid()] += abs(term)
        end
    end
    return sum(partialSum)/(alpha + 1)
end

"""
	surface(P::Lar.LAR, signedInt::Bool=false)::Float64
`surface` integral on polyhedron `P`.
# Example unit 3D triangle
```julia
julia> V = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0]
3×3 Array{Float64,2}:
 0.0  1.0  0.0
 0.0  0.0  1.0
 0.0  0.0  0.0
 
julia> FV = [[1,2,3]]
1-element Array{Array{Int64,1},1}:
 [1, 2, 3]
 
julia> P = V,FV
([0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0], Array{Int64,1}[[1, 2, 3]])
julia> Integrals.surface(P)
0.5
```
"""
function surface(P::LAR, signedInt::Bool=false)::Float64
    return II(P, 0, 0, 0, signedInt)
end

"""
	volume(P::LAR, signedInt::Bool=false)::Float64
`volume` integral on polyhedron `P`.
# Example # unit 3D tetrahedron
```julia
julia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
3×4 Array{Float64,2}:
 0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0
 0.0  0.0  0.0  1.0
julia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]]
4-element Array{Array{Int64,1},1}:
 [1, 2, 4]
 [1, 3, 2]
 [4, 3, 1]
 [2, 3, 4]
julia> P = V,FV
([0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0], 
Array{Int64,1}[[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]])
julia> Integrals.volume(P)
0.16666666666666674
```
"""
function volume(P::LAR, signedInt::Bool=false)::Float64
    return III(P, 0, 0, 0, signedInt)
end

""" 
	firstMoment(P::Lar.LAR)::Array{Float64,1}
First moments as terms of the Euler tensor. Remember that the integration algorithm is a boundary integration. Hence the model must be a boundary model. In this case, a 2-complex of triangles. 
# Example # unit 3D tetrahedron
```julia
julia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
julia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
julia> P = V,FV;
julia> Integrals.firstMoment(P)
3-element Array{Float64,1}:
 0.0416667
 0.0416667
 0.0416667
```
"""
function firstMoment(P::LAR)::Array{Float64,1}
    out = zeros(3)
    @async begin
		out[1] = III(P, 1, 0, 0)
		out[2] = III(P, 0, 1, 0)
		out[3] = III(P, 0, 0, 1)
    end
    return fetch(out)
end

""" 
	secondMoment(P::Lar.LAR)::Array{Float64,1}
Second moments as terms of the Euler tensor.
# Example # unit 3D tetrahedron
```julia
julia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
julia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
julia> P = V,FV;
julia> Integrals.secondMoment(P)
3-element Array{Float64,1}:
 0.0166667
 0.0166667
 0.0166667
```
"""
function secondMoment(P::LAR)::Array{Float64,1}
    out = zeros(3)
    @async begin
		out[1] = III(P, 2, 0, 0)
		out[2] = III(P, 0, 2, 0)
		out[3] = III(P, 0, 0, 2)
    end
    return fetch(out)
end

""" 
	inertiaProduct(P::Lar.LAR)::Array{Float64,1}
Inertia products as terms of the Euler tensor.
# Example # unit 3D tetrahedron
```julia
julia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
julia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
julia> P = V,FV;
julia> Integrals.inertiaProduct(P)
3-element Array{Float64,1}:
 0.00833333
 0.00833333
 0.00833333
```
"""
function inertiaProduct(P::LAR)::Array{Float64,1}
    out = zeros(3)
    @async begin
		out[1] = III(P, 0, 1, 1)
		out[2] = III(P, 1, 0, 1)
		out[3] = III(P, 1, 1, 0)
    end
    return fetch(out)
end

""" 
	centroid(P::Lar.LAR)::Array{Float64,1}
Barycenter or `centroid` of polyhedron `P`.
# Example # unit 3D tetrahedron
```julia
julia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
julia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
julia> P = V,FV;
julia> Integrals.centroid(P)
3-element Array{Float64,1}:
 0.25
 0.25
 0.25
```
"""
function centroid(P::LAR)::Array{Float64,1}
    return firstMoment(P)./volume(P)
end

""" 
	inertiaMoment(P::Lar.LAR)::Array{Float64,1}
Inertia moments  of polyhedron `P`.
# Example # unit 3D tetrahedron
```julia
julia> V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
julia> FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
julia> P = V,FV;
julia> Integrals.inertiaMoment(P)
3-element Array{Float64,1}:
 0.0333333
 0.0333333
 0.0333333
```
"""
function inertiaMoment(P::LAR)::Array{Float64,1}
    out = zeros(3)
    result = secondMoment(P)
    out[1] = result[2] + result[3]
    out[2] = result[3] + result[1]
    out[3] = result[1] + result[2]
    return out
end

#=
function chainAreas(V::Array{Float64,2},EV::Array{Int64,2},chains::Array{Int64,2})
	FE = [chains[:,f] for f=1:size(chains,2)]
	return chainAreas(V,EV,FE)
end
function chainAreas(V::Array{Float64,2}, EV::Array{Int64,2}, chains::Array{Array{Int64,1},1})
	if size(V,1) == 2
		V = vcat(V,zeros(1,size(V,2)))
	end	
	pivots = [EV[:,abs(chain[1])][1] for chain in chains]
	out = zeros(length(pivots))
	for k=1:length(chains)
		area = 0
		triangles = [[] for h=1:length(chains[k])]
		for h=1:length(chains[k])
			edge = chains[k][h]
			v1,v2 = EV[:,abs(edge)]
			if sign(edge) == -1
				v1,v2=v2,v1
			end
			triangles[h] = Int[pivots[k],v1,v2]
		end
		P = V,hcat(triangles...)
		out[k] = Surface(P,true)
	end
	return out
end
=#
