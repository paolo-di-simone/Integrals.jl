
using LinearAlgebra

const Points = Matrix

const Cells = Array{Array{Int,1},1}

const LAR = Union{ Tuple{Points, Cells},Tuple{Points, Cells, Cells} }


function M(alpha::Int, beta::Int)::Float64
    a = 0
    for l=0:(alpha + 1)
        a += binomial(alpha+1,l) * (-1)^l/(l+beta+1)
    end
    return a/(alpha + 1)
end

# Quanto diminuisce la velocità di calcolo se invece di un'approssimazione Float64 si usa una Float128?
# Documentare molto i test
# Utilizzare librerie ad-hoc per l'hardware sottostante


# P = modello poliedrale, nel quale si va ad integrare (descritto mediante la sua superficie fatta di triangoli)
# alpha, beta e gamma esponenti di un monomio in x, y e z (si sta integrando un monomio del tipo x^alpha * y ^ beta * z ^ gamma su un insieme di triangoli)

# Funzione principale di integrazione sul singolo triangolo
# Un triangolo una matrice 3x3 di coordinate
# Da quello che ho capito viene fatto l'integrale di superficie del triangolo, quindi viene calcolata l'area del triangolo (e non un volume)
# L'algoritmo porta ad un errore di approssimazione alla terza cifra decimale (trascurabile?)

# In generale l'algoritmo non l'ho capito
function TT(tau::Array{Float64,2}, alpha::Int, beta::Int, gamma::Int, signedInt::Bool=false)
	# tau = triangolo
	vo,va,vb = tau[:,1],tau[:,2],tau[:,3]
	a = va - vo
	b = vb - vo
	s1 = 0.0
	for h=0:alpha
		for k=0:beta
			for m=0:gamma
				s2 = 0.0
				for i=0:h 
					s3 = 0.0
					for j=0:k
						s4 = 0.0
						for l=0:m
							s4 += binomial(m,l) * a[3]^(m-l) * b[3]^l * M(h+k+m-i-j-l, i+j+l)
						end
						s3 += binomial(k,j) * a[2]^(k-j) * b[2]^j * s4
					end
					s2 += binomial(h,i) * a[1]^(h-i) * b[1]^i * s3;
				end
				s1 += binomial(alpha,h) * binomial(beta,k) * binomial(gamma,m) * vo[1]^(alpha-h) * vo[2]^(beta-k) * vo[3]^(gamma-m) * s2
			end
		end
	end
	c = cross(a,b)
	if signedInt == true
		return s1 * norm(c) * sign(c[3])
	else
		return s1 * norm(c)
	end	
end

# Da quello che ho capito abbbiamo V e FV che rappresentano (vedi PDF):
#		- FV due celle 
#		- V uno celle

function II(P::LAR, alpha::Int, beta::Int, gamma::Int, signedInt=false)::Float64
    w = 0
    V, FV = P
    if typeof(FV) == Array{Int64,2}
    	FV = [FV[:,k] for k=1:size(FV,2)]
    end
    for i=1:length(FV)
        tau = hcat([V[:,v] for v in FV[i]]...)
        if size(tau,2) == 3
        	term = TT(tau, alpha, beta, gamma, signedInt)
        	if signedInt
        		w += term
        	else
        		w += abs(term)
        	end
        elseif size(tau,2) > 3
        	println("ERROR: FV[$(i)] is not a triangle")
        else
        	println("ERROR: FV[$(i)] is degenerate")
        end
    end    
    return w
end

function III(P::LAR, alpha::Int, beta::Int, gamma::Int)::Float64
    w = 0
    V, FV = P
    for i=1:length(FV)
        tau = hcat([V[:,v] for v in FV[i]]...)
        vo,va,vb = tau[:,1],tau[:,2],tau[:,3]
        a = va - vo
        b = vb - vo
        c = cross(a,b)
        w += c[1]/norm(c) * TT(tau, alpha+1, beta, gamma)
    end
    return w/(alpha + 1)
end




function surface(P::LAR, signedInt::Bool=false)::Float64
    return II(P, 0, 0, 0, signedInt)
end


function volume(P::LAR)::Float64
    return III(P, 0, 0, 0)
end

function firstMoment(P::LAR)::Array{Float64,1}
    out = zeros(3)
    out[1] = III(P, 1, 0, 0)
    out[2] = III(P, 0, 1, 0)
    out[3] = III(P, 0, 0, 1)
    return out
end

function secondMoment(P::LAR)::Array{Float64,1}
    out = zeros(3)
    out[1] = III(P, 2, 0, 0)
    out[2] = III(P, 0, 2, 0)
    out[3] = III(P, 0, 0, 2)
    return out
end

function inertiaProduct(P::LAR)::Array{Float64,1}
    out = zeros(3)
    out[1] = III(P, 0, 1, 1)
    out[2] = III(P, 1, 0, 1)
    out[3] = III(P, 1, 1, 0)
    return out
end


function centroid(P::LAR)::Array{Float64,1}
	return firstMoment(P)./volume(P)
end


function inertiaMoment(P::LAR)::Array{Float64,1}
    out = zeros(3)
    result = secondMoment(P)
    out[1] = result[2] + result[3]
    out[2] = result[3] + result[1]
    out[3] = result[1] + result[2]
    return out
end


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

