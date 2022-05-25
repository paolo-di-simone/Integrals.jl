using LinearAlgebraicRepresentation
using Plasm
Lar = LinearAlgebraicRepresentation
include("obj2lar.jl")
include("integrals.jl")

function obj2lar(path)
    vs = Array{Float64, 2}(undef, 0, 3)
    edges = Array{Array{Array{Int, 1}, 1}, 1}()
    faces = Array{Array{Array{Int, 1}, 1}, 1}()
    push!(edges, Array{Array{Int, 1}, 1}[])
    push!(faces, Array{Array{Int, 1}, 1}[])
    g = 1

    open(path, "r") do fd
        for line in eachline(fd)
            elems = split(line)
            if length(elems) > 0
                if elems[1] == "v"
                    
                    x = parse(Float64, elems[2])
                    y = parse(Float64, elems[3])
                    z = parse(Float64, elems[4])

                    vs = [vs; x y z]
                    
                elseif elems[1] == "f"  
                    
                    v1 = parse(Int, split(elems[2], "/")[1])
                    v2 = parse(Int, split(elems[3], "/")[1])
                    v3 = parse(Int, split(elems[4], "/")[1])
                    
                    e1 = sort([v1, v2])
                    e2 = sort([v2, v3])
                    e3 = sort([v3, v1])
                    

                    push!(edges[g], e1)
                    push!(edges[g], e2)
                    push!(edges[g], e3)

                    push!(faces[g], sort([v1, v2, v3]))

                elseif elems[1] == "g"

                    g += 1
                    push!(edges, Array{Array{Int, 1}, 1}[])
                    push!(faces, Array{Array{Int, 1}, 1}[])
                    
                end
            end
        end
    end

    return convert(Lar.Points, vs'), edges[1:end], faces[1:end]
end

V, EV, FV = obj2lar("bunny.obj")
EV = EV[1]
FV = FV[1]
VV = [[k] for k=1:size(V,2)]
model = (V, (VV, EV, FV))
P = V, FV
# Plasm.view(V, FV)
volume(P)

function makeSurface(s)
    W,(WW,EW,FW) = Lar.cuboid([s,s], true)
    main_square = (W,EW)
    Q,(QQ,EQ,FQ) = Lar.cuboid([1,1], true)
    square = (Q,EQ)
    args = []
    append!(args, [main_square])
    append!(args, [square])
    for i in 1:(s - 1)
        append!(args, [Lar.t(1,1)])
        append!(args, [square])
    end
    model = Lar.Struct(args)
    V,EV = Lar.struct2lar(model)
    VV = [[k] for k=1:size(V,2)]
    FV = Lar.triangulate2d(V, EV)
    a = zeros(1, size(V, 2))
    V = vcat(V, a)
    return V, VV, EV, FV
end
V, VV, EV, FV = makeSurface(10)
model = (V, (VV, EV, FV))
P = V, FV
# Plasm.view(Plasm.hpc_exploded(model)(1.1,1.1,1.1))
surface(P)
