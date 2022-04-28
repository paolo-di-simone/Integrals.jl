using LinearAlgebraicRepresentation
using Plasm
Lar = LinearAlgebraicRepresentation
include("obj2lar.jl")
include("integrals.jl")

V, EV, FV = obj2lar("bunny.obj")
EV = EV[1]
FV = FV[1]
VV = [[k] for k=1:size(V,2)]
model = (V, (VV, EV, FV))
P = V, FV
Plasm.view(V, FV)
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
Plasm.view(Plasm.hpc_exploded(model)(1.1,1.1,1.1))
surface(P)
