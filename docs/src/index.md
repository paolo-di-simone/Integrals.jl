# Integrals.jl

Questo modulo implementa un metodo di integrazione finita di polinomi. 

## Interfacce principali
- II(P::LAR, alpha::Int, beta::Int, gamma::Int, signedInt=false)::Float64
- III(P::LAR, alpha::Int, beta::Int, gamma::Int, signedInt::Bool=false)::Float64
- surface(P::LAR, signedInt::Bool=false)::Float64
- volume(P::LAR, signedInt::Bool=false)::Float64
- firstMoment(P::LAR)::Array{Float64,1}
- secondMoment(P::Lar.LAR)::Array{Float64,1}
- inertiaProduct(P::LAR)::Array{Float64,1}
- centroid(P::LAR)::Array{Float64,1}
- inertiaMoment(P::LAR)::Array{Float64,1}


## Documentazione funzioni
```@docs
Integrals.M
```

```@docs
Integrals.TT
```

```@docs
Integrals.II
```

```@docs
Integrals.III
```

```@docs
Integrals.surface
```

```@docs
Integrals.volume
```

```@docs
Integrals.firstMoment
```

```@docs
Integrals.secondMoment
```

```@docs
Integrals.inertiaProduct
```

```@docs
Integrals.centroid
```

```@docs
Integrals.inertiaMoment
```


