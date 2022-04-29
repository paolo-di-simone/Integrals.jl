# Lar Integrals

<!-- https://khan.github.io/KaTeX/ -->

### Introduzione

Questo modulo implementa un metodo di integrazione finita di polinomi del tipo: <img alt="x^{\alpha}y^{\beta}z^{\gamma}" src="https://render.githubusercontent.com/render/math?math=x%5E%7B%5Calpha%7Dy%5E%7B%5Cbeta%7Dz%5E%7B%5Cgamma%7D" style="transform: translateY(20%);" />

Vengono implementate due funzioni `II` e `III` che permettono di fare rispettivamente l'integrale di superficie e di volume di polinomi del tipo specificato.

L'integrale di superficie viene calcolato facendo facendo la somma  degli integrali dei triangoli. I triangoli devono essere ottenuti triangolando opportunamente la superficie.

La funzione `TT` permette di calcolare l'integrale sul singolo triangolo.

L'integrale di volume si ottiene facilmente grazie al Teorema della Divergenza; questo teorema permette di calcolare un integrale di volume da un integrale di superficie.

È importante notare che tutti i domini, sia 2D che 3D, sono definiti in 3 dimensioni.

### Descrizione funzioni di integrazione

La funzione `function M(alpha::Int, beta::Int)::Float64` calcola la seguente formula:

<img alt="II^{\alpha\beta} = \frac{1}{\alpha + 1}\sum_{h=0}^{\alpha + 1}{\binom{\alpha + 1}{h}\frac{(-1)^h}{h+\beta+1}}" src="https://render.githubusercontent.com/render/math?math=II%5E%7B%5Calpha%5Cbeta%7D%20%3D%20%5Cfrac%7B1%7D%7B%5Calpha%20%2B%201%7D%5Csum_%7Bh%3D0%7D%5E%7B%5Calpha%20%2B%201%7D%7B%5Cbinom%7B%5Calpha%20%2B%201%7D%7Bh%7D%5Cfrac%7B%28-1%29%5Eh%7D%7Bh%2B%5Cbeta%2B1%7D%7D" style="transform: translateY(20%);" />

Che con <img alt="\alpha=0" src="https://render.githubusercontent.com/render/math?math=%5Calpha%3D0" style="transform: translateY(20%);" /> e <img alt="\beta=0" src="https://render.githubusercontent.com/render/math?math=%5Cbeta%3D0" style="transform: translateY(20%);" /> si riduce al calcolo dell'area del triangolo con vertici <img alt="w_o = (0, 0)" src="https://render.githubusercontent.com/render/math?math=w_o%20%3D%20%280%2C%200%29" style="transform: translateY(20%);" />, <img alt="w_a = (1, 0)" src="https://render.githubusercontent.com/render/math?math=w_a%20%3D%20%281%2C%200%29" style="transform: translateY(20%);" /> e <img alt="w_b = (0, 1)" src="https://render.githubusercontent.com/render/math?math=w_b%20%3D%20%280%2C%201%29" style="transform: translateY(20%);" />, pari a <img alt="\frac{1}{2}" src="https://render.githubusercontent.com/render/math?math=%5Cfrac%7B1%7D%7B2%7D" style="transform: translateY(20%);" />.

La seguente funzione:

function TT(tau::Array{Float64,2}, alpha::Int, beta::Int, gamma::Int, signedInt::Bool=false)

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

Permette di calcolare l'integrale di un triangolo implementando la seguente formula:

<img alt="II_{\tau}^{\alpha\beta\gamma}=II^{uv}|a \times b|\sum_{h=0}^{\alpha}\binom{\alpha}{h}x_0^{\alpha-h}\sum_{k=0}^{\beta}\binom{\beta}{k}y_0^{\beta-k}\sum_{m=0}^{\gamma}\binom{\gamma}{m}z_0^{\gamma-m}\sum_{i=0}^{h}\binom{h}{i}a_x^{h-i}b_x^i\sum_{j=0}^{k}\binom{k}{j}a_v^{k-j}b_y^j\sum_{l=0}^{m}\binom{m}{l}a_z^{m-l}b_z^l" src="https://render.githubusercontent.com/render/math?math=II_%7B%5Ctau%7D%5E%7B%5Calpha%5Cbeta%5Cgamma%7D%3DII%5E%7Buv%7D%7Ca%20%5Ctimes%20b%7C%5Csum_%7Bh%3D0%7D%5E%7B%5Calpha%7D%5Cbinom%7B%5Calpha%7D%7Bh%7Dx_0%5E%7B%5Calpha-h%7D%5Csum_%7Bk%3D0%7D%5E%7B%5Cbeta%7D%5Cbinom%7B%5Cbeta%7D%7Bk%7Dy_0%5E%7B%5Cbeta-k%7D%5Csum_%7Bm%3D0%7D%5E%7B%5Cgamma%7D%5Cbinom%7B%5Cgamma%7D%7Bm%7Dz_0%5E%7B%5Cgamma-m%7D%5Csum_%7Bi%3D0%7D%5E%7Bh%7D%5Cbinom%7Bh%7D%7Bi%7Da_x%5E%7Bh-i%7Db_x%5Ei%5Csum_%7Bj%3D0%7D%5E%7Bk%7D%5Cbinom%7Bk%7D%7Bj%7Da_v%5E%7Bk-j%7Db_y%5Ej%5Csum_%7Bl%3D0%7D%5E%7Bm%7D%5Cbinom%7Bm%7D%7Bl%7Da_z%5E%7Bm-l%7Db_z%5El" style="transform: translateY(20%);" />

Dato il triangolo <img alt="\tau" src="https://render.githubusercontent.com/render/math?math=%5Ctau" style="transform: translateY(20%);" />,  come array di array di vertici di  tre dimensioni <img alt="v_o = (x_o, y_o, z_o)" src="https://render.githubusercontent.com/render/math?math=v_o%20%3D%20%28x_o%2C%20y_o%2C%20z_o%29" style="transform: translateY(20%);" />, estraiamo i vertici <img alt="v_o, v_a, v_b" src="https://render.githubusercontent.com/render/math?math=v_o%2C%20v_a%2C%20v_b" style="transform: translateY(20%);" />. Estratti i vertici calcoliamo i vettori <img alt="a" src="https://render.githubusercontent.com/render/math?math=a" style="transform: translateY(20%);" /> e <img alt="b" src="https://render.githubusercontent.com/render/math?math=b" style="transform: translateY(20%);" /> e l’equazione parametrica del piano di inclusione del triangolo vale <img alt="p = v_o + u * a + v * b" src="https://render.githubusercontent.com/render/math?math=p%20%3D%20v_o%20%2B%20u%20%2a%20a%20%2B%20v%20%2a%20b" style="transform: translateY(20%);" />.

Possiamo  notare come <img alt="s4" src="https://render.githubusercontent.com/render/math?math=s4" style="transform: translateY(20%);" /> corrisponda all’ultima sommatoria e l’ultimo termine della sommatoria corrisponda alla funzione `M`(con  <img alt="u = (h+k+m)–(i+j+l)" src="https://render.githubusercontent.com/render/math?math=u%20%3D%20%28h%2Bk%2Bm%29%E2%80%93%28i%2Bj%2Bl%29" style="transform: translateY(20%);" />  e <img alt="v=(i+j+l)" src="https://render.githubusercontent.com/render/math?math=v%3D%28i%2Bj%2Bl%29" style="transform: translateY(20%);" />), <img alt="s3" src="https://render.githubusercontent.com/render/math?math=s3" style="transform: translateY(20%);" /> alla penultima sommatoria, <img alt="s2" src="https://render.githubusercontent.com/render/math?math=s2" style="transform: translateY(20%);" /> alla terz’ultima sommatoria e <img alt="s1" src="https://render.githubusercontent.com/render/math?math=s1" style="transform: translateY(20%);" /> alle prime tre sommatorie. Alla fine viene eseguito il prodotto vettoriale fra <img alt="a" src="https://render.githubusercontent.com/render/math?math=a" style="transform: translateY(20%);" /> e <img alt="b" src="https://render.githubusercontent.com/render/math?math=b" style="transform: translateY(20%);" />.

Le funzioni:
* `function II(P::LAR, alpha::Int, beta::Int, gamma::Int, signedInt=false)::Float64`
* `function III(P::LAR, alpha::Int, beta::Int, gamma::Int)::Float64`
Permettono di calcolare rispettivamente l'integrale doppio e l'integrale triplo di un polinomio che ha come dominio un poligono; l'integrale viene calcolato attraverso la somma delle aree dei triangoli (funzione `TT`) ottenuti dalla suddivisione della sua superficie.

### Esempi

Un esempio di calcolo della superficie può essere il seguente (da cui sono stati omessi gli output):
       
        julia> function makeSurface(s)
               
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
       
       julia> V, VV, EV, FV = makeSurface(10)
       
       julia> model = (V, (VV, EV, FV))
       
       julia> P = V, FV
       
       julia> Plasm.view(Plasm.hpc_exploded(model)(1.1,1.1,1.1))
       
       julia> surface(P)
       
       90.0

<div align="center"><img src="./immagini_markdown/quadrato.png" width="700"></div>

Come riferimento per il calcolo del volume si è preso il modello 3D *Stanford bunny*:

   julia> V, EV, FV = obj2lar("bunny.obj")
    
    julia> EV = EV[1]
    
    julia> FV = FV[1]
    
    julia> VV = [[k] for k=1:size(V,2)]
    
    julia> model = (V, (VV, EV, FV))
    
    julia> P = V, FV
    
    julia> Plasm.view(V, FV)
    
    julia> volume(P)
    -1.8818908358435506e-5

<div align="center"><img src="./immagini_markdown/bunny.png" width="700"></div>

### Interfacce principali

Le interfacce principali del modulo sono:
#### `Lar.surface`
#### `Lar.volume`
#### `Lar.centroid`
#### `Lar.inertiaProduct`
#### `Lar.inertiaMoment`
