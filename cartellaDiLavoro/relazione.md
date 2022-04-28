# Lar Integrals
<!-- https://khan.github.io/KaTeX/ -->
### Introduzione

Questo modulo implementa un metodo di integrazione finita di polinomi del tipo: $x^{\alpha}y^{\beta}z^{\gamma}$

Vengono implementate due funzioni `II` e `III` che permettono di fare rispettivamente l'integrale di superficie e di volume di polinomi del tipo specificato.

L'integrale di superficie viene calcolato facendo facendo la somma  degli integrali dei triangoli. I triangoli devono essere ottenuti triangolando opportunamente la superficie.

La funzione `TT` permette di calcolare l'integrale sul singolo triangolo.

L'integrale di volume si ottiene facilmente grazie al Teorema della Divergenza; questo teorema permette di calcolare un integrale di volume da un integrale di superficie.

È importante notare che tutti i domini, sia 2D che 3D, sono definiti in 3 dimensioni.

### Funzioni

La funzione `function M(alpha::Int, beta::Int)::Float64` calcola la seguente formula:

![formula2](./immagini_markdown/formula2.png)

Che con $\alpha=0$ e $\beta=0$ si riduce al calcolo dell'area del triangolo con vertici $w_o = (0, 0)$, $w_a = (1, 0)$ e $w_b = (0, 1)$, pari a $\frac{1}{2}$.

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

![formula10](./immagini_markdown/formula10.png)

Dato il triangolo $\tau$,  come array di array di vertici di  tre dimensioni $v_o = (x_o, y_o, z_o)$, estraiamo i vertici $v_o, v_a, v_b$. Estratti i vertici calcoliamo i vettori $a$ e $b$ e l’equazione parametrica del piano di inclusione del triangolo vale $p = v_o + u * a + v * b$.  

Possiamo  notare come $s4$ corrisponda all’ultima sommatoria e l’ultimo termine della sommatoria corrisponda alla funzione `M`(con  $u = (h+k+m)–(i+j+l)$  e $v=(i+j+l)$), $s3$ alla penultima sommatoria, $s2$ alla terz’ultima sommatoria e $s1$ alle prime tre sommatorie. Alla fine viene eseguito il prodotto vettoriale fra $a$ e $b$.


Le funzioni:
* `function II(P::LAR, alpha::Int, beta::Int, gamma::Int, signedInt=false)::Float64`
* `function III(P::LAR, alpha::Int, beta::Int, gamma::Int)::Float64`

Permettono di calcolare rispettivamente l'integrale doppio e l'integrale triplo di un polinomio che ha come dominio un poligono; l'integrale viene calcolato attraverso la somma delle aree dei triangoli (funzione `TT`) ottenuti dalla suddivisione della sua superficie.

### Esempi

Un esempio di calcolo della superficie può essere il seguente:

<center><img src="./immagini_markdown/quadrato.png" width="1000"></center>

La cui superficie vale: 90.0

Come riferimento per il calcolo del volume si è preso il modello 3D *Stanford bunny*:

<center><img src="./immagini_markdown/bunny.png" width="1000"></center>

Il cui volume vale: -1.8818908358435506e-5

