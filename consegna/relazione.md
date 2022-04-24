# Lar Integrals

Questo modulo implementa un metodo di integrazione finita di polinomi del tipo:

$x^\alphay^\betaz^\gamma$

	
Vengono implementate due funzioni II e III che permettono di fare
rispettivamente l'integrale di superficie e di volume di polinomi del tipo
specificato.

L'integrale di superficie viene calcolato facendo facendo la somma 
degli integrali dei triangoli. I triangoli devono essere ottenuti triangolando 
opportunamente la superficie.

La funzione TT permette appunto di calcolare l'integrale sul singolo triangolo.

L'integrale di volume si ottiene facilmente grazie al Teorema della Divergenza;
questo teorema permette di trasformare un integrale di volume in un integrale di
superficie.

È importante notare che tutti i domini, sia 2D che 3D, sono definiti in 3 dimensioni.
