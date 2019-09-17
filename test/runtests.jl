using FractionalGaussianFields, AsyPlots 
using Test

# write your own tests here

zeroboundary(torus_gff(100))

P = [(1+cos(8θ)/8)*cis(θ) for θ=range(0,2π,length=500)]
loop = Path(P).points
n = 40

G = domainapprox(loop,n)
