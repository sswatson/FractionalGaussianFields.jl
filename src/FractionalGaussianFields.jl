__precompile__(true)

module FractionalGaussianFields

import Contour,
       Interpolations,
       LightGraphs,
       AsyPlots

export torus_gff,
       torus_fgf,
       gff,
       fgf,
       solvelaplace,
       laplacesurface, 
       SquareGrid,
       rectangulargrid,
       domainapprox,
       zeroboundary, 
       neighbors,
       vertices,
       edges,
       draw,
       boundary,
       flowline,
       interpolate


include("gridgraphs.jl")

"""
    torus_fgf(m,n,s,rng=GLOBAL_RNG)

    torus_fgf(n,s,rng=GLOBAL_RNG)

Return a sample of the FGF on an `m` × `n` square grid
with periodic boundary conditions
"""
function torus_fgf(m::Integer,
                  n::Integer,
                  s::Real,
                  rng=Base.Random.GLOBAL_RNG)
    h = complex(zeros(m,n))
    for j = 1:m
        for k = 1:n
            h[j,k] = 1/sqrt(2.0) *
                (j+k == 2 ? 0 :
                 randn(rng) + im*randn(rng) *
                 1/(sin((j-1)*pi/m)^2+(sin((k-1)*pi/n))^2)^(s/2))
        end
    end
    return real(n/sqrt(2)*ifft(h))
end

torus_fgf(n::Integer,s::Real,rng=Base.Random.GLOBAL_RNG) = torus_gff(n,n,s)

torus_gff(m::Integer,n::Integer,rng=Base.Random.GLOBAL_RNG) =
    torus_fgf(m,n,1,rng)

torus_gff(n::Integer,rng=Base.Random.GLOBAL_RNG) = torus_gff(n,n)

using DataStructures

function laplacian(G::SquareGrid,bvertices::Set)
    function newentry!(I,J,V,i,j,v)
        push!(I,i); push!(J,j); push!(V,v);
    end
    m,n = size(G)
    vertexmap = OrderedDict(map(reverse,enumerate([(i,j) for j=1:n
                                        for i=1:m if G.vertices[i,j]])))
    I = Int64[]
    J = Int64[]
    V = Float64[]
    for (v,k) in vertexmap
        if v in bvertices
            newentry!(I,J,V,k,k,1)
        else
            nbs = neighbors(G,v...)
            l = length(nbs)
            for nb in nbs
                j = vertexmap[nb]
                newentry!(I,J,V,k,j,-1)
                newentry!(I,J,V,k,k,1)
            end
        end
    end
    sparse(I,J,V), vertexmap
end

function fgf(G::SquareGrid,s::Real,bvertices::Set)
    Δ,vertexmap = laplacian(G,bvertices)
    h_vector = real(full(Δ)^(-s/2)*randn(size(Δ,1)))
    h = zeros(size(G))
    for (v,k) in vertexmap
        h[v...] = h_vector[k]
    end
    return h
end

gff(G::SquareGrid,bvertices::Set) =
    fgf(G,1,bvertices)

import Interpolations.interpolate
function interpolate(h::Array{<:Real,2})
    interpolate(h,
                Interpolations.BSpline(Interpolations.Linear()),
                Interpolations.OnGrid())
end

import Base.getindex
getindex(h::Interpolations.BSplineInterpolation,
         z0::Complex) = h[real(z0),imag(z0)]

"""
Find the flow line of angle θ in the field h.
In other words, solve the DE η'(t) = exp(ih(η(t))/χ + iθ)".
"""
function flowline(h::Interpolations.AbstractInterpolation,
                  z0::Complex,
                  χ::Real,
                  θ::Real;
                  δ::Real=0.01,
                  tol::Real=1e-3,
                  S::Set{Complex}=Set{Complex}())
    a,b = size(h.coefs)
    η = [z0]
    while 1.0 ≤ real(η[end]) ≤ a && 1.0 ≤ imag(η[end]) ≤ b
        push!(η, η[end] + δ * exp(im*h[η[end]]/χ + im*θ))
        w = η[end]
        for z in S
            if abs2(z-w) < tol
                return η
            end
        end
    end
    return η
end

flowline(h::Array{<:Real,2},
         z₀::Complex,
         χ::AbstractFloat,
         θ::AbstractFloat;
         kwargs...) = flowline(interpolate(h),z₀,χ,θ;kwargs...)

function zeroboundary(h::Array{<:Real,2})
    m,n = size(h)
    G = rectangulargrid(m,n)
    harmonic = solvelaplace(G,(i,j)->0,(i,j)->h[i,j])
    h - harmonic
end

end # module
