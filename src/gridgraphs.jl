
mutable struct SquareNeighbors
    E::Bool
    N::Bool
    W::Bool
    S::Bool
end

struct SquareGrid
    vertices::Array{Bool,2}
    neighbors::Array{SquareNeighbors,2}
    location::Tuple{<:Real,<:Real}
    spacing::Real
end

import Base.size
Base.size(G::SquareGrid) = size(G.vertices)

import Base.show
function Base.show(io::IO,G::SquareGrid)
    n = sum(G.vertices)
    E = sum([nb.N+nb.E for nb in G.neighbors])
    vert = n == 1 ? "vertex" : "vertices"
    edge = E == 1 ? "edge" : "edges"
    print(io,"SquareGrid($n $vert, $E $edge)")
end

function vertices(G::SquareGrid)
    m,n = size(G)
    [(i,j) for i=1:m for j=1:n if G.vertices[i,j]]
end

function edges(G::SquareGrid)
    m,n = size(G)
    ([((i,j),(i,j+1)) for i=1:m,j=1:n-1 if G.neighbors[i,j].N] ∪
     [((i,j),(i+1,j)) for i=1:m-1,j=1:n if G.neighbors[i,j].E])
end


import Base.position
function position(G::SquareGrid,i,j)
    return G.location .+ G.spacing.*(i-1,j-1)
end

function rectangulargrid(m::Integer,n::Integer;location=(0,0),spacing=1)
    V = [true for i=1:m,j=1:n]
    N = [SquareNeighbors(true,true,true,true) for i=1:m,j=1:n]
    for i=1:m
        N[i,1].S = false
        N[i,end].N = false
    end
    for j=1:n
        N[1,j].W = false
        N[end,j].E = false
    end
    return SquareGrid(V,N,location,spacing)
end

const STEP = Dict(
    :E => (1,0),
    :N => (0,1),
    :W => (-1,0),
    :S => (0,-1)
)

const OPP = Dict(
    :E => :W,
    :W => :E,
    :N => :S,
    :S => :N
)

function neighbors(G::SquareGrid,i,j)
    [(i,j) .+ STEP[s] for s in (:E,:N,:W,:S)
                            if getfield(G.neighbors[i,j],s)]
end

function between(a::Real,b::Real)
    if a == b
        return 1:0 # empty range
    elseif a < b
        return ceil(Int,a):floor(Int,b)
    else
        return ceil(Int,b):floor(Int,a)
    end
end

function domainapprox(loop::Array{<:AsyPlots.Vec2,1},N::Integer)
    mx = minimum(v.x for v in loop)
    my = minimum(v.y for v in loop)
    Mx = maximum(v.x for v in loop)
    My = maximum(v.y for v in loop)
    location = (mx,my)
    δ = max(Mx-mx,My-my)/(N-1)
    xvals = mx:δ:Mx
    yvals = my:δ:My
    m,n = map(length,(xvals,yvals))
    V = [AsyPlots.iswellinside(AsyPlots.Vec2(xvals[i],yvals[j]),
                loop;epsilon=1e-3) for i=1:m,j=1:n]
    northEdgesHit = NTuple{2,Int64}[]
    eastEdgesHit = NTuple{2,Int64}[]
    for i=1:length(loop)
        P,Q = i < length(loop) ? loop[i:i+1] : [loop[end],loop[1]]
        a,b,c,d = [t/δ for t in (P.x-mx,P.y-my,Q.x-mx,Q.y-my)]
        for i = between(a,c)
            j = floor(Int,(d-b)/(c-a)*(i-a)+b)
            push!(northEdgesHit,(i+1,j+1))
        end
        for j = between(b,d)
            i = floor(Int,(c-a)/(d-b)*(j-b)+a)
            push!(eastEdgesHit,(i+1,j+1))
        end
    end
    northHitSet, eastHitSet = map(Set,(northEdgesHit,eastEdgesHit))
    N = [SquareNeighbors(false,false,false,false) for i=1:m,j=1:n]
    for i=1:m
        for j=1:n
            for dir in (:E,:N)
                nb = (i,j) .+ STEP[dir]
                if (nb[1] ≤ m && nb[2] ≤ n &&
                       V[i,j] && V[nb...]  &&
                        ~(dir == :N && (i,j) in northHitSet) &&
                        ~(dir == :E && (i,j) in eastHitSet))
                    setfield!(N[i,j],dir,true)
                    setfield!(N[nb...],OPP[dir],true)
                end
            end
        end
    end
    return SquareGrid(V,N,location,δ)
end

function isboundary(G::SquareGrid,i,j)
    m,n = size(G)
    if ~G.vertices[i,j]
        return false
    elseif i == 1 || i == m || j == 1 || j == n
        return true
    elseif ~all(G.vertices[i+k,j+l] for (k,l) in values(STEP))
        return true
    else
        return false
    end
end

function boundary(G::SquareGrid)
    m,n = size(G)
    [(i,j) for i=1:m,j=1:n if isboundary(G,i,j)]
end

function solvelaplace(G::SquareGrid,
                      Δvals::Function,
                      bvals::Function,
                      bvertices::Set=Set(boundary(G)))
    m,n = size(G)
    function newentry!(I,J,V,i,j,v)
        push!(I,i); push!(J,j); push!(V,v);
    end
    # number the vertices in the graph, for purposes of
    # setting up a linear system
    vertexmap = Dict(map(reverse,enumerate([(i,j) for j=1:n
                                for i=1:m if G.vertices[i,j]])))
    # I, J are index lists and V a value list for a sparse matrix
    # used to construct the linear system describing the Laplace
    # equation
    I = Int64[]
    J = Int64[]
    V = Float64[]
    # b is a vector representing the right-hand side of the
    # linear system describing the Laplace equation
    b = zeros(length(vertexmap))
    for (v,k) in vertexmap
        if v in bvertices
            # ensure boundary values are satisfied
            newentry!(I,J,V,k,k,1)
            b[k] = bvals(v...)
        else
            nbs = neighbors(G,v...)
            l = length(nbs)
            for nb in nbs
                # row for Laplacian of vertex v
                j = vertexmap[nb]
                newentry!(I,J,V,k,j,1/l)
                newentry!(I,J,V,k,k,-1/l)
            end
            b[k] = Δvals(v...)
        end
    end
    A = SparseArrays.sparse(I,J,V)
    soln = LinearAlgebra.lu(A) \ b
    h = zeros(m,n)
    for (k,v) in vertexmap
        h[k...] = soln[v]
    end
    return h
end

function draw(G::SquareGrid;kwargs...)
    grlist = AsyPlots.GraphicElement2D[]
    for v in vertices(G)
        push!(grlist,AsyPlots.Point(position(G,v...);kwargs...))
    end
    for e in edges(G)
        push!(grlist,AsyPlots.Path([position(G,v...) for v in e];kwargs...))
    end
    AsyPlots.Plot(grlist)
end


function laplacesurface(G::SquareGrid,
                        Δvals::Function,
                        bvals::Function,
                        bvertices::Set=Set(boundary(G)))
    h = solvelaplace(G,Δvals,bvals,bvertices)
    AsyPlots.Plot(AsyPlots.Surface(h,
                                  clip=G.vertices,
                                  spline=false),
                  ignoreaspect=true,zmax=2)
end
