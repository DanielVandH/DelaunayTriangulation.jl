############################################
##
## MAIN TRIANGULATION FUNCTIONS
##
############################################
"""
    initialise_triangulation(pts; 
    IntegerType=Int64,
    TriangleType=Triangle{IntegerType},
    EdgeType=Edge{IntegerType})

This function returns the initial form of a [`Triangulation`](@ref) data structure, storing 
the points in `pts` (converted into a `Points` type). You can specify custom integer, 
triangle, and edge types using the keywords `IntegerType`, `TriangleType`, and 
`EdgeType`, respectively.
"""
function initialise_triangulation(pts;
    IntegerType::Type{I}=Int64,
    TriangleType=Triangle{IntegerType},
    EdgeType=Edge{IntegerType}) where {I}
    # The data structures
    root = TriangleType(I(LowerRightBoundingIndex),
        I(UpperBoundingIndex),
        I(LowerLeftBoundingIndex))
    T = Triangles{I,TriangleType}(Set{TriangleType}([root]))
    HG = HistoryDAG{I,TriangleType}()
    adj = Adjacent{I,EdgeType}()
    adj2v = Adjacent2Vertex{I,EdgeType}()
    DG = DelaunayGraph{I}()
    # Add the root to the DAG
    add_triangle!(HG, root)
    # Add the initial adjacencies 
    add_edge!(adj, I(LowerRightBoundingIndex), I(UpperBoundingIndex), I(LowerLeftBoundingIndex))
    add_edge!(adj, I(UpperBoundingIndex), I(LowerLeftBoundingIndex), I(LowerRightBoundingIndex))
    add_edge!(adj, I(LowerLeftBoundingIndex), I(LowerRightBoundingIndex), I(UpperBoundingIndex))
    add_edge!(adj, I(LowerRightBoundingIndex), I(LowerLeftBoundingIndex), I(BoundaryIndex))
    add_edge!(adj, I(LowerLeftBoundingIndex), I(UpperBoundingIndex), I(BoundaryIndex))
    add_edge!(adj, I(UpperBoundingIndex), I(LowerRightBoundingIndex), I(BoundaryIndex))
    add_edge!(adj2v, I(LowerLeftBoundingIndex), I(LowerRightBoundingIndex), I(UpperBoundingIndex))
    add_edge!(adj2v, I(LowerRightBoundingIndex), I(UpperBoundingIndex), I(LowerLeftBoundingIndex))
    add_edge!(adj2v, I(UpperBoundingIndex), I(LowerLeftBoundingIndex), I(LowerRightBoundingIndex))
    # Add the initial neighbours 
    add_point!(DG, I(LowerLeftBoundingIndex), I(LowerRightBoundingIndex), I(UpperBoundingIndex))
    add_neighbour!(DG, I(LowerLeftBoundingIndex), I(LowerRightBoundingIndex), I(UpperBoundingIndex))
    add_neighbour!(DG, I(LowerRightBoundingIndex), I(LowerLeftBoundingIndex), I(UpperBoundingIndex))
    add_neighbour!(DG, I(UpperBoundingIndex), I(LowerLeftBoundingIndex), I(LowerRightBoundingIndex))
    return Triangulation(adj, adj2v, DG, HG, T, Points(pts), root)
end

"""
    locate_triangle(HG::HistoryDAG, pts, p, init=find_root(HG; method=:rng))

Given the point location data structure `HG` and a set of `pts`, finds the triangle in 
the current triangulation such that `p` is in its interior. The point location starts at `init`.
The function is recursive, and returns a tuple `(tri, flag)`:

    - `tri`: This is the triangle that `p` is in.
    - `flag`: If `flag == 0`, then `p` is on an edge of `tri`. Otherwise, it is in the open interior.
"""
function locate_triangle(HG::HistoryDAG, pts, p, init=find_root(HG; method=:rng))
    if out_deg(HG, init) == 0
        return init, intriangle(init, pts, p)
    end
    out = out_neighbors(HG, init)
    for T in out
        intriangle(T, pts, p) ≥ 0 && return locate_triangle(HG, pts, p, T)
    end
    throw("Failed to find triangle.")
end

"""
    add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)

Given a triangulation `T`, adds the `r`th point of the point set into the triangulation.

# Arguments 
- `T`: The current triangulation.
- `HG`: The point location data structure.
- `adj`: The adjacency list.
- `adj2v`: The adjacent-to-vertex list.
- `DG`: The vertex-neighbour data structure.
-` Tᵢⱼₖ`: The triangle that the `r`th point is inside of. Must be positively oriented.
- `r`: The index of the point in the original point set that is being introduced.

# Outputs 
`T`, `HG`, `adj`, `adj2v`, and `DG` are all updated in-place.
"""
function add_point!(T::Triangles{I,V}, HG::HistoryDAG,
    adj::Adjacent, adj2v::Adjacent2Vertex,
    DG::DelaunayGraph, Tᵢⱼₖ::V, r) where {I,V<:AbstractTriangle{I}}
    i, j, k = Tᵢⱼₖ # The triangle to be split into three
    delete_triangle!(T, Tᵢⱼₖ) # Now that we've split the triangle, we can remove the triangle
    T₁, T₂, T₃ = V(i, j, r), V(j, k, r), V(k, i, r) # New triangles to add. Note that these triangles are all positively oriented.
    add_triangle!(T, T₁, T₂, T₃) # The three new triangles
    add_triangle!(HG, T₁, T₂, T₃) # Add the new triangles into DAG
    add_edge!(HG, Tᵢⱼₖ, T₁, T₂, T₃) # Add edges from the old triangle to the new triangles
    add_triangle!(adj, T₁, T₂, T₃) # Add the new edges into the adjacency list
    update_after_insertion!(adj2v, i, j, k, r)
    add_neighbour!(DG, r, i, j, k)
    add_neighbour!(DG, i, r)
    add_neighbour!(DG, j, r)
    add_neighbour!(DG, k, r)
    return nothing
end

"""
    flip_edge!(T, HG, adj, adj2v, DG, i, j, k, r)

Performs an edge flip, flipping the edge `(i, j)` into the edge `(k, r)`.

# Arguments
- `T`: The current triangulation.
- `HG`: The point location data structure.
- `adj`: The adjacency list.
- `adj2v`: The adjacent-to-vertex list.
- `DG`: The vertex-neighbour data structure.
- `i, j`: The current edge.
- `k, r`: Indices for the points the edge is flipped onto.

It is assumed that `(i, k, j)` and `(i, j, r)` are positively oriented triangles.

# Outputs 
`T`, `HG`, `adj`, `adj2v`, and `DG` are all updated in-place.
"""
function flip_edge!(T::Triangles{I,V}, HG::HistoryDAG,
    adj::Adjacent, adj2v::Adjacent2Vertex, DG::DelaunayGraph, i, j, k, r) where {I,V<:AbstractTriangle{I}}
    # The old triangles
    Tᵢₖⱼ = V(i, k, j)
    Tᵢⱼᵣ = V(i, j, r)
    delete_triangle!(T, Tᵢₖⱼ, Tᵢⱼᵣ)
    delete_edge!(adj, i, j)
    delete_neighbour!(DG, i, j) #delete_neighbour!(DG, j, i)
    # The new triangles 
    Tᵣₖⱼ = V(r, k, j)
    Tᵣᵢₖ = V(r, i, k)
    # Add the new triangles to the data structure
    add_triangle!(T, Tᵣₖⱼ, Tᵣᵢₖ)
    add_triangle!(HG, Tᵣₖⱼ, Tᵣᵢₖ)
    add_triangle!(adj, Tᵣₖⱼ, Tᵣᵢₖ)
    update_after_flip!(adj2v, i, j, k, r)
    # Connect the new triangles to the replaced triangles in the DAG
    add_edge!(HG, Tᵢₖⱼ, Tᵣₖⱼ, Tᵣᵢₖ)
    add_edge!(HG, Tᵢⱼᵣ, Tᵣₖⱼ, Tᵣᵢₖ)
    # Add the new neighbours 
    add_neighbour!(DG, r, k) # add_neighbour!(𝒟𝒢, k, r)
    return nothing
end

"""
    legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
    
Legalises the edge `(i, j)` if it is illegal.

# Arguments 
- `T`: The current triangulation.
- `HG`: The point location data structure.
- `adj`: The adjacency list.
- `adj2v`: The adjacent-to-vertex list.
- `DG`: The vertex-neighbour data structure.
- `i, j`: The edge to make legal. Nothing happens if `is_legal(i, j, 𝒜, pts)`.
- `r`: The point being added into the triangulation. 
- `pts`: The point set of the triangulation.

# Outputs 
`T`, `HG`, `adj`, `adj2v`, and `DG` are all updated in-place.
"""
function legalise_edge!(T::Triangles, HG::HistoryDAG,
    adj::Adjacent, adj2v::Adjacent2Vertex,
    DG::DelaunayGraph, i, j, r, pts::Points)
    if !islegal(i, j, adj, pts)
        e = get_edge(adj, j, i)
        flip_edge!(T, HG, adj, adj2v, DG, i, j, e, r)
        legalise_edge!(T, HG, adj, adj2v, DG, i, e, r, pts)
        legalise_edge!(T, HG, adj, adj2v, DG, e, j, r, pts)
    end
    return nothing
end

"""
    remove_bounding_triangle!(DT::Triangulation)
    remove_bounding_triangle!(T, adj, adj2v, DG)

Remove the bounding triangle from the triangulation, where 

- `T`: The current triangulation.
- `adj`: The adjacency list.
- `adj2v`: The adjacent-to-vertex list.
- `DG`: The vertex-neighbour data structure.

These structures are updated in-place.
"""
function remove_bounding_triangle!(T::Triangles{I,V}, adj::Adjacent,
    adj2v::Adjacent2Vertex, DG::DelaunayGraph) where {I,V<:AbstractTriangle{I}}
    for w in BoundingTriangle
        neighbours = get_edge(adj2v, w)
        for (u, v) in neighbours # (u, v, w) is a triangle..
            delete_edge!(adj, w, u; protect_boundary=false)
            delete_edge!(adj, w, v; protect_boundary=false)
            delete_edge!(adj2v, u, v, w)
            delete_edge!(adj2v, v, w, u)
            if u ≥ FirstPointIndex && v ≥ FirstPointIndex # This can only be a boundary edge
                add_edge!(adj2v, BoundaryIndex, u, v)
                add_edge!(adj, u, v, BoundaryIndex)
            end
            delete_triangle!(T, V(u, v, w))
        end
        delete_point!(DG, w)
        delete_point!(adj2v, w)
    end
    return nothing
end
@doc (@doc remove_bounding_triangle!(::Triangles{I,V}, ::Adjacent,
    ::Adjacent2Vertex, ::DelaunayGraph) where {I,V<:AbstractTriangle{I}})
function remove_bounding_triangle!(DT::Triangulation)
    remove_bounding_triangle!(triangles(DT),
        adjacent(DT), adjacent2vertex(DT), graph(DT))
    return nothing
end

"""
    add_point!(DT::Triangulation, r::Integer)
    add_point!(T::Triangles, HG::HistoryDAG,
     adj::Adjacent, adj2v::Adjacent2Vertex, 
     DG::DelaunayGraph, root, pts, r)

Adds `get_points(pts, r)` to the triangulation. 
"""
function add_point!(T::Triangles, HG::HistoryDAG,
    adj::Adjacent, adj2v::Adjacent2Vertex,
    DG::DelaunayGraph, root, pts, r)
    pᵣ = get_point(pts, r)
    Tᵢⱼₖ, interior_flag = locate_triangle(HG, pts, pᵣ, root)
    i, j, k = Tᵢⱼₖ
    if interior_flag == 1
        add_point!(T, HG, adj, adj2v, DG, Tᵢⱼₖ, r)
        legalise_edge!(T, HG, adj, adj2v, DG, i, j, r, pts)
        legalise_edge!(T, HG, adj, adj2v, DG, j, k, r, pts)
        legalise_edge!(T, HG, adj, adj2v, DG, k, i, r, pts)
    else
        eᵢⱼ, pₖ = locate_edge(pᵣ, 𝒯ᵢⱼₖ)
        𝒯ᵢₗⱼ = adjacent_triangles(𝒯, 𝒯ᵢⱼₖ, eᵢⱼ)
        pₗ = select_adjacent_vertex(𝒯, eᵢⱼ, 𝒯ᵢₗⱼ)
        add_edges!(𝒯, 𝒟, pᵣ, pₖ)
        add_edges!(𝒯, 𝒟, pᵣ, pₗ)
        legalise_edge!(𝒯, 𝒟, 𝒯ᵢₗⱼ, pᵣ, pᵢ, pₗ)
        legalise_edge!(𝒯, 𝒟, 𝒯ᵢₗⱼ, pᵣ, pₗ, pⱼ)
        legalise_edge!(𝒯, 𝒟, 𝒯ᵢⱼₖ, pᵣ, pⱼ, pₖ)
        legalise_edge!(𝒯, 𝒟, 𝒯ᵢⱼₖ, pᵣ, pₖ, pᵢ)
    end
end
@doc (@doc add_point!(::Triangles, ::HistoryDAG, ::Adjacent, ::Adjacent2Vertex, ::DelaunayGraph, ::Any, ::Any, ::Any))
function add_point!(DT::Triangulation, r::Integer)
    add_point!(triangles(DT), history(DT), adjacent(DT),
        adjacent2vertex(DT), graph(DT), root(DT), points(DT), r)
    return nothing
end

"""
    add_point!(DT::Triangulation, p)

Adds the point `p` into the triangulation.
"""
function add_point!(DT::Triangulation, p)
    add_point!(points(DT), Point(p))
    add_point!(DT, lastindex(points(DT)))
    return nothing
end

"""
    triangulate(pts; shuffle_pts=true, trim=true)

Computes the Delaunay triangulation of the points in `pts` using randomised incremental 
insertion. The points are shuffled in-place, but this shuffling can be disabled by 
setting `shuffle_pts=false`. The bounding triangle of the triangulation can be retained 
by setting `trim=true`.

You can use custom integer, triangle, and edge types by using the keyword arguments 
`IntegerType`, `TriangleType`, and `EdgeType`. See their definitions in 
[`initialise_triangulation`](@ref).
"""
function triangulate(pts; shuffle_pts=true, trim=true, method = :berg,
    IntegerType=Int64,TriangleType=Triangle{IntegerType},EdgeType=Edge{IntegerType})
    DT = initialise_triangulation(pts;IntegerType,TriangleType,EdgeType)
    shuffle_pts && @views shuffle!(points(DT))
    for r in eachindex(points(DT))
        add_point!(DT, r)
    end
    trim && remove_bounding_triangle!(DT)
    return DT
end