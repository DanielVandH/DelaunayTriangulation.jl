"""
    add_triangle!(T::Ts, i, j, r, k) where Ts

Adds the triangle `(i, j, r)` into `T`, deleting `(i, j, k)`.
"""
function add_triangle!(T::Ts, i, j, r, k) where {Ts}
    V = triangle_type(Ts)
    newT = construct_triangle(V, i, j, r)
    oldT = construct_triangle(V, i, j, k)
    add_triangle!(T, newT)
    delete_triangle!(T, oldT)
    return nothing
end
"""
    add_triangle!(adj::Adjacent, i, j, r)

Adds the triangle `(i, j, r)` into `T`, deleting `(i, j, k)`. This method updates 
the adjacent map `adj` accordingly.
"""
function add_triangle!(adj::Adjacent, i, j, r)
    add_edge!(adj, i, j, r)
    add_edge!(adj, j, r, i)
    add_edge!(adj, r, i, j)
    return nothing
end
"""
    add_triangle!(adj::Adjacent2Vertex, i, j, r)

Adds the triangle `(i, j, r)` into `T`, deleting `(i, j, k)`. This method updates 
the adjacent-to-vertex map `adj2v` accordingly.
"""
function add_triangle!(adj2v::Adjacent2Vertex, i, j, r, k)
    add_edge!(adj2v, i, j, r)
    add_edge!(adj2v, j, r, i)
    add_edge!(adj2v, r, i, j)
    delete_edge!(adj2v, i, j, k)
    delete_edge!(adj2v, j, k, i)
    delete_edge!(adj2v, k, i, j)
    return nothing
end
"""
    add_triangle!(DG::DelaunayGraph, i, j, r, k; delete_adjacent_neighbours = true)

Adds the triangle `(i, j, r)` into `T`, deleting `(i, j, k)`. This method updates 
the graph representation `DG` accordingly. If `(i, j, k)` should all remain 
neighbours, set `delete_adjacent_neighbours = false`.
"""
function add_triangle!(DG::DelaunayGraph, i, j, r, k; delete_adjacent_neighbours=true)
    add_neighbour!(DG, r, i, j)
    delete_adjacent_neighbours && delete_neighbour!(DG, k, i, j)
    return nothing
end
"""
    add_triangle!(HG::HistoryGraph, V, i, j, r, k)

Adds the triangle `(i, j, r)` into `T`, deleting `(i, j, k)`. This method updates 
the history graph `HG` accordingly. `V` is the triangle type.
"""
function add_triangle!(HG::HistoryGraph, V, i, j, r, k)
    Tᵢⱼₖ = construct_triangle(V, i, j, k)
    Tᵢⱼᵣ = construct_triangle(V, i, j, r)
    add_triangle!(HG, Tᵢⱼᵣ)
    add_edge!(HG, Tᵢⱼₖ, Tᵢⱼᵣ)
    return nothing
end
"""
    add_triangle!(T, adj::Adjacent, adj2v::Adjacent2Vertex, DG::DelaunayGraph, i, j, r, k; delete_adjacent_neighbours = true)

Adds the triangle `(i, j, r)` into the triangulation, assuming the `r`th point is to the left 
of the directed edge `(i, j)`. and `k = get_edge(adj, i, j)`. If `(i, j, k)` should all remain 
neighbours, set `delete_adjacent_neighbours = false`.
"""
function add_triangle!(T::Ts, adj::Adjacent, adj2v::Adjacent2Vertex,
    DG::DelaunayGraph, HG::HistoryGraph,
    i, j, r, k; delete_adjacent_neighbours=true) where Ts
    add_triangle!(T, i, j, r, k)
    add_triangle!(adj, i, j, r)
    add_triangle!(adj2v, i, j, r, k)
    add_triangle!(DG, i, j, r, k; delete_adjacent_neighbours)
    V = triangle_type(Ts)
    add_triangle!(HG, V, i, j, r, k)
    return nothing
end

"""
    split_point!(T, HG::HistoryGraph, adj, adj2v, DG, Tᵢⱼₖ::V, r) where V

Given a triangulation `T`, adds the `r`th point of the point set into the triangulation, assuming 
that `r` is in the interior of the triangle `Tᵢⱼₖ`.

# Arguments 
- `T`: The current set of triangles defining the triangulation.
- `HG`: The point location data structure.
- `adj`: The adjacency list.
- `adj2v`: The adjacent-to-vertex list.
- `DG`: The vertex-neighbour data structure.
-` Tᵢⱼₖ`: The triangle that the `r`th point is inside of. Must be positively oriented.
- `r`: The index of the point in the original point set that is being introduced.

# Outputs 
`T`, `HG`, `adj`, `adj2v`, and `DG` are all updated in-place.
"""
function split_triangle!(T, HG::HistoryGraph, adj, adj2v, DG, Tᵢⱼₖ::V, r) where {V}
    i, j, k = indices(Tᵢⱼₖ)
    add_triangle!(T, adj, adj2v, DG, HG, i, j, r, k; delete_adjacent_neighbours=false)
    add_triangle!(T, adj, adj2v, DG, HG, j, k, r, i; delete_adjacent_neighbours=false)
    add_triangle!(T, adj, adj2v, DG, HG, k, i, r, j; delete_adjacent_neighbours=false)
    return nothing
end