"""
    edge_radius_ratio(T, pts)

Given a triangle `T` and a set of points `pts`, computes the edge-radius triangle of `T`, 
defined by `ρ = R/ℓ`, where `R = abc/(4A)` is the triangle's circumradius, the edge lengths 
are `a`, `b`, and `c`, `A` is the triangle's area, and `ℓ` is the shortest edge length.
"""
function edge_radius_ratio(T, pts)
    F = number_type(pts)
    shortest_length = typemax(F)
    circumradius = one(F) # abc/(4A)
    triangle_area = zero(F)
    for (u, v) in edges(T)
        pᵤ, pᵥ = get_point(pts, u, v)
        xᵤ, yᵤ = getx(pᵤ), gety(pᵤ)
        xᵥ, yᵥ = getx(pᵥ), gety(pᵥ)
        ℓ = sqrt((xᵤ - xᵥ)^2 + (yᵤ - yᵥ)^2)
        circumradius *= ℓ
        triangle_area += 0.5(xᵤ * yᵥ - xᵥ * yᵤ) # shoelace formula 
        shortest_length = ℓ < shortest_length ? ℓ : shortest_length
    end
    circumradius /= (4triangle_area)
    return circumradius / shortest_length
end

"""
    encroaches_upon(q, e, pts)

Given a point `q`, an edge `e`, and a set of points `pts` for getting the coordinates 
of the vertices of `e`, checks if `q` encroaches upon `e`.

(A point `q` is said to encroach upon `e` if, provided `q` is none of the endpoints 
of `e`, it is inside the closed diametric ball of `e`.)
"""
function encroaches_upon(q, e, pts)
    mid_x, mid_y = edge_midpoint(e, pts)
    ball_radius = 0.5edge_length(e, pts)
    qx, qy = getx(q), gety(q)
    dist_to_midpoint_sq = (qx - mid_x)^2 + (qy - mid_y)^2
    return dist_to_midpoint_sq < ball_radius^2
end

"""
    segment_encroached(e, adj, pts)

Tests if the edge `e` is encroached, using the fact that an edge `e` is encroached 
if `e ∉ Del(S)`, and if `e ∈ Del(S)` then `e` is encroached if there exists a triangle 
with e as an edge that has an angle opposite e exceeding 90°.

(A segment is encroached if there is a point, not including the edge's endpoints, inside 
its closed diametric ball.))
"""
function segment_encroached(e, adj, pts)
    if !edge_exists(e, adj)
        return true
    else
        # If e ∈ Del(S), it is encroached iff there exists a triangle 
        # with e as an edge that has an angle exceeding 90° opposite e
        i, j = e
        k = get_edge(adj, i, j)
        ℓ = get_edge(adj, j, i)
        pᵢ, pⱼ = get_point(pts, i, j)
        pᵢx, pᵢy = getx(pᵢ), gety(pᵢ)
        pⱼx, pⱼy = getx(pⱼ), gety(pⱼ)
        ij_length = sqrt((pᵢx - pⱼx)^2 + (pᵢy - pⱼy)^2)
        if k ≠ BoundaryIndex && k ≠ DefaultAdjacentValue
            pₖ = get_point(pts, k)
            pₖx, pₖy = getx(pₖ), gety(pₖ)
            jk_length = sqrt((pⱼx - pₖx)^2 + (pⱼy - pₖy)^2)
            ik_length = sqrt((pᵢx - pₖx)^2 + (pᵢy - pₖy)^2)
            if opposite_angle_is_obtuse(jk_length, ik_length, ij_length)
                return true
            end
        end
        if ℓ ≠ BoundaryIndex && ℓ ≠ DefaultAdjacentValue
            pₗ = get_point(pts, ℓ)
            pₗx, pₗy = getx(pₗ), gety(pₗ)
            jℓ_length = sqrt((pⱼx - pₗx)^2 + (pⱼy - pₗy)^2)
            iℓ_length = sqrt((pᵢx - pₗx)^2 + (pᵢy - pₗy)^2)
            if opposite_angle_is_obtuse(jℓ_length, iℓ_length, ij_length)
                return true
            end
        end
        return false
    end
end

