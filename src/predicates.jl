function ExactPredicates.orient(pᵢx, pᵢy, pⱼx, pⱼy, pₖx, pₖy)
    return orient((pᵢx, pᵢy), (pⱼx, pⱼy), (pₖx, pₖy))
end

"""
    isoriented(T, pts)

Tests if the triangle `T` is positively oriented, returning
`1` if the triangle is positively oriented, `0` 
if the triangle is degenerate (i.e. the points are all 
collinear), and `-1` if the triangle is negatively 
oriented.

The argument `pts` should be a collection of points 
such that the `i`th vertex of `T` corresponds to 
`get_point(pts, i)`.

The check is exact, making use of `ExactPredicates.orient` (unless 
you define a new method for it for your point type, obviously).
"""
function isoriented(T, pts)
    i, j, k = indices(T)
    pᵢ, pⱼ, pₖ = get_point(pts, i, j, k)
    if !is_ghost_triangle(i, j, k)
        return orient(pᵢ, pⱼ, pₖ)
    else
        return orient(pₖ, pⱼ, pᵢ) # Ghost triangles have the 0 point represented at a centroid which is inwards, so flip the orientation
    end
end

"""
    isincircle(pts, i, j, k, ℓ)

Tests if the `ℓ`th point is in the circle through the 
`i`th, `j`th, and `k` points (as obtained via `get_point(pts, i)`, etc.).
Returns `1` if the point is inside, `0` if the point is on the circle,
and `-1` if the point is outside. It is assumed that the 
`i`th, `j`th`, and `k`th points are provided counter-clockwise. 

The check is exact, making use of `ExactPredicates.incircle` (unless 
you define a new method for it for your point type, obviously).

Checks are made for ghost triangles, which are triangles of the form `(i, j, $BoundaryIndex)`. We say 
that a point is in the circumcircle of `(i, j, $BoundaryIndex)` if it is to the left of `(i, j)`.
"""
function isincircle(pts, i::I, j::I, k::I, ℓ::I) where {I} # clean this up later to reduce redundant code
    pᵢ, pⱼ, pₖ, pₗ = get_point(pts, i, j, k, ℓ)
    if i == I(BoundaryIndex)
        in_open_halfplane = I(isleftofline(pts, j, k, ℓ))
        if in_open_halfplane == I(0)
            on_edge = I(-sameside(pₗ, pⱼ, pₖ))
            if on_edge == I(1)
                return I(0)
            else
                return I(-1)
            end
        else
            return in_open_halfplane
        end
        return in_open_halfplane == I(0) ? I(-sameside(pₗ, pⱼ, pₖ)) : in_open_halfplane # If in_in_open_halfplane = 0, then the point is collinear with the edge, so we need to check if the point is on the edge or away from it.
    elseif j == I(BoundaryIndex)
        in_open_halfplane = I(isleftofline(pts, k, i, ℓ))
        if in_open_halfplane == I(0)
            on_edge = I(-sameside(pₗ, pₖ, pᵢ))
            if on_edge == I(1)
                return I(0)
            else
                return I(-1)
            end
        else
            return in_open_halfplane
        end
    elseif k == I(BoundaryIndex)
        in_open_halfplane = I(isleftofline(pts, i, j, ℓ))
        if in_open_halfplane == I(0)
            on_edge = I(-sameside(pₗ, pᵢ, pⱼ))
            if on_edge == I(1)
                return I(0)
            else
                return I(-1)
            end
        else
            return in_open_halfplane
        end
    end
    return I(incircle(pᵢ, pⱼ, pₖ, pₗ))
end
"""
    isincircle(T, pts, ℓ) 

Tests if the `ℓ`th point in `pts` is inside the circumdisk of the 
triangle `T`. 
"""
function isincircle(T, pts, ℓ)
    i, j, k = indices(T)
    return isincircle(pts, i, j, k, ℓ)
end

"""
    isleftofline(pts, i, j, k)

Tests if the `k`th point if left of the line from the `i`th point 
to the `j`th point (as obtained via `get_point(pts, i)`, etc.). 
Returns `1` if the point is to the left, `0` if the points are 
collinear, and `-1` if the point is to the right.

The check is exact, making use of `ExactPredicates.orient` (unless 
you define a new method for it for your point type, obviously).
"""
function isleftofline(pts, i::I, j::I, k::I) where {I}
    if i == I(LowerRightBoundingIndex) && j == I(LowerLeftBoundingIndex)
        return I(-1)
    elseif i == I(LowerRightBoundingIndex) && j == I(UpperBoundingIndex)
        return I(1)
    elseif i == I(LowerLeftBoundingIndex) && j == I(LowerRightBoundingIndex)
        return I(1)
    elseif i == I(LowerLeftBoundingIndex) && j == I(UpperBoundingIndex)
        return I(-1)
    elseif i == I(UpperBoundingIndex) && j == I(LowerLeftBoundingIndex)
        return I(1)
    elseif i == I(UpperBoundingIndex) && j == I(LowerRightBoundingIndex)
        return I(-1)
    elseif i == I(BoundaryIndex) || j == I(BoundaryIndex)
        j, i = i, j # If the line is ℓ₀ⱼ, then this is parallel with the line ℓⱼₚ, where ₚ is the centroid. But 0 maps to p (`get_point(pts, 0)` is the centroid), so let's just swap the coordinates. Similarly, for the line ℓᵢ₀, we just consider ℓₚᵢ, so again we swap.
    end
    pᵢ, pⱼ, pₖ = get_point(pts, i, j, k)
    return I(orient(pᵢ, pⱼ, pₖ))
end

"""
    isintriangle(T, pts, ℓ)

Tests if the `ℓ`th point of `pts` is in the triangle `T`. Returns `1`
if the point is inside, `0` if a point is on an edge of `T`, and `-1` 
if the point is outside. It is assumed that `T` is positively oriented. 

The check is exact, making use of `isleftofline`. 
"""
function isintriangle(T, pts, ℓ::I) where {I}
    i, j, k = indices(T)
    if i < I(FirstPointIndex) && j < I(FirstPointIndex) && k < I(FirstPointIndex)
        return I(1)
    end
    is_left_of_edge_ij = isleftofline(pts, i, j, ℓ)
    is_left_of_edge_jk = isleftofline(pts, j, k, ℓ)
    is_left_of_edge_ki = isleftofline(pts, k, i, ℓ)
    return isintriangle(
        is_left_of_edge_ij,
        is_left_of_edge_jk,
        is_left_of_edge_ki
    )
end
function isintriangle(e1::I, e2::I, e3::I) where {I}
    isinexterior = e1 == I(-1) || e2 == I(-1) || e3 == I(-1)
    if isinexterior
        return I(-1)
    end
    if e1 == I(0) || e2 == I(0) || e3 == I(0)
        return I(0)
    end
    isininterior = e1 == I(1) && e2 == I(1) && e3 == I(1)
    if isininterior
        return I(1)
    else
        return I(-1)
    end
end
function isintriangle(p::P, q::P, r::P, s::P) where {P}
    o1 = orient(p, q, s)
    o2 = orient(q, r, s)
    o3 = orient(r, p, s)
    if any(==(-1), (o1, o2, o3))
        return -1
    elseif any(==(0), (o1, o2, o3))
        return 0
    else
        return 1
    end
end
"""
    find_edge(T, pts, ℓ)

Given a triangle `T` and a set of points `pts`, with the `ℓ`th point of `pts` 
on an edge of `T`, find the edge of `T` that the point is on.
"""
function find_edge(T, pts, ℓ)
    i, j, k = indices(T)
    isleftofline(pts, i, j, ℓ) == 0 && return (i, j)
    isleftofline(pts, j, k, ℓ) == 0 && return (j, k)
    isleftofline(pts, k, i, ℓ) == 0 && return (k, i)
    throw("The point $p is not on an edge of $T.")
end

"""
    edge_on_bounding_triangle(i, j)

Returns true if `(i, j)` is an edge of bounding triangle $(BoundingTriangle).
"""
edge_on_bounding_triangle(i, j) = i < FirstPointIndex && j < FirstPointIndex

"""
    islegal(i, j, adj, pts)

Tests if the edge `(i, j)` is legal. `adj` is the adjacent map of the 
triangulation, and `pts` is the point set.
"""
function islegal(i, j, adj, pts)
    edge_on_bounding_triangle(i, j) && return true
    is_boundary_edge(i, j, adj) && return true
    is_boundary_edge(j, i, adj) && return true
    !edge_exists(i, j, adj) && return true
    !edge_exists(j, i, adj) && return true
    k = get_edge(adj, i, j)
    ℓ = get_edge(adj, j, i)
    return islegal(i, j, k, ℓ, pts)
end
function islegal(i::I, j::I, k::I, ℓ::I, pts) where {I}
    incirc = isincircle(pts, i, j, k, ℓ)
    return incirc == I(-1) || incirc == I(0)
end

"""
    isinconvexhull(pts, BN, q)

Tests if the point `q` is inside the convex hull `pts[BN]`, assuming that the `pts[BN]` are in counter-clockwise 
order. Returns `1` if the point is inside, `0` if the point is on, or `-1` if the point is 
outside of the convex hull.
"""
function isinconvexhull(pts, BN, q)
    for i in eachindex(BN)
        u = BN[i]
        j = i == lastindex(BN) ? firstindex(BN) : i + 1
        v = BN[j]
        r, s = get_point(pts, u, v)
        o = orient(r, s, q)
        if o == 0
            return 0
        elseif o == -1
            return -1
        end
    end
    return 1
end