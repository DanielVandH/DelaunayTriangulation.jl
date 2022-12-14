function example_triangulation()
    p1 = @SVector[0.0, 1.0]
    p2 = @SVector[3.0, -1.0]
    p3 = @SVector[2.0, 0.0]
    p4 = @SVector[-1.0, 2.0]
    p5 = @SVector[4.0, 2.0]
    p6 = @SVector[-2.0, -1.0]
    pts = [p1, p2, p3, p4, p5, p6]
    T = Set{NTuple{3,Int64}}([
        (6, 3, 1),
        (3, 2, 5),
        (4, 1, 5),
        (4, 6, 1),
        (5, 1, 3)
    ])
    A = [
        0 0 1 1 1 1 1
        0 0 0 1 1 1 1
        1 0 0 1 0 1 0
        1 1 1 0 0 1 1
        1 1 0 0 0 1 1
        1 1 1 1 1 0 0
        1 1 0 1 1 0 0
    ]
    DG = DT.DelaunayGraph(relabel(UndirectedGraph(A), Dict(1:7 .=> 0:6)))
    adj = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
            (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
            (4, 1) => 5, (1, 5) => 4, (5, 4) => 1,
            (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
            (5, 1) => 3, (1, 3) => 5, (3, 5) => 1,
            (4, 5) => DT.BoundaryIndex, (5, 2) => DT.BoundaryIndex,
            (2, 3) => DT.BoundaryIndex, (3, 6) => DT.BoundaryIndex,
            (6, 4) => DT.BoundaryIndex
        )
    ))
    adj2v = DT.Adjacent2Vertex(Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 3), (3, 6), (6, 4)]),
        1 => Set{NTuple{2,Int64}}([(5, 4), (3, 5), (6, 3), (4, 6)]),
        2 => Set{NTuple{2,Int64}}([(5, 3)]),
        3 => Set{NTuple{2,Int64}}([(1, 6), (5, 1), (2, 5)]),
        4 => Set{NTuple{2,Int64}}([(1, 5), (6, 1)]),
        5 => Set{NTuple{2,Int64}}([(4, 1), (1, 3), (3, 2)]),
        6 => Set{NTuple{2,Int64}}([(1, 4), (3, 1)])
    ))
    all(all(DT.get_edge(adj, i, j) == k for (i, j) in DT.get_edge(adj2v, k)) for k in keys(adj2v.adjacent2vertex)) &&
        return pts, T, DG, adj, adj2v
end

function example_empty_triangulation()
    p1 = @SVector[0.0, 1.0]
    p2 = @SVector[3.0, -1.0]
    p3 = @SVector[2.0, 0.0]
    pts = [p1, p2, p3]
    T = Set{NTuple{3,Int64}}([])
    A = zeros(Int64, 0, 0)
    DG = DT.DelaunayGraph(UndirectedGraph(A))
    adj = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue, Dict{NTuple{2,Int64},Int64}()))
    adj2v = DT.Adjacent2Vertex(Dict(DT.BoundaryIndex => Set{NTuple{2,Int64}}()))
    return pts, T, DG, adj, adj2v
end

function circle_three_points(a, b, c)
    ax, ay = a
    bx, by = b
    cx, cy = c
    d1x, d1y = [by - ay, ax - bx]
    d2x, d2y = [cy - ay, ax - cx]
    k = d2x * d1y - d2y * d1x
    s1x, s1y = [(ax + bx) / 2, (ay + by) / 2]
    s2x, s2y = [(ax + cx) / 2, (ay + cy) / 2]
    ??? = d1x * (s2y - s1y) - d1y * (s2x - s1x)
    m = ??? / k
    centx, centy = [s2x + m * d2x, s2y + m * d2y]
    dx = centx - ax
    dy = centy - ay
    r = sqrt(dx^2 + dy^2)
    return [centx, centy], r
end

function slow_get_voronoi_cells(T::Ts, DG, i, pts, tri_to_idx) where {Ts}
    V = DT.triangle_type(Ts)
    nghs = DT.get_neighbour(DG, i) |> collect
    ref = DT.get_point(pts, i)
    DT.sort_boundary!(pts, nghs, ref)
    cell_idx = Int64[]
    for j in eachindex(nghs)
        _i, _j, _k = i, nghs[j], nghs[j == length(nghs) ? 1 : j + 1]
        if !DT.is_ghost_triangle(_i, _j, _k)
            ?? = construct_triangle(V, i, nghs[j], nghs[j == length(nghs) ? 1 : j + 1])
            if ?? ??? T
                ?? = DT.shift_triangle_1(??)
            end
            if ?? ??? T
                ?? = DT.shift_triangle_1(??)
            end
            push!(cell_idx, tri_to_idx[??])
        else
            push!(cell_idx, DT.BoundaryIndex)
        end
    end
    return cell_idx
end