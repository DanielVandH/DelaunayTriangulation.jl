@testset "Can we identify boundary edges" begin
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
    for (k, S) in adjacent2vertex(adj2v)
        for (i, j) in S
            if k == DT.BoundaryIndex
                @test DT.is_boundary_edge((i, j), adj)
                @test DT.is_boundary_edge(i, j, adj)
                @test DT.is_boundary_edge((i, j), adj2v)
                @test DT.is_boundary_edge(i, j, adj2v)
            else
                @test !DT.is_boundary_edge((i, j), adj)
                @test !DT.is_boundary_edge(i, j, adj)
                @test !DT.is_boundary_edge((i, j), adj2v)
                @test !DT.is_boundary_edge(i, j, adj2v)
            end
        end
    end
end

@testset "Can we check edge validity?" begin
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
    for (i, j) in edges(adj)
        @test DT.edge_exists(i, j, adj)
    end
    for _ in 1:1000
        i, j = abs.(rand(Int64, 3))
        @test !DT.edge_exists(i, j, adj)
    end
end

@testset "Is choose_uvw working correctly?" begin
    for _ in 1:99818
        i, j, k = abs.(rand(Int64, 3))
        @test DT.choose_uvw(true, false, false, i, j, k) == (i, j, k)
        @test DT.choose_uvw(false, true, false, i, j, k) == (j, k, i)
        @test DT.choose_uvw(false, false, true, i, j, k) == (k, i, j)
    end
end

@testset "Can we correctly clear the empty keys in the adjacent map?" begin
    p1 = (5.0, 6.0)
    p2 = (9.0, 6.0)
    p3 = (13.0, 5.0)
    p4 = (10.38, 0.0)
    p5 = (12.64, -1.69)
    p6 = (2.0, -2.0)
    p7 = (3.0, 4.0)
    p8 = (7.5, 3.53)
    p9 = (4.02, 1.85)
    p10 = (4.26, 0.0)
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]
    Random.seed!(928881)
    T, adj, adj2v, DG, HG = DT.triangulate_berg(pts; trim_empty_features=false)
    p11 = (6.0, 2.5)
    push!(pts, p11)
    r = 11
    V = DT.locate_triangle(T, pts, r)
    i, j, k = indices(V)
    DT.delete_triangle!(i, j, k, T, adj, adj2v, DG)
    @test adjacent(adj) == DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (8, 3) => 2, (3, 2) => 8, (2, 8) => 3,
            (1, 7) => 8, (7, 8) => 1, (8, 1) => 7,
            (8, 7) => 9, (7, 9) => 8, (9, 8) => 7,
            (4, 10) => 6, (10, 6) => 4, (6, 4) => 10,
            (4, 6) => 5, (6, 5) => 4, (5, 4) => 6,
            (8, 4) => 3, (4, 3) => 8, (3, 8) => 4,
            (9, 6) => 10, (6, 10) => 9, (10, 9) => 6,
            (1, 8) => 2, (8, 2) => 1, (2, 1) => 8,
            (9, 7) => 6, (7, 6) => 9, (6, 9) => 7,
            (5, 3) => 4, (3, 4) => 5, (4, 5) => 3,
            (8, 10) => 4, (10, 4) => 8, (4, 8) => 10,
            (1, 2) => DT.BoundaryIndex,
            (2, 3) => DT.BoundaryIndex,
            (3, 5) => DT.BoundaryIndex,
            (5, 6) => DT.BoundaryIndex,
            (6, 7) => DT.BoundaryIndex,
            (7, 1) => DT.BoundaryIndex,
            (((8, -1), (2, 7), (10, 7), (7, 5), (-3, 8), (-3, 10),
                (-3, 2), (5, 10), (2, -1), (-3, 5), (8, 5)) .=> DT.DefaultAdjacentValue)...
        )
    )
    DT.clear_empty_keys!(adj)
    @test adjacent(adj) == DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (8, 3) => 2, (3, 2) => 8, (2, 8) => 3,
            (1, 7) => 8, (7, 8) => 1, (8, 1) => 7,
            (8, 7) => 9, (7, 9) => 8, (9, 8) => 7,
            (4, 10) => 6, (10, 6) => 4, (6, 4) => 10,
            (4, 6) => 5, (6, 5) => 4, (5, 4) => 6,
            (8, 4) => 3, (4, 3) => 8, (3, 8) => 4,
            (9, 6) => 10, (6, 10) => 9, (10, 9) => 6,
            (1, 8) => 2, (8, 2) => 1, (2, 1) => 8,
            (9, 7) => 6, (7, 6) => 9, (6, 9) => 7,
            (5, 3) => 4, (3, 4) => 5, (4, 5) => 3,
            (8, 10) => 4, (10, 4) => 8, (4, 8) => 10,
            (1, 2) => DT.BoundaryIndex,
            (2, 3) => DT.BoundaryIndex,
            (3, 5) => DT.BoundaryIndex,
            (5, 6) => DT.BoundaryIndex,
            (6, 7) => DT.BoundaryIndex,
            (7, 1) => DT.BoundaryIndex,
        )
    )
end

@testset "Can we correctly compare collections of triangles?" begin
    T = [(1, 5, 7), (10, 5, 3), (1, 2, 3), (3, 2, 1), (7, 10, 0)]
    V = [(1, 5, 7), (10, 5, 3), (1, 2, 3), (3, 2, 1), (7, 10, 0)]
    @test DT.compare_triangle_sets(T, V)
    @test DT.compare_triangle_sets(V, T)
    V = [(1, 5, 7), (10, 5, 3), (1, 2, 3), (3, 2, 1)]
    @test !DT.compare_triangle_sets(T, V)
    @test !DT.compare_triangle_sets(V, T)
    V = [(1, 5, 7), (10, 5, 3), (1, 2, 3), (1, 3, 2), (0, 7, 10)]
    @test DT.compare_triangle_sets(T, V)
    @test DT.compare_triangle_sets(V, T)
    V = [(1, 5, 7), (10, 5, 3), (1, 2, 3), (3, 2, 1), (7, 6, 3)]
    @test !DT.compare_triangle_sets(T, V)
    @test !DT.compare_triangle_sets(V, T)
    V = [(5, 1, 7), (10, 5, 3), (1, 2, 3), (3, 2, 1), (7, 6, 3)]
    @test !DT.compare_triangle_sets(T, V)
    @test !DT.compare_triangle_sets(V, T)
    T = Set{NTuple{3,Int64}}([(1, 5, 7), (10, 5, 3), (1, 2, 3), (3, 2, 1), (7, 10, 0)])
    V = Set{NTuple{3,Int64}}([(1, 5, 7), (10, 5, 3), (1, 2, 3), (3, 2, 1), (7, 10, 0)])
    @test DT.compare_triangle_sets(T, V)
    @test DT.compare_triangle_sets(V, T)
    V = Set{NTuple{3,Int64}}([(1, 5, 7), (10, 5, 3), (1, 2, 3), (3, 2, 1)])
    @test !DT.compare_triangle_sets(T, V)
    @test !DT.compare_triangle_sets(V, T)
    V = Set{NTuple{3,Int64}}([(1, 5, 7), (10, 5, 3), (1, 2, 3), (1, 3, 2), (0, 7, 10)])
    @test DT.compare_triangle_sets(T, V)
    @test DT.compare_triangle_sets(V, T)
    V = Set{NTuple{3,Int64}}([(1, 5, 7), (10, 5, 3), (1, 2, 3), (3, 2, 1), (7, 6, 3)])
    @test !DT.compare_triangle_sets(T, V)
    @test !DT.compare_triangle_sets(V, T)
    V = Set{NTuple{3,Int64}}([(5, 1, 7), (10, 5, 3), (1, 2, 3), (3, 2, 1), (7, 6, 3)])
    @test !DT.compare_triangle_sets(T, V)
    @test !DT.compare_triangle_sets(V, T)
end

@testset "Can we correctly compare triangulations?" begin
    p1 = (5.0, 6.0)
    p2 = (9.0, 6.0)
    p3 = (13.0, 5.0)
    p4 = (10.38, 0.0)
    p5 = (12.64, -1.69)
    p6 = (2.0, -2.0)
    p7 = (3.0, 4.0)
    p8 = (7.5, 3.53)
    p9 = (4.02, 1.85)
    p10 = (4.26, 0.0)
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]
    Random.seed!(928881)
    T, adj, adj2v, DG, HG = DT.triangulate_berg(pts)
    @test DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, T, adj, adj2v, DG)
    @test !DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, Set{NTuple{3,Int64}}([(5, 1, 7), (10, 5, 3), (1, 2, 3), (3, 2, 1), (7, 6, 3)]), adj, adj2v, DG)
    adj2 = deepcopy(adj)
    adj2.adjacent[(6, 10)] = 120
    @test !DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, T, adj2, adj2v, DG)
    adj2v2 = deepcopy(adj2v)
    DT.delete_point!(adj2v2, 3)
    @test !DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, T, adj, adj2v2, DG)
    DG2 = deepcopy(DG)
    DT.delete_point!(DG2, 7)
    @test !DT.compare_unconstrained_triangulations(T, adj, adj2v, DG, T, adj, adj2v, DG2)
end

@testset "Can we correctly validate the relationship between Adjacent and Adjacent2Vertex?" begin
    adj = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (1, 2) => 3, (2, 3) => 1, (3, 1) => 2,
            (3, 2) => 4, (2, 4) => 3, (4, 3) => 2,
            (3, 4) => 5, (4, 5) => 3, (5, 3) => 4,
            (4, 2) => 6, (2, 6) => 4, (6, 4) => 2,
            (5, 4) => 6, (4, 6) => 5, (6, 5) => 4,
            (2, 1) => DT.BoundaryIndex, (1, DT.BoundaryIndex) => 2, (DT.BoundaryIndex, 2) => 1,
            (1, 3) => DT.BoundaryIndex, (3, DT.BoundaryIndex) => 1, (DT.BoundaryIndex, 1) => 3,
            (3, 5) => DT.BoundaryIndex, (5, DT.BoundaryIndex) => 3, (DT.BoundaryIndex, 3) => 5,
            (5, 6) => DT.BoundaryIndex, (6, DT.BoundaryIndex) => 5, (DT.BoundaryIndex, 5) => 6,
            (6, 2) => DT.BoundaryIndex, (2, DT.BoundaryIndex) => 6, (DT.BoundaryIndex, 6) => 2
        )))
    adj2v = DT.Adjacent2Vertex(
        Dict(
            DT.BoundaryIndex => Set{NTuple{2,Int64}}([(1, 3), (3, 5), (5, 6), (6, 2), (2, 1)]),
            1 => Set{NTuple{2,Int64}}([(2, 3), (3, DT.BoundaryIndex), (DT.BoundaryIndex, 2)]),
            2 => Set{NTuple{2,Int64}}([(3, 1), (4, 3), (6, 4), (1, DT.BoundaryIndex), (DT.BoundaryIndex, 6)]),
            3 => Set{NTuple{2,Int64}}([(1, 2), (2, 4), (4, 5), (DT.BoundaryIndex, 1), (5, DT.BoundaryIndex)]),
            4 => Set{NTuple{2,Int64}}([(3, 2), (2, 6), (6, 5), (5, 3)]),
            5 => Set{NTuple{2,Int64}}([(3, 4), (4, 6), (6, DT.BoundaryIndex), (DT.BoundaryIndex, 3)]),
            6 => Set{NTuple{2,Int64}}([(4, 2), (5, 4), (DT.BoundaryIndex, 5), (2, DT.BoundaryIndex)])
        )
    )
    @test DT.check_adjacent_is_adjacent2vertex_inverse(adj, adj2v)
    adj2v = DT.Adjacent2Vertex(
        Dict(
            DT.BoundaryIndex => Set{NTuple{2,Int64}}([(1, 3), (3, 5), (5, 6), (10, 11), (6, 2), (2, 1)]),
            1 => Set{NTuple{2,Int64}}([(2, 3), (3, DT.BoundaryIndex), (DT.BoundaryIndex, 2)]),
            2 => Set{NTuple{2,Int64}}([(3, 1), (4, 3), (6, 4), (1, DT.BoundaryIndex), (DT.BoundaryIndex, 6)]),
            3 => Set{NTuple{2,Int64}}([(1, 2), (2, 4), (4, 5), (DT.BoundaryIndex, 1), (5, DT.BoundaryIndex)]),
            4 => Set{NTuple{2,Int64}}([(3, 2), (2, 6), (6, 5), (5, 3)]),
            5 => Set{NTuple{2,Int64}}([(3, 4), (4, 6), (6, DT.BoundaryIndex), (DT.BoundaryIndex, 3)]),
            6 => Set{NTuple{2,Int64}}([(4, 2), (5, 4), (DT.BoundaryIndex, 5), (2, DT.BoundaryIndex)])
        )
    )
    @test !DT.check_adjacent_is_adjacent2vertex_inverse(adj, adj2v)
    adj2v = DT.Adjacent2Vertex(
        Dict(
            DT.BoundaryIndex => Set{NTuple{2,Int64}}([(1, 3), (3, 5), (5, 6), (6, 2), (2, 1)]),
            1 => Set{NTuple{2,Int64}}([(2, 3), (3, DT.BoundaryIndex), (DT.BoundaryIndex, 2)]),
            2 => Set{NTuple{2,Int64}}([(3, 1), (4, 3), (6, 4), (1, DT.BoundaryIndex), (DT.BoundaryIndex, 6)]),
            3 => Set{NTuple{2,Int64}}([(1, 2), (2, 4), (4, 5), (DT.BoundaryIndex, 1), (5, DT.BoundaryIndex)]),
            4 => Set{NTuple{2,Int64}}([(3, 2), (2, 6), (6, 5), (5, 3)]),
            5 => Set{NTuple{2,Int64}}([(3, 4), (4, 6), (6, DT.BoundaryIndex), (DT.BoundaryIndex, 3)]),
            6 => Set{NTuple{2,Int64}}([(4, 2), (5, 4), (DT.BoundaryIndex, 5), (2, DT.BoundaryIndex)])
        )
    )
    DT.get_edge(adj, 92918, 29991) # Should just skip over the empty keys
    @test DT.check_adjacent_is_adjacent2vertex_inverse(adj, adj2v)
end

@testset "Can we correctly delete empty sets from Adjacent2Vertex?" begin
    adj2v = DT.Adjacent2Vertex(
        Dict(
            DT.BoundaryIndex => Set{NTuple{2,Int64}}([(1, 3), (3, 5), (5, 6), (6, 2), (2, 1)]),
            1 => Set{NTuple{2,Int64}}([(2, 3), (3, DT.BoundaryIndex), (DT.BoundaryIndex, 2)]),
            2 => Set{NTuple{2,Int64}}([(3, 1), (4, 3), (6, 4), (1, DT.BoundaryIndex), (DT.BoundaryIndex, 6)]),
            3 => Set{NTuple{2,Int64}}([(1, 2), (2, 4), (4, 5), (DT.BoundaryIndex, 1), (5, DT.BoundaryIndex)]),
            4 => Set{NTuple{2,Int64}}([(3, 2), (2, 6), (6, 5), (5, 3)]),
            5 => Set{NTuple{2,Int64}}([(3, 4), (4, 6), (6, DT.BoundaryIndex), (DT.BoundaryIndex, 3)]),
            6 => Set{NTuple{2,Int64}}([(4, 2), (5, 4), (DT.BoundaryIndex, 5), (2, DT.BoundaryIndex)]),
            7 => Set{NTuple{2,Int64}}(),
            10 => Set{NTuple{2,Int64}}()
        )
    )
    DT.clear_empty_keys!(adj2v)
    @test adj2v.adjacent2vertex == Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(1, 3), (3, 5), (5, 6), (6, 2), (2, 1)]),
        1 => Set{NTuple{2,Int64}}([(2, 3), (3, DT.BoundaryIndex), (DT.BoundaryIndex, 2)]),
        2 => Set{NTuple{2,Int64}}([(3, 1), (4, 3), (6, 4), (1, DT.BoundaryIndex), (DT.BoundaryIndex, 6)]),
        3 => Set{NTuple{2,Int64}}([(1, 2), (2, 4), (4, 5), (DT.BoundaryIndex, 1), (5, DT.BoundaryIndex)]),
        4 => Set{NTuple{2,Int64}}([(3, 2), (2, 6), (6, 5), (5, 3)]),
        5 => Set{NTuple{2,Int64}}([(3, 4), (4, 6), (6, DT.BoundaryIndex), (DT.BoundaryIndex, 3)]),
        6 => Set{NTuple{2,Int64}}([(4, 2), (5, 4), (DT.BoundaryIndex, 5), (2, DT.BoundaryIndex)])
    )
    DT.clear_empty_keys!(adj2v)
    @test adj2v.adjacent2vertex == Dict(
        DT.BoundaryIndex => Set{NTuple{2,Int64}}([(1, 3), (3, 5), (5, 6), (6, 2), (2, 1)]),
        1 => Set{NTuple{2,Int64}}([(2, 3), (3, DT.BoundaryIndex), (DT.BoundaryIndex, 2)]),
        2 => Set{NTuple{2,Int64}}([(3, 1), (4, 3), (6, 4), (1, DT.BoundaryIndex), (DT.BoundaryIndex, 6)]),
        3 => Set{NTuple{2,Int64}}([(1, 2), (2, 4), (4, 5), (DT.BoundaryIndex, 1), (5, DT.BoundaryIndex)]),
        4 => Set{NTuple{2,Int64}}([(3, 2), (2, 6), (6, 5), (5, 3)]),
        5 => Set{NTuple{2,Int64}}([(3, 4), (4, 6), (6, DT.BoundaryIndex), (DT.BoundaryIndex, 3)]),
        6 => Set{NTuple{2,Int64}}([(4, 2), (5, 4), (DT.BoundaryIndex, 5), (2, DT.BoundaryIndex)])
    )
end

@testset "Can we correctly clear degree 0 points from a DelaunayGraph?" begin
    DG = DT.DelaunayGraph(UndirectedGraph(
        [
            0 0 0
            0 1 1
            0 1 1
        ]
    ))
    DT.clear_empty_points!(DG)
    @test DG.graph == relabel(UndirectedGraph([1 1; 1 1]), Dict((1, 2) .=> (2, 3)))
    DT.clear_empty_points!(DG)
    @test DG.graph == relabel(UndirectedGraph([1 1; 1 1]), Dict((1, 2) .=> (2, 3)))
    DG = DT.DelaunayGraph(UndirectedGraph(
        [
            0 1 0 0
            1 1 1 0
            0 1 1 0
            0 0 0 0
        ]
    ))
    DT.clear_empty_points!(DG)
    @test DG.graph == relabel(UndirectedGraph([0 1 0; 1 1 1; 0 1 1]))
    DT.clear_empty_points!(DG)
    @test DG.graph == relabel(UndirectedGraph([0 1 0; 1 1 1; 0 1 1]))
    DG = DT.DelaunayGraph(UndirectedGraph(
        [
            0 1 0 0
            1 1 1 0
            0 1 1 0
            0 0 0 0
        ]
    ))
    DT.clear_empty_points!(DG)
    @test DG.graph == relabel(UndirectedGraph([0 1 0; 1 1 1; 0 1 1]))
    DT.clear_empty_points!(DG)
    @test DG.graph == relabel(UndirectedGraph([0 1 0; 1 1 1; 0 1 1]))
    DG = DT.DelaunayGraph(UndirectedGraph(
        [
            1 0 1
            0 0 0
            1 0 0
        ]
    ))
    DT.clear_empty_points!(DG)
    @test DG.graph == relabel(UndirectedGraph([1 1; 1 0]), Dict((1, 2) .=> (1, 3)))
end

@testset "Can we correctly identify a ghost edge?" begin
    @test DT.is_ghost_edge(2, 0)
    @test DT.is_ghost_edge(0, 2)
    @test DT.is_ghost_edge(0, 5)
    @test !DT.is_ghost_edge(5, 3)
    @test !DT.is_ghost_edge(3, 5)
    @test !DT.is_ghost_edge(7, 13)
end

@testset "Can we correctly identify a boundary point?" begin
    # Testing when ∂ ∈ DG.graph.V
    p1 = @SVector[-3.32, 3.53]
    p2 = @SVector[-5.98, 2.17]
    p3 = @SVector[-6.36, -1.55]
    p4 = @SVector[-2.26, -4.31]
    p5 = @SVector[6.34, -3.23]
    p6 = @SVector[-3.24, 1.01]
    p7 = @SVector[0.14, -1.51]
    p8 = @SVector[0.2, 1.25]
    p9 = @SVector[1.0, 4.0]
    p10 = @SVector[4.74, 2.21]
    p11 = @SVector[2.32, -0.27]
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
    T = Set{NTuple{3,Int64}}()
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DG = DT.DelaunayGraph{Int64}()
    for (i, j, k) in (
        (1, 2, 6),
        (1, 6, 8),
        (9, 1, 8),
        (9, 8, 10),
        (10, 8, 11),
        (8, 7, 11),
        (8, 6, 7),
        (6, 2, 3),
        (6, 3, 4),
        (6, 4, 7),
        (7, 4, 5),
        (11, 7, 5),
        (10, 11, 5)
    )
        DT.add_triangle!(i, j, k, T, adj, adj2v, DG; update_ghost_edges=true)
    end
    k = [9, 10, 5, 4, 3, 2, 1]
    @test all(DT.is_boundary_point(k, adj, DG) for k in k)
    k = [6, 7, 8, 11]
    @test all(!DT.is_boundary_point(k, adj, DG) for k in k)

    T, adj, adj2v, DG, HG = DT.triangulate_berg(pts)
    k = [9, 10, 5, 4, 3, 2, 1]
    @test all(DT.is_boundary_point(k, adj, DG) for k in k)
    k = [6, 7, 8, 11]
    @test all(!DT.is_boundary_point(k, adj, DG) for k in k)
end

@testset "Can we correctly identify if a given triangulation has ghost triangles?" begin
    p1 = @SVector[-3.32, 3.53]
    p2 = @SVector[-5.98, 2.17]
    p3 = @SVector[-6.36, -1.55]
    p4 = @SVector[-2.26, -4.31]
    p5 = @SVector[6.34, -3.23]
    p6 = @SVector[-3.24, 1.01]
    p7 = @SVector[0.14, -1.51]
    p8 = @SVector[0.2, 1.25]
    p9 = @SVector[1.0, 4.0]
    p10 = @SVector[4.74, 2.21]
    p11 = @SVector[2.32, -0.27]
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
    T = Set{NTuple{3,Int64}}()
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DG = DT.DelaunayGraph{Int64}()
    for (i, j, k) in (
        (1, 2, 6),
        (1, 6, 8),
        (9, 1, 8),
        (9, 8, 10),
        (10, 8, 11),
        (8, 7, 11),
        (8, 6, 7),
        (6, 2, 3),
        (6, 3, 4),
        (6, 4, 7),
        (7, 4, 5),
        (11, 7, 5),
        (10, 11, 5)
    )
        DT.add_triangle!(i, j, k, T, adj, adj2v, DG; update_ghost_edges=true)
    end
    @test DT.triangulation_has_ghost_triangles(adj, adj2v)
    p1 = (5.0, 6.0)
    p2 = (9.0, 6.0)
    p3 = (13.0, 5.0)
    p4 = (10.38, 0.0)
    p5 = (12.64, -1.69)
    p6 = (2.0, -2.0)
    p7 = (3.0, 4.0)
    p8 = (7.5, 3.53)
    p9 = (4.02, 1.85)
    p10 = (4.26, 0.0)
    pts = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]
    Random.seed!(928881)
    T, adj, adj2v, DG, HG = DT.triangulate_berg(pts)
    @test !DT.triangulation_has_ghost_triangles(adj, adj2v)
end

@testset "Can we correctly identify ghost triangles?" begin
    T = (3, 7, 5)
    @test !DT.is_ghost_triangle(T)
    T = (1, 2, 0)
    @test DT.is_ghost_triangle(T)
    @test DT.is_ghost_triangle(DT.shift_triangle_1(T))
    @test DT.is_ghost_triangle(DT.shift_triangle_2(T))
    @test DT.is_ghost_triangle(1, 2, 0)
    @test DT.is_ghost_triangle(0, 5, 3)
    @test DT.is_ghost_triangle(5, 0, 1)
    @test !DT.is_ghost_triangle(1, 3, 7)
end

@testset "Can we correctly rotate ghost triangles into the standard boundary configuration?" begin
    T = (3, 5, 0)
    newT = DT.rotate_ghost_triangle_to_boundary_form(T)
    @test newT == T
    T = (0, 4, 5)
    newT = DT.rotate_ghost_triangle_to_boundary_form(T)
    @test newT == (4, 5, 0)
    T = (1, 0, 7)
    newT = DT.rotate_ghost_triangle_to_boundary_form(T)
    @test newT == (7, 1, 0)
    i, j, k = 3, 5, 0
    u, v, w = DT.rotate_ghost_triangle_to_boundary_form(i, j, k)
    @test (u, v, w) == (3, 5, 0)
    i, j, k = 0, 4, 5
    u, v, w = DT.rotate_ghost_triangle_to_boundary_form(i, j, k)
    @test (u, v, w) == (4, 5, 0)
    i, j, k = 1, 0, 7
    u, v, w = DT.rotate_ghost_triangle_to_boundary_form(i, j, k)
    @test (u, v, w) == (7, 1, 0)
end

@testset "Can we correctly reset the centroid?" begin
    DT.CentroidCoordinates.x = 0.291
    DT.CentroidCoordinates.y = 0.4991
    DT.CentroidCoordinates.n = 99
    DT.clear_centroid_coordinates!()
    @test DT.CentroidCoordinates.x == 0.0
    @test DT.CentroidCoordinates.y == 0.0
    @test DT.CentroidCoordinates.n == 0
end

@testset "Can we correctly update the centroid incrementally?" begin
    pts = rand(SVector{2,Float64}, 9371)
    idx = shuffle(1:length(pts))
    DT.clear_centroid_coordinates!()
    for (j, i) in enumerate(idx)
        DT.update_centroid_after_new_point!(pts, i)
        pxy = sum(pts[idx[1:j]]) / j
        @test DT.CentroidCoordinates.x ≈ pxy[1]
        @test DT.CentroidCoordinates.y ≈ pxy[2]
        @test DT.CentroidCoordinates.n == j
    end
    idx = shuffle(1:(length(pts)-1)) # Don't want to go down to dividing by zero
    for (j, i) in enumerate(idx)
        DT.update_centroid_after_deleted_point!(pts, i)
        pxy = (sum(pts) - sum(pts[idx[1:j]])) / (length(pts) - j)
        @test DT.CentroidCoordinates.x ≈ pxy[1] rtol = 1e-7
        @test DT.CentroidCoordinates.y ≈ pxy[2] rtol = 1e-7
        @test DT.CentroidCoordinates.n == length(pts) - j
    end
end

@testset "Can we correctly access a set of points at the boundary index?" begin
    pts = rand(SVector{2,Float64}, 1827301)
    for i in shuffle(Int64[zeros(5000)..., (1:length(pts))...])
        if i == 0
            @test DT.get_point(pts, i) == (DT.CentroidCoordinates.x, DT.CentroidCoordinates.y)
        else
            @test DT.get_point(pts, i) == NTuple{2,Float64}(pts[i])
        end
    end
end

@testset "Can we correctly compute the centroid of some given points?" begin
    DT.CentroidCoordinates.x = 0.291
    DT.CentroidCoordinates.y = 0.4991
    DT.CentroidCoordinates.n = 99
    pts = rand(SVector{2,Float64}, 1827301)
    DT.compute_centroid!(pts)
    @test DT.CentroidCoordinates.x ≈ sum(pts)[1] / length(pts)
    @test DT.CentroidCoordinates.y ≈ sum(pts)[2] / length(pts)
    @test DT.CentroidCoordinates.n == length(pts)

    _pts = reduce(hcat, pts)
    DT.compute_centroid!(_pts)
    @test DT.CentroidCoordinates.x ≈ sum(pts)[1] / length(pts)
    @test DT.CentroidCoordinates.y ≈ sum(pts)[2] / length(pts)
    @test DT.CentroidCoordinates.n == size(_pts, 2)
end

@testset "Removing duplicate triangles" begin
    T = (1, 5, 3)
    @test DT.sort_triangle(T) == (1, 5, 3)
    T = (7, 1, 3)
    @test DT.sort_triangle(T) == (1, 3, 7)
    T = (2, 3, 1)
    @test DT.sort_triangle(T) == (1, 2, 3)

    T = Set{NTuple{3,Int64}}([
        (1, 2, 3),
        (10, 9, 3),
        (11, 5, 1),
        (193, 12, 10),
        (5, 3, 1),
        (19, 18, 17),
        (17, 5, 23),
        (20, 50, 72),
        (30, 31, 32),
        (20, 13, 37)
    ])
    V = DT.sort_triangles(T)
    @test V == Set{NTuple{3,Int64}}([
        (1, 2, 3),
        (3, 10, 9),
        (1, 11, 5),
        (10, 193, 12),
        (1, 5, 3),
        (17, 19, 18),
        (5, 23, 17),
        (20, 50, 72),
        (30, 31, 32),
        (13, 37, 20)
    ])
    @test DT.compare_triangle_sets(T, V)

    T = Set{NTuple{3,Int64}}([
        (1, 2, 3),
        (2, 3, 1),
        (3, 1, 2),
        (4, 5, 6),
        (11, 8, 3),
        (3, 11, 8),
        (18, 13, 1)
    ])
    V = DT.remove_duplicate_triangles(T)
    @test V == Set{NTuple{3,Int64}}([
        (1, 2, 3),
        (4, 5, 6),
        (1, 18, 13),
        (3, 11, 8),
    ])
    T = Set{NTuple{3,Int64}}([
        (11, 9, 1),
        (1, 11, 9)
    ])
    V = DT.remove_duplicate_triangles(T)
    @test V == Set{NTuple{3,Int64}}([
        (1, 11, 9)
    ])
end

@testset "Testing if a triangle is a boundary triangle" begin
    a = 0.0
    b = 10.0
    c = 0.0
    d = 10.0
    nx = 11
    ny = 11
    function DT._get_point(pts::AbstractMatrix, i)
        return @view pts[:, i]
    end
    function DT._eachindex(pts::AbstractMatrix)
        return axes(pts, 2)
    end
    T, adj, adj2v, DG, pts = triangulate_structured(a, b, c, d, nx, ny)
    @test DT.is_boundary_triangle((1, 2, 12), adj)
    @test !DT.is_boundary_triangle((12, 2, 1), adj)
    @test DT.is_boundary_triangle((2, 3, 13), adj)
end

@testset "Setdiff for triangles" begin
    V = Set{NTuple{3,Int64}}([
        (1, 2, 3),
        (4, 5, 6),
        (7, 8, 9)
    ])
    _V = DT.setdiff_triangles(V, V)
    @test _V == Set{NTuple{3,Int64}}()
    T = Set{NTuple{3,Int64}}([
        (1, 2, 3),
        (13, 11, 3),
        (7, 10, 11),
        (5, 6, 23),
        (5, 10, 173),
    ])
    _T = Set{NTuple{3,Int64}}([
        (11, 3, 13),
        (5, 10, 200),
        (23, 5, 6),
        (7, 10, 11)
    ])
    @test DT.setdiff_triangles(T, _T) == Set{NTuple{3,Int64}}([
        (1, 2, 3),
        (5, 10, 173)
    ])
end

@testset "Sorting a boundary" begin
    Random.seed!(228818181)
    A = (0.4, 2.6)#1
    B = (1.1428403421814, 2.9209248186714)#2
    C = (2.0927803855902, 2.381117619401)#3
    D = (2.011356953298, 1.6271969500289)#4
    E = (1.4172674658328, 1.0300917798863)#5
    F = (0.4160608169067, 1.4402246240247)#6
    G = (0.0571945782857, 2.1700198319768)#7
    ref = (1.1790285343113, 2.121768909137)
    pts = [A, B, C, D, E, F, G]
    θ = zeros(length(pts))
    idx_vec = zeros(Int64, length(pts))
    boundary_idx = shuffle(collect(1:7))
    DelaunayTriangulation.sort_boundary!(θ, idx_vec, pts, boundary_idx, ref)
    @test boundary_idx == [6, 5, 4, 3, 2, 1, 7]

    for _ in 1:250
        θ = LinRange(0, 2π - 1 / 200, 250)
        x = cos.(θ)
        y = sin.(θ)
        pts = [[x, y] for (x, y) in zip(x, y)]
        ref = [0.0, 0.0]
        ψ = zeros(length(θ))
        idx_vec = zeros(Int64, length(pts))
        boundary_idx = collect(1:length(pts))
        shuffle!(boundary_idx)
        DelaunayTriangulation.sort_boundary!(ψ, idx_vec, pts, boundary_idx, ref)
        shift_idx = findfirst(boundary_idx .== 1)
        circshift!(boundary_idx, shift_idx - 1)
        @test boundary_idx == collect(1:length(pts))
    end

    for _ in 1:250
        θ = LinRange(0, 2π - 1 / 200, 250)
        x = cos.(θ)
        y = sin.(θ)
        pts = [[x, y] for (x, y) in zip(x, y)]
        ref = [0.0, 0.0]
        boundary_idx = collect(1:length(pts))
        shuffle!(boundary_idx)
        DelaunayTriangulation.sort_boundary!(pts, boundary_idx, ref)
        shift_idx = findfirst(boundary_idx .== 1)
        circshift!(boundary_idx, shift_idx - 1)
        @test boundary_idx == collect(1:length(pts))
    end

    for _ in 1:250
        pts = rand(SVector{2,Float64}, 500)
        T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
        ref = get_point(pts, DT.BoundaryIndex)
        BN = convex_hull(DG, pts)
        DG_BN = DT.get_neighbour(DG, DT.BoundaryIndex) |> collect
        θ = zeros(length(DG_BN))
        idx_vec = zeros(Int64, length(DG_BN))
        DT.sort_boundary!(θ, idx_vec, pts, DG_BN, ref)
        shift_idx = findfirst(DG_BN .== BN[1])
        circshift!(DG_BN, shift_idx - 1)
        @test DG_BN == BN
    end

    for _ in 1:2500
        pts = rand(SVector{2,Float64}, 500)
        T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
        ref = get_point(pts, DT.BoundaryIndex)
        BN = convex_hull(DG, pts)
        DG_BN = DT.get_neighbour(DG, DT.BoundaryIndex) |> collect
        DT.sort_boundary!(pts, DG_BN, ref)
        shift_idx = findfirst(DG_BN .== BN[1])
        circshift!(DG_BN, shift_idx - 1)
        @test DG_BN == BN
    end
end

@testset "Sorting a boundary for a point on the boundary" begin
    Random.seed!(2992991)
    qa = (-1.4, 0.61)
    qb = (1.0, 1.0)
    qc = (0.0, -2.0)
    qd = (1.0, 3.0)
    qe = (0.0, 0.0)
    qf = (-0.56, 1.53)
    qg = (2.22, 1.43)
    qh = (2.46, -0.57)
    qi = (-0.68, -1.07)
    pts = [qa, qb, qc, qd, qe, qf, qg, qh, qi]
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts; trim=false)
    i = 6
    BN = [5, 2, 4, 0, 1]
    DT.sort_boundary!(pts, BN, pts[i])
    @test BN == [1, 5, 2, 4, 0]
    i = 1
    BN = [6, 9, 5, 0]
    DT.sort_boundary!(pts, BN, pts[i])
    @test BN == [9, 5, 6, 0]
    i = 9
    BN = DT.get_neighbour(DG, i) |> collect
    DT.sort_boundary!(pts, BN, pts[i])
    @test BN == [0, 3, 5, 1]
    i = 4
    BN = DT.get_neighbour(DG, i) |> collect
    DT.sort_boundary!(pts, BN, pts[i])
    @test BN == [6, 2, 7, 0]
end

@testset "Area of a polygon" begin
    Random.seed!(2992991)
    qa = (-1.4, 0.61)
    qb = (1.0, 1.0)
    qc = (0.0, -2.0)
    qd = (1.0, 3.0)
    qe = (0.0, 0.0)
    qf = (-0.56, 1.53)
    qg = (2.22, 1.43)
    qh = (2.46, -0.57)
    qi = (-0.68, -1.07)
    pts = [qa, qb, qc, qd, qe, qf, qg, qh, qi]
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
    BN = convex_hull(DG, pts)
    pts_BN = pts[BN]
    @test DT.area(pts_BN) ≈ 11.6082
    reverse!(pts_BN)
    @test DT.area(pts_BN) ≈ -11.6082
end

@testset "Finding the first boundary index" begin
    v = [1, 3, 0, 0]
    @test DT.find_first_boundary_index(v) == 3
    v = [0, 1, 3, 0]
    @test DT.find_first_boundary_index(v) == 4
    v = [0, 0, 8, 10]
    @test DT.find_first_boundary_index(v) == 1
    v = [9, 1, 3, 0, 0, 5]
    @test DT.find_first_boundary_index(v) == 4
    v = [1, 3, 5, 7]
    @test DT.find_first_boundary_index(v) === nothing
end

@testset "Two-dimensional cross product" begin
    u = [2.0, 3.0]
    v = [6.3, 7.1]
    @test DT.cross_2d(u, v) == cross([u..., 0.0], [v..., 0.0])[3]
    for _ in 1:500
        u = 15rand(2)
        v = 15rand(2)
        @test DT.cross_2d(u, v) == cross([u..., 0.0], [v..., 0.0])[3]
    end
end

@testset "Intersections" begin
    p = (-8.04, 4.74)
    q = (-1.56, -5.32)
    r = (3.64, 4.7)
    s = (-12.78, -4.28)
    @test DT.intersection_of_two_line_segments(p, q, r, s) |> collect ≈
          DT.intersection_of_two_line_segments(r, s, p, q) |> collect ≈
          DT.intersection_of_two_line_segments(q, p, r, s) |> collect ≈
          DT.intersection_of_two_line_segments(q, p, s, r) |> collect ≈
          [-4.9782513757098, -0.0132702407962]
    p = (-3.6776328459745, 2.7776569098864)
    q = (-6.3950281994953, 2.5340463209513)
    r = (-2.70866118888, 1.9910799176567)
    s = (-5.7635834535008, 2.977401142607)
    @test DT.intersection_of_two_line_segments(p, q, r, s) |> collect ≈
          DT.intersection_of_two_line_segments(r, s, p, q) |> collect ≈
          DT.intersection_of_two_line_segments(q, p, r, s) |> collect ≈
          DT.intersection_of_two_line_segments(q, p, s, r) |> collect ≈
          [-4.8260418530713, 2.6747036924995]
    p = (3.0, 6.0)
    q = (3.0, 0.0)
    r = (5.0, 2.0)
    s = (1.0, 2.0)
    @test DT.intersection_of_two_line_segments(p, q, r, s) |> collect ≈
          DT.intersection_of_two_line_segments(r, s, p, q) |> collect ≈
          DT.intersection_of_two_line_segments(q, p, r, s) |> collect ≈
          DT.intersection_of_two_line_segments(q, p, s, r) |> collect ≈
          [3.0, 2.0]
end

@testset "Circular indices" begin
    v = [1, 2, 3, 4, 5]
    @test DT.nextindex_circular(v, 1) == 2
    @test DT.nextindex_circular(v, 2) == 3
    @test DT.nextindex_circular(v, 3) == 4
    @test DT.nextindex_circular(v, 4) == 5
    @test DT.nextindex_circular(v, 5) == 1
    @test DT.previndex_circular(v, 1) == 5
    @test DT.previndex_circular(v, 2) == 1
    @test DT.previndex_circular(v, 3) == 2
    @test DT.previndex_circular(v, 4) == 3
    @test DT.previndex_circular(v, 5) == 4
    for _ in 1:209
        v = rand(Float64, rand(1:3000))
        for i in eachindex(v)
            ni = i == length(v) ? 1 : i + 1
            pri = i == 1 ? length(v) : i - 1
            @test DT.nextindex_circular(v, i) == ni
            @test DT.previndex_circular(v, i) == pri
        end
    end
end

@testset "Testing if an index is a vertex of triangle" begin
    T = (1, 5, 7)
    @test DT.is_vertex_of(T, 1)
    @test DT.is_vertex_of(T, 5)
    @test DT.is_vertex_of(T, 7)
    @test !DT.is_vertex_of(T, 3)
end

@testset "Comparing individual triangles" begin
    T = (1, 2, 3)
    V = (2, 3, 1)
    @test DelaunayTriangulation.compare_triangles(T, V)
    @test DelaunayTriangulation.compare_triangles(V, T)
    T = (2, 3, 1)
    @test DelaunayTriangulation.compare_triangles(T, V)
    @test DelaunayTriangulation.compare_triangles(V, T)
    T = (3, 1, 2)
    @test DelaunayTriangulation.compare_triangles(T, V)
    @test DelaunayTriangulation.compare_triangles(V, T)
    T = (1, 3, 2)
    @test !DelaunayTriangulation.compare_triangles(T, V)
    @test !DelaunayTriangulation.compare_triangles(V, T)
    T = (5, 7, 1)
    @test !DelaunayTriangulation.compare_triangles(T, V)
    @test !DelaunayTriangulation.compare_triangles(V, T)
end

@testset "Testing if a point set is unique" begin
    pts = [[1.0, 2.0], [3.0, 4.0]]
    @test DT.all_points_are_unique(pts)
    push!(pts, [3.0, 4.0])
    @test !DT.all_points_are_unique(pts)
    pts = [1.0 2.0 3.0 1.37; 3.5 1.2 0.0 7.1]
    @test DT.all_points_are_unique(pts)
    pts = [1.0 2.3 4.5 1.0; 3.5 2.3 4.5 3.5]
    @test !DT.all_points_are_unique(pts)
end

@testset "Midpoint computation" begin
    pts = rand(2, 500)
    T, adj, adj2v, DG = triangulate_bowyer(pts)
    for e in edges(DG)
        u, v = e
        pu, pv = get_point(pts, u, v)
        mid = (pu .+ pv) ./ 2
        midx, midy = DelaunayTriangulation.edge_midpoint(e, pts)
        @test midx ≈ mid[1]
        @test midy ≈ mid[2]
    end
end

@testset "Length computation" begin
    pts = rand(2, 500)
    T, adj, adj2v, DG = triangulate_bowyer(pts)
    for e in edges(DG)
        u, v = e
        pu, pv = get_point(pts, u, v)
        ℓ = norm(pv .- pu)
        ℓs = DelaunayTriangulation.edge_length(e, pts)
        @test ℓ ≈ ℓs
    end
end

@testset "Testing if opposite angle is obtuse" begin
    for _ in 1:500
        pts = rand(2, 500)
        T, adj, adj2v, DG = triangulate_bowyer(pts)
        for V in T
            i, j, k = indices(V)
            pᵢ, pⱼ, pₖ = get_point(pts, i, j, k)
            e1 = pⱼ .- pᵢ
            e2 = pₖ .- pⱼ
            e3 = pᵢ .- pₖ
            θ1 = acos(dot(e1, e2) / (norm(e1) * norm(e2)))
            θ2 = acos(dot(e2, e3) / (norm(e2) * norm(e3)))
            θ3 = acos(dot(e3, e1) / (norm(e3) * norm(e1)))
            angle_opposite_1 = θ2
            angle_opposite_2 = θ3
            angle_opposite_3 = θ1
            angle_opposite_1_obtuse = θ2 > π / 2
            angle_opposite_2_obtuse = θ3 > π / 2
            angle_opposite_3_obtuse = θ1 > π / 2
            @test angle_opposite_1_obtuse == DT.opposite_angle_is_obtuse(norm(e2), norm(e3), norm(e1))
            @test angle_opposite_2_obtuse == DT.opposite_angle_is_obtuse(norm(e3), norm(e1), norm(e2))
            @test angle_opposite_3_obtuse == DT.opposite_angle_is_obtuse(norm(e1), norm(e2), norm(e3))
        end
    end
end