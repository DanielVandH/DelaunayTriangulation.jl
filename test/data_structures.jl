###################################################
#/
#/
#/ Constructors for the triangulation data structures
#/
#/
###################################################
@testset "Adjacent" begin
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    @test adj.adjacent == DefaultDict{NTuple{2,Int64},Int64}(DT.DefaultAdjacentValue)
    @test typeof(adj) == DT.Adjacent{Int64,NTuple{2,Int64}}
    @test isconcretetype(typeof(adj))
    adj = DT.Adjacent{Int16,NTuple{2,Int16}}()
    @test adj.adjacent == DefaultDict{NTuple{2,Int16},Int16}(Int16(DT.DefaultAdjacentValue))
    @test typeof(adj) == DT.Adjacent{Int16,NTuple{2,Int16}}
    @test isconcretetype(typeof(adj))
    @test DT.adjacent(adj) == DefaultDict{NTuple{2,Int16},Int16}(Int16(DT.DefaultAdjacentValue))

    adj = DT.Adjacent(DefaultDict(DT.DefaultAdjacentValue,
        Dict(
            (6, 3) => 1,
            (3, 1) => 6,
            (1, 6) => 3,
            (6, 2) => 3,
            (2, 3) => 6,
            (3, 6) => 2,
            (3, 2) => 5,
            (2, 5) => 3,
            (5, 3) => 2,
            (4, 1) => 5,
            (1, 5) => 4,
            (5, 4) => 1,
            (4, 6) => 1,
            (6, 1) => 4,
            (1, 4) => 6,
            (5, 1) => 3,
            (1, 3) => 5,
            (3, 5) => 1,
            (4, 5) => DT.BoundaryIndex,
            (5, 2) => DT.BoundaryIndex,
            (2, 6) => DT.BoundaryIndex,
            (6, 4) => DT.BoundaryIndex
        )
    ))
    @test DT.get_edge(adj, 6, 3) == 1
    @test DT.get_edge(adj, 3, 6) == 2
    @test DT.get_edge(adj, 6, 10) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 6, 4) == DT.BoundaryIndex
end

@testset "Adjacent2Vertex" begin
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    @test adj2v.adjacent2vertex == Dict{NTuple{2,16},Int64}()
    @test typeof(adj2v) == DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}
    @test isconcretetype(typeof(adj2v))
    adj2v = DT.Adjacent2Vertex{Int16,Set{NTuple{2,Int16}},NTuple{2,Int16}}()
    @test adj2v.adjacent2vertex == Dict{NTuple{2,Int16},Int16}()
    @test typeof(adj2v) == DT.Adjacent2Vertex{Int16,Set{NTuple{2,Int16}},NTuple{2,Int16}}
    @test isconcretetype(typeof(adj2v))
    @test DT.adjacent2vertex(adj2v) == Dict{Int64,Set{NTuple{2,Int16}}}()

    adj2v = DT.Adjacent2Vertex(Dict(
        DT.BoundaryIndex => [(4, 5), (5, 2), (2, 6), (6, 4)],
        1 => [(5, 4), (3, 5), (6, 3), (4, 6)],
        2 => [(5, 3), (3, 6)],
        3 => [(1, 6), (5, 1), (2, 5), (6, 2)],
        4 => [(1, 5), (6, 1)],
        5 => [(4, 1), (1, 3), (3, 2)],
        6 => [(1, 4), (3, 1), (2, 3)]
    ))
    @test adj2v isa DT.Adjacent2Vertex{Int64,Vector{NTuple{2,Int64}},NTuple{2,Int64}}
    @test DT.get_edge(adj2v, 1) == [(5, 4), (3, 5), (6, 3), (4, 6)]
    @test DT.get_edge(adj2v, 2) == [(5, 3), (3, 6)]
    @test DT.get_edge(adj2v, 3) == [(1, 6), (5, 1), (2, 5), (6, 2)]
    @test DT.get_edge(adj2v, 4) == [(1, 5), (6, 1)]
    @test DT.get_edge(adj2v, 5) == [(4, 1), (1, 3), (3, 2)]
    @test DT.get_edge(adj2v, 6) == [(1, 4), (3, 1), (2, 3)]
    @test DT.get_edge(adj2v, DT.BoundaryIndex) == [(4, 5), (5, 2), (2, 6), (6, 4)]
end

@testset "DelaunayGraph" begin
    dg = DT.DelaunayGraph()
    @test dg.graph == UndirectedGraph{Int64}()
    @test typeof(dg) == DT.DelaunayGraph{Int64}
    @test isconcretetype(typeof(dg))
    dg = DT.DelaunayGraph{Int16}()
    @test dg.graph == UndirectedGraph{Int16}()
    @test typeof(dg) == DT.DelaunayGraph{Int16}
    @test isconcretetype(typeof(dg))
    @test DT.graph(dg) == UndirectedGraph{Int16}()
    tvg = UndirectedGraph{Int64}()
    add!(tvg, 7)
    add!(tvg, 9)
    add!(tvg, 9, 7)
    dg = DT.DelaunayGraph(tvg)
    @test DT.graph(dg) == tvg
end

@testset "HistoryGraph" begin
    hg = DT.HistoryGraph{NTuple{3,Int64}}()
    @test hg.graph == DirectedGraph{NTuple{3,Int64}}()
    @test typeof(hg) == DT.HistoryGraph{NTuple{3,Int64}}
    @test isconcretetype(typeof(hg))
    hg = DT.HistoryGraph{NTuple{3,Int64}}()
    @test hg.graph == DirectedGraph{NTuple{3,Int64}}()
    @test typeof(hg) == DT.HistoryGraph{NTuple{3,Int64}}
    @test DT.graph(hg) == DirectedGraph{NTuple{3,Int64}}()
    dg = DirectedGraph{NTuple{3,Int64}}()
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9))
end

###################################################
#/
#/
#/ Indexing and updating of the triangulation data structures
#/
#/
###################################################
@testset "Adjacent" begin
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    i, j, k, r = 1, 2, 3, 4
    DT.add_edge!(adj, i, j, r)
    DT.add_edge!(adj, j, r, i)
    DT.add_edge!(adj, r, i, j)
    DT.add_edge!(adj, i, k, j)
    DT.add_edge!(adj, k, j, i)
    DT.add_edge!(adj, j, i, k)
    @test DT.get_edge(adj, i, j) == r
    @test DT.get_edge(adj, j, r) == i
    @test DT.get_edge(adj, r, i) == j
    @test DT.get_edge(adj, i, k) == j
    @test DT.get_edge(adj, k, j) == i
    @test DT.get_edge(adj, j, i) == k
    @test DT.get_edge(adj, (j, r)) == i
    DT.add_edge!(adj, (5, 7), 9)
    DT.add_edge!(adj, (7, 5), 3)
    @test DT.get_edge(adj, 5, 7) == DT.get_edge(adj, (5, 7)) == 9
    DT.delete_edge!(adj, 5, 7)
    DT.delete_edge!(adj, 7, 5)
    @test DT.get_edge(adj, 5, 7) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, (5, 7)) == DT.DefaultAdjacentValue
    DT.delete_edge!(adj, i, j)
    DT.delete_edge!(adj, j, i)
    @test DT.get_edge(adj, i, j) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, j, i) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, (i, j)) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, (j, i)) == DT.DefaultAdjacentValue

    DT.add_edge!(adj, 16, 13, 5)
    DT.add_edge!(adj, 13, 16, DT.BoundaryIndex)
    DT.delete_edge!(adj, 16, 13)
    @test DT.get_edge(adj, (16, 13)) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, 13, 16) == DT.BoundaryIndex
    DT.add_edge!(adj, 16, 13, 5)
    DT.delete_edge!(adj, 16, 13)
    DT.delete_edge!(adj, 13, 16)
    @test DT.get_edge(adj, (16, 13)) == DT.DefaultAdjacentValue
    @test DT.get_edge(adj, (13, 16)) == DT.DefaultAdjacentValue
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    T1 = DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3)
    T2 = DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6)
    T3 = DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9)
    DT.add_triangle!(adj, T1)
    @test DT.get_edge(adj, 1, 2) == 3
    @test DT.get_edge(adj, 2, 3) == 1
    @test DT.get_edge(adj, 3, 1) == 2
    DT.add_triangle!(adj, T2, T3)
    @test DT.get_edge(adj, 4, 5) == 6
    @test DT.get_edge(adj, 5, 6) == 4
    @test DT.get_edge(adj, 6, 4) == 5
    @test DT.get_edge(adj, 7, 8) == 9
    @test DT.get_edge(adj, 8, 9) == 7
    @test DT.get_edge(adj, 9, 7) == 8
    @test length(DT.adjacent(adj)) == 9

    p1 = (-1.0, 4.0)
    p2 = (4.0, 6.0)
    p3 = (2.0, 2.0)
    p4 = (6.0, -3.0)
    p5 = (7.0, 3.0)
    pts = [p1, p2, p3, p4, p5]
    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 1, 3, 2))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 1, 4, 3))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 3, 4, 5))
    DT.add_triangle!(adj, DT.construct_triangle(NTuple{3,Int64}, 3, 5, 2))
    @test DT.get_edge(adj, 1, 3) == 2
    @test DT.get_edge(adj, 3, 2) == 1
    @test DT.get_edge(adj, 2, 1) == 3
    @test DT.get_edge(adj, 1, 4) == 3
    @test DT.get_edge(adj, 4, 3) == 1
    @test DT.get_edge(adj, 3, 1) == 4
    @test DT.get_edge(adj, 3, 4) == 5
    @test DT.get_edge(adj, 4, 5) == 3
    @test DT.get_edge(adj, 5, 3) == 4
    @test DT.get_edge(adj, 3, 5) == 2
    @test DT.get_edge(adj, 5, 2) == 3
    @test DT.get_edge(adj, 2, 3) == 5
    DT.delete_edge!(adj, 2, 3)
    @test DT.get_edge(adj, 3, 2) == 1
    @test DT.get_edge(adj, (2, 3)) == DT.DefaultAdjacentValue

    adj = DT.Adjacent{Int64,NTuple{2,Int64}}()
    DT.add_edge!(adj, 2, 3, 7)
    DT.add_edge!(adj, 5, 7, 17)
    @test collect(DT.edges(adj)) == [(2, 3), (5, 7)] || collect(DT.edges(adj)) == [(5, 7), (2, 3)]

    @inferred DT.get_edge(adj, 2, 3)
    @inferred DT.get_edge(adj, (2, 3))
    @inferred DT.get_edge(adj, -3, 50)
end

@testset "Adjacent2Vertex" begin
    adj2v = DT.Adjacent2Vertex{Int64,Set{NTuple{2,Int64}},NTuple{2,Int64}}()
    DT.add_edge!(adj2v, 3, (1, 2))
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([(1, 2)])
    DT.add_edge!(adj2v, 10, 5, 7)
    @test DT.get_edge(adj2v, 10) == Set{NTuple{2,Int64}}([(5, 7)])
    DT.add_edge!(adj2v, 3, (1, 4))
    DT.add_edge!(adj2v, 3, (10, 11))
    DT.add_edge!(adj2v, 3, (13, 5))
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([(1, 2), (1, 4), (10, 11), (13, 5)])
    DT.delete_edge!(adj2v, 3, (1, 4))
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([(1, 2), (10, 11), (13, 5)])
    DT.delete_edge!(adj2v, 3, 1, 2)
    @test DT.get_edge(adj2v, 3) == Set{NTuple{2,Int64}}([(10, 11), (13, 5)])
    DT.delete_point!(adj2v, 10)
    @test_throws KeyError DT.get_edge(adj2v, 10)
end

@testset "DelaunayGraph" begin
    DG = DT.DelaunayGraph{Int64}()
    DT.add_point!(DG, 5)
    @test 5 ??? DT.graph(DG).V
    DT.add_point!(DG, 17, 53, 32)
    @test all(v ??? DT.graph(DG).V for v in [5, 17, 53, 32])
    DT.add_neighbour!(DG, 5, 17)
    @test DT.graph(DG).N[5] == Set(17)
    DT.add_neighbour!(DG, 5, 53, 32)
    @test DT.graph(DG).N[5] == Set([17, 53, 32])
    @test DT.get_neighbour(DG, 5) == Set([17, 53, 32])
    @test DT.get_neighbour(DG, 53) == Set(5)
    @test DT.get_neighbour(DG, 32) == Set(5)
    DT.delete_neighbour!(DG, 32, 5)
    @test DT.get_neighbour(DG, 5) == Set([17, 53])
    @test DT.get_neighbour(DG, 32) == Set()
    DT.delete_point!(DG, 5)
    @test_throws KeyError DT.get_neighbour(DG, 5)
    DT.add_point!(DG, 17, 101, 291)
    @test 17 ??? DT.graph(DG).V
    @test 101 ??? DT.graph(DG).V
    @test 291 ??? DT.graph(DG).V
    DT.delete_point!(DG, 17, 101, 291)
    @test 17 ??? DT.graph(DG).V
    @test 101 ??? DT.graph(DG).V
    @test 291 ??? DT.graph(DG).V
    @test edges(DG) == DG.graph.E

    tvn = UndirectedGraph{Int64}()
    add!(tvn, 1)
    add!(tvn, 2)
    add!(tvn, 3)
    [add!(tvn, 1, i) for i in [4, 5, 6, 9, 3]]
    [add!(tvn, 2, i) for i in [11, 9, 8, 15, 16]]
    add!(tvn, 3, 1)
    DG = DT.DelaunayGraph(tvn)
    @test DT.graph(DG) == tvn
    @test DT.get_neighbour(DG, 1) == Set([4, 5, 6, 9, 3])
    @test DT.get_neighbour(DG, 2) == Set([11, 9, 8, 15, 16])
    @test DT.get_neighbour(DG, 3) == Set([1])
    DT.add_point!(DG, 40)
    DT.add_neighbour!(DG, 40, 50)
    @test DT.get_neighbour(DG, 40) == Set([50])
    DT.add_neighbour!(DG, 60, 70)
    @test DT.get_neighbour(DG, 60) == Set([70])
    @test DT.get_neighbour(DG, 70) == Set([60])
    DT.add_point!(DG, 70, -1, -2, -3, -4, -5, -6, -7)
    @test all(has(DT.graph(DG), u) for u ??? [70, -1, -2, -3, -4, -5, -6, -7])
    DT.add_neighbour!(DG, 70, 4, 5, 9, 3)
    @test DT.get_neighbour(DG, 70) == Set([4, 60, 5, 9, 3])
    DT.add_neighbour!(DG, 1, 10)
    @test DT.get_neighbour(DG, 1) == Set([4, 5, 6, 9, 3, 10])
    DT.add_neighbour!(DG, 1, 10, 11, 12, 13)
    @test DT.get_neighbour(DG, 1) == Set([4, 5, 6, 9, 3, 10, 11, 12, 13])
    DG_empty = DT.DelaunayGraph()
    @test isempty(DG_empty.graph.V)
    DT.add_neighbour!.(Ref(DG_empty), 2, [3, 4, 5, 1])
    @test DT.get_neighbour(DG_empty, 2) == Set([3, 4, 5, 1])
    @test edges(DG) == DG.graph.E

    tvn = UndirectedGraph{Int64}()
    add!(tvn, 1)
    add!(tvn, 2)
    add!(tvn, 3)
    [add!(tvn, 1, i) for i in [4, 5, 6, 9, 3]]
    [add!(tvn, 2, i) for i in [11, 9, 8, 15, 16]]
    add!(tvn, 3, 1)
    DG = DT.DelaunayGraph(tvn)
    DT.delete_neighbour!(DG, 1, 6)
    @test DT.get_neighbour(DG, 1) == Set([3, 4, 5, 9])
    @test DT.get_neighbour(DG, 2) == Set([11, 9, 8, 15, 16])
    @test DT.get_neighbour(DG, 3) == Set([1])
    @test edges(DG) == DG.graph.E

    pts = rand(SVector{2,Float64}, 250)
    T, adj, adj2v, DG = DT.triangulate_bowyer(pts)
    for i in eachindex(pts)
        @test deg(graph(DG), i) == DT.num_neighbours(DG, i)
    end
end

@testset "HistoryGraph" begin
    dg = DirectedGraph{NTuple{3,Int64}}()
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6))
    add!(dg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9))
    hg = DT.HistoryGraph(dg)
    @test out_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3)) == Set([DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6), DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9)])
    @test in_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9)) == [DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3)]
    @test in_deg(hg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3)) == 0
    @test in_deg(hg, DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9)) == 1
    @test out_deg(hg, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3)) == 2
    @test out_deg(hg, DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9)) == 0
    DT.add_triangle!(hg, DT.construct_triangle(NTuple{3,Int64}, 11, 13, 15))
    @test DT.construct_triangle(NTuple{3,Int64}, 11, 13, 15) ??? hg.graph.V
    DT.add_triangle!(hg, DT.construct_triangle(NTuple{3,Int64}, 17, 18, 19), DT.construct_triangle(NTuple{3,Int64}, 22, 25, 29))
    @test DT.construct_triangle(NTuple{3,Int64}, 17, 18, 19) ??? hg.graph.V
    @test DT.construct_triangle(NTuple{3,Int64}, 22, 25, 29) ??? hg.graph.V
    DT.add_triangle!(hg, DT.construct_triangle(NTuple{3,Int64}, 18, 19, 17))
    DT.add_edge!(hg, DT.construct_triangle(NTuple{3,Int64}, 17, 18, 19), DT.construct_triangle(NTuple{3,Int64}, 22, 25, 29))
    @test DT.construct_triangle(NTuple{3,Int64}, 22, 25, 29) ??? out_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, 17, 18, 19))
    DT.add_edge!(hg, DT.construct_triangle(NTuple{3,Int64}, 18, 19, 17), DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3))
    @test DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3) ??? out_neighbors(hg, DT.construct_triangle(NTuple{3,Int64}, 17, 18, 19))
    DT.add_edge!(hg, DT.construct_triangle(NTuple{3,Int64}, 9, 7, 8), DT.construct_triangle(NTuple{3,Int64}, 13, 15, 11))
    og = out_deg(hg, DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9))
    DT.add_edge!(hg, DT.construct_triangle(NTuple{3,Int64}, 9, 7, 8), DT.construct_triangle(NTuple{3,Int64}, 13, 15, 11))
    DT.add_edge!(hg, DT.construct_triangle(NTuple{3,Int64}, 9, 7, 8), DT.construct_triangle(NTuple{3,Int64}, 13, 15, 11))
    DT.add_edge!(hg, DT.construct_triangle(NTuple{3,Int64}, 9, 7, 8), DT.construct_triangle(NTuple{3,Int64}, 13, 15, 11))
    DT.add_edge!(hg, DT.construct_triangle(NTuple{3,Int64}, 9, 7, 8), DT.construct_triangle(NTuple{3,Int64}, 13, 15, 11))
    @test og == out_deg(hg, DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9))

    H = DT.HistoryGraph{NTuple{3,Int64}}()
    DT.add_triangle!(H, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3))
    DT.add_triangle!(H, DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6))
    DT.add_edge!(H, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6))
    Htrue = deepcopy(H)
    DT.add_edge!(H, DT.construct_triangle(NTuple{3,Int64}, 2, 3, 1), DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6))
    @test Htrue.graph == H.graph
    DT.add_edge!(H, DT.construct_triangle(NTuple{3,Int64}, 2, 3, 1), DT.construct_triangle(NTuple{3,Int64}, 5, 6, 4))
    @test Htrue.graph == H.graph
    DT.add_edge!(H, DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3), DT.construct_triangle(NTuple{3,Int64}, 6, 4, 5))
    @test Htrue.graph == H.graph
    DT.add_triangle!(H, DT.construct_triangle(NTuple{3,Int64}, 3, 1, 2))
    @test Htrue.graph == H.graph

    H = DT.HistoryGraph{NTuple{3,Int64}}()
    HH = DT.HistoryGraph{NTuple{3,Int64}}()
    T??? = DT.construct_triangle(NTuple{3,Int64}, 1, 2, 3)
    T??? = DT.construct_triangle(NTuple{3,Int64}, 4, 5, 6)
    T??? = DT.construct_triangle(NTuple{3,Int64}, 7, 8, 9)
    DT.add_triangle!(H, T???, T???, T???)
    DT.add_triangle!(HH, T???)
    DT.add_triangle!(HH, T???)
    DT.add_triangle!(HH, T???)
    @test DT.graph(H) == DT.graph(HH)
    DT.add_triangle!(H, DT.construct_triangle(NTuple{3,Int64}, 2, 3, 1)) # same triangle as T???
    @test DT.graph(H) == DT.graph(HH)
    DT.add_triangle!(H, DT.construct_triangle(NTuple{3,Int64}, 2, 3, 1), DT.construct_triangle(NTuple{3,Int64}, 9, 7, 8)) # same as T??? and T???
    @test DT.graph(H) == DT.graph(HH)
    DT.add_edge!(H, T???, T???, T???)
    DT.add_edge!(HH, T???, T???)
    DT.add_edge!(HH, T???, T???)
    @test DT.graph(H) == DT.graph(HH)
end

###################################################
#/
#/
#/ Storing a Delaunay triangulation 
#/
#/
###################################################
@testset "Testing from the Bowyer-Watson algorithm" begin
    pts = rand(2, 500)
    Random.seed!(2929292)
    T, adj, adj2v, DG = triangulate_bowyer(pts)
    tri = Triangulation(T, adj, adj2v, DG, pts)

    @test DT.get_triangles(tri) === T
    @test DT.get_adjacent(tri) === adj
    @test DT.get_adjacent2vertex(tri) === adj2v
    @test DT.get_graph(tri) === DG
    @test DT.get_points(tri) === pts

    _T, _adj, _adj2v, _DG, _pts = tri
    @test _T === T
    @test _adj === adj
    @test _adj2v === adj2v
    @test _DG === DG
    @test _pts === pts

    Random.seed!(2929292)
    _tri = triangulate_bowyer(pts)
    @test _tri isa Triangulation
    @test _tri.triangles == T
    @test _tri.adjacent.adjacent == adj.adjacent
    @test _tri.adjacent2vertex.adjacent2vertex == adj2v.adjacent2vertex
    @test _tri.graph.graph == DG.graph
    @test _tri.points == pts
end

@testset "Testing from de Berg's method" begin
    pts = rand(2, 500)
    Random.seed!(29292222292)
    (T, adj, adj2v, DG), HG = triangulate_berg(pts)
    @test HG isa DT.HistoryGraph
    tri = Triangulation(T, adj, adj2v, DG, pts)

    @test DT.get_triangles(tri) === T
    @test DT.get_adjacent(tri) === adj
    @test DT.get_adjacent2vertex(tri) === adj2v
    @test DT.get_graph(tri) === DG
    @test DT.get_points(tri) === pts

    _T, _adj, _adj2v, _DG, _pts = tri
    @test _T === T
    @test _adj === adj
    @test _adj2v === adj2v
    @test _DG === DG
    @test _pts === pts

    Random.seed!(29292222292)
    _tri, _HG = triangulate_berg(pts)
    @test _tri isa Triangulation
    @test _tri.triangles == T
    @test _tri.adjacent.adjacent == adj.adjacent
    @test _tri.adjacent2vertex.adjacent2vertex == adj2v.adjacent2vertex
    @test _tri.graph.graph == DG.graph
    @test _tri.points == pts
    @test HG.graph == _HG.graph
end