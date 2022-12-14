include("functions.jl")

pts, T, DG, adj, adj2v = example_triangulation()
p1, p2, p3, p4, p5, p6 = pts
p7 = @SVector[2.0, 1.0]
DT.add_point!(pts, p7)
DT.add_triangle!(6, 2, 3, T, adj, adj2v, DG)
DT.split_triangle!(1, 3, 5, 7, T, adj, adj2v, DG)
DT.flip_edge!(1, 5, T, adj, adj2v, DG)
p8 = @SVector[1.0, 1.0]
DT.add_point!(pts, p8)
p9 = @SVector[2.5, 0.5]
DT.add_point!(pts, p9)
p10 = @SVector[1.5, 2.0]
DT.add_point!(pts, p10)

fig = Figure(fontsize=55)
ax = Axis(fig[1, 1])
xlims!(ax, -2.5, 4.5)
ylims!(ax, -1.5, 2.5)
triplot!(ax, T, pts; strokewidth=2, color=(:white, 0))
t1 = text!(ax, [(-0.05, 0.55), (3.1, -1.2), (2.2, -0.15),
        (-1.0, 2.02), (4.0, 2.02), (-2.4, -1.15), (2.0, 1.0),
        (0.7, 0.6), (2.65, 0.3), (1.5, 2.0)];
    text=[L"1", L"2", L"3", L"4", L"5", L"6", L"7", L"8", L"9", L"10"], textsize=55)
hidedecorations!(ax)
save("writeups/figures/split_edge_1.pdf", fig)

ℓ1 = lines!(ax, [p4, p8, p3], color=:blue, linewidth=4)
save("writeups/figures/split_edge_2.pdf", fig)

ℓ1.color = :black
ℓ1.linewidth = 2
ℓ2 = lines!(ax, [p7, p9, p2], color=:blue, linewidth=4)
save("writeups/figures/split_edge_3.pdf", fig)

ℓ2.color = :black
ℓ2.linewidth = 2
ℓ3 = lines!(ax, [p10, p7], color=:blue, linewidth=4)
save("writeups/figures/split_edge_4.pdf", fig)

## Setup 
pts, T, DG, adj, adj2v = example_triangulation()
p1, p2, p3, p4, p5, p6 = pts
p7 = @SVector[2.0, 1.0]
DT.add_point!(pts, p7)
DT.add_triangle!(6, 2, 3, T, adj, adj2v, DG)
DT.split_triangle!(1, 3, 5, 7, T, adj, adj2v, DG)
DT.flip_edge!(1, 5, T, adj, adj2v, DG)
p8 = @SVector[1.0, 1.0]
DT.add_point!(pts, p8)
p9 = @SVector[2.5, 0.5]
DT.add_point!(pts, p9)
p10 = @SVector[1.5, 2.0]
DT.add_point!(pts, p10)

## Interior edge 
i, j, r = 1, 7, 8
DT.split_edge!(i, j, r, T, adj, adj2v, DG)
DT.split_edge!(j, i, r, T, adj, adj2v, DG)
true_T = Set{NTuple{3,Int64}}([
    (3, 2, 5),
    (1, 8, 4),
    (7, 8, 3),
    (3, 5, 7),
    (6, 3, 1),
    (4, 6, 1),
    (4, 7, 5),
    (6, 2, 3),
    (8, 7, 4),
    (8, 1, 3),
])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (3, 2) => 5, (2, 5) => 3, (5, 3) => 2,
        (1, 8) => 4, (8, 4) => 1, (4, 1) => 8,
        (7, 8) => 3, (8, 3) => 7, (3, 7) => 8,
        (3, 5) => 7, (5, 7) => 3, (7, 3) => 5,
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
        (4, 7) => 5, (7, 5) => 4, (5, 4) => 7,
        (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
        (8, 7) => 4, (7, 4) => 8, (4, 8) => 7,
        (8, 1) => 3, (1, 3) => 8, (3, 8) => 1,
        (4, 5) => DT.BoundaryIndex,
        (5, 2) => DT.BoundaryIndex,
        (2, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex,
        (5, 1) => DT.DefaultAdjacentValue,
        (1, 7) => DT.DefaultAdjacentValue
    ))
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
    1 => Set{NTuple{2,Int64}}([(4, 6), (6, 3), (3, 8), (8, 4)]),
    2 => Set{NTuple{2,Int64}}([(5, 3), (3, 6)]),
    3 => Set{NTuple{2,Int64}}([(6, 2), (2, 5), (5, 7), (7, 8), (8, 1), (1, 6)]),
    4 => Set{NTuple{2,Int64}}([(6, 1), (1, 8), (8, 7), (7, 5)]),
    5 => Set{NTuple{2,Int64}}([(4, 7), (7, 3), (3, 2)]),
    6 => Set{NTuple{2,Int64}}([(2, 3), (3, 1), (1, 4)]),
    7 => Set{NTuple{2,Int64}}([(3, 5), (5, 4), (4, 8), (8, 3)]),
    8 => Set{NTuple{2,Int64}}([(1, 3), (3, 7), (7, 4), (4, 1)])
)
true_DG = UndirectedGraph(
    [
        0 0 1 1 0 1 0 1
        0 0 1 0 1 1 0 0
        1 1 0 0 1 1 1 1
        1 0 0 0 1 1 1 1
        0 1 1 1 0 0 1 0
        1 1 1 1 0 0 0 0
        0 0 1 1 1 0 0 1
        1 0 1 1 0 0 1 0
    ]
)
@test T == true_T
@test adjacent(adj) == true_adj
@test adjacent2vertex(adj2v) == true_adj2v
@test graph(DG) == true_DG

## Interior edge for a triangle that is on the boundary
i, j, r = 5, 3, 9
DT.split_edge!(i, j, r, T, adj, adj2v, DG)
DT.split_edge!(j, i, r, T, adj, adj2v, DG)
true_T = Set{NTuple{3,Int64}}([
    (1, 8, 4),
    (7, 8, 3),
    (6, 3, 1),
    (4, 6, 1),
    (4, 7, 5),
    (6, 2, 3),
    (8, 7, 4),
    (8, 1, 3),
    (9, 5, 7),
    (3, 9, 7),
    (9, 3, 2),
    (5, 9, 2)
])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (1, 8) => 4, (8, 4) => 1, (4, 1) => 8,
        (7, 8) => 3, (8, 3) => 7, (3, 7) => 8,
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
        (4, 7) => 5, (7, 5) => 4, (5, 4) => 7,
        (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
        (8, 7) => 4, (7, 4) => 8, (4, 8) => 7,
        (8, 1) => 3, (1, 3) => 8, (3, 8) => 1,
        (7, 9) => 5, (9, 5) => 7, (5, 7) => 9,
        (3, 9) => 7, (9, 7) => 3, (7, 3) => 9,
        (9, 3) => 2, (3, 2) => 9, (2, 9) => 3,
        (5, 9) => 2, (9, 2) => 5, (2, 5) => 9,
        (4, 5) => DT.BoundaryIndex,
        (5, 2) => DT.BoundaryIndex,
        (2, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex,
        (5, 1) => DT.DefaultAdjacentValue,
        (1, 7) => DT.DefaultAdjacentValue,
        (5, 3) => DT.DefaultAdjacentValue
    ))
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 5), (5, 2), (2, 6), (6, 4)]),
    1 => Set{NTuple{2,Int64}}([(4, 6), (6, 3), (3, 8), (8, 4)]),
    2 => Set{NTuple{2,Int64}}([(5, 9), (9, 3), (3, 6)]),
    3 => Set{NTuple{2,Int64}}([(2, 9), (9, 7), (7, 8), (8, 1), (1, 6), (6, 2)]),
    4 => Set{NTuple{2,Int64}}([(6, 1), (1, 8), (8, 7), (7, 5)]),
    5 => Set{NTuple{2,Int64}}([(4, 7), (7, 9), (9, 2)]),
    6 => Set{NTuple{2,Int64}}([(2, 3), (3, 1), (1, 4)]),
    7 => Set{NTuple{2,Int64}}([(9, 5), (5, 4), (4, 8), (8, 3), (3, 9)]),
    8 => Set{NTuple{2,Int64}}([(1, 3), (3, 7), (7, 4), (4, 1)]),
    9 => Set{NTuple{2,Int64}}([(2, 5), (5, 7), (7, 3), (3, 2)])
)
true_DG = UndirectedGraph(
    [
        0 0 1 1 0 1 0 1 0
        0 0 1 0 1 1 0 0 1
        1 1 0 0 0 1 1 1 1
        1 0 0 0 1 1 1 1 0
        0 1 0 1 0 0 1 0 1
        1 1 1 1 0 0 0 0 0
        0 0 1 1 1 0 0 1 1
        1 0 1 1 0 0 1 0 0
        0 1 1 0 1 0 1 0 0
    ]
)
@test T == true_T
@test adjacent(adj) == true_adj
@test adjacent2vertex(adj2v) == true_adj2v
@test graph(DG) == true_DG


## Boundary edge 
i, j, r = 5, 4, 10
DT.split_edge!(i, j, r, T, adj, adj2v, DG)
true_T = Set{NTuple{3,Int64}}([
    (1, 8, 4),
    (7, 8, 3),
    (6, 3, 1),
    (4, 6, 1),
    (6, 2, 3),
    (8, 7, 4),
    (8, 1, 3),
    (9, 5, 7),
    (3, 9, 7),
    (9, 3, 2),
    (5, 9, 2),
    (5, 10, 7),
    (10, 4, 7)
])
true_adj = DefaultDict(DT.DefaultAdjacentValue,
    Dict(
        (1, 8) => 4, (8, 4) => 1, (4, 1) => 8,
        (7, 8) => 3, (8, 3) => 7, (3, 7) => 8,
        (6, 3) => 1, (3, 1) => 6, (1, 6) => 3,
        (4, 6) => 1, (6, 1) => 4, (1, 4) => 6,
        (6, 2) => 3, (2, 3) => 6, (3, 6) => 2,
        (8, 7) => 4, (7, 4) => 8, (4, 8) => 7,
        (8, 1) => 3, (1, 3) => 8, (3, 8) => 1,
        (7, 9) => 5, (9, 5) => 7, (5, 7) => 9,
        (3, 9) => 7, (9, 7) => 3, (7, 3) => 9,
        (9, 3) => 2, (3, 2) => 9, (2, 9) => 3,
        (5, 9) => 2, (9, 2) => 5, (2, 5) => 9,
        (5, 10) => 7, (10, 7) => 5, (7, 5) => 10,
        (10, 4) => 7, (4, 7) => 10, (7, 10) => 4,
        (4, 10) => DT.BoundaryIndex,
        (10, 5) => DT.BoundaryIndex,
        (5, 2) => DT.BoundaryIndex,
        (2, 6) => DT.BoundaryIndex,
        (6, 4) => DT.BoundaryIndex,
        (5, 1) => DT.DefaultAdjacentValue,
        (1, 7) => DT.DefaultAdjacentValue,
        (5, 3) => DT.DefaultAdjacentValue
    ))
true_adj2v = Dict(
    DT.BoundaryIndex => Set{NTuple{2,Int64}}([(4, 10),(10,5), (5, 2), (2, 6), (6, 4)]),
    1 => Set{NTuple{2,Int64}}([(4, 6), (6, 3), (3, 8), (8, 4)]),
    2 => Set{NTuple{2,Int64}}([(5, 9), (9, 3), (3, 6)]),
    3 => Set{NTuple{2,Int64}}([(2, 9), (9, 7), (7, 8), (8, 1), (1, 6), (6, 2)]),
    4 => Set{NTuple{2,Int64}}([(6, 1), (1, 8), (8, 7), (7, 10)]),
    5 => Set{NTuple{2,Int64}}([(10, 7), (7, 9), (9, 2)]),
    6 => Set{NTuple{2,Int64}}([(2, 3), (3, 1), (1, 4)]),
    7 => Set{NTuple{2,Int64}}([(5, 10), (10, 4), (4, 8), (8, 3), (3, 9), (9, 5)]),
    8 => Set{NTuple{2,Int64}}([(1, 3), (3, 7), (7, 4), (4, 1)]),
    9 => Set{NTuple{2,Int64}}([(2, 5), (5, 7), (7, 3), (3, 2)]),
    10 => Set{NTuple{2,Int64}}([(4, 7), (7, 5)])
)
true_DG = UndirectedGraph(
    [
        0 0 1 1 0 1 0 1 0 0
        0 0 1 0 1 1 0 0 1 0 
        1 1 0 0 0 1 1 1 1 0
        1 0 0 0 0 1 1 1 0 1
        0 1 0 0 0 0 1 0 1 1
        1 1 1 1 0 0 0 0 0 0
        0 0 1 1 1 0 0 1 1 1
        1 0 1 1 0 0 1 0 0 0 
        0 1 1 0 1 0 1 0 0 0
        0 0 0 1 1 0 1 0 0 0
    ]
)
@test T == true_T
@test adjacent(adj) == true_adj
@test adjacent2vertex(adj2v) == true_adj2v
@test graph(DG) == true_DG

