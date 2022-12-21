using LinearAlgebra

@testset "Computation of the edge-radius ratio" begin
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
            θ = min(θ1, θ2, θ3)
            ρ = DelaunayTriangulation.edge_radius_ratio(V, pts)
            ρℓ = norm(e1) * norm(e2) * norm(e3) / (4 * DelaunayTriangulation.area((pᵢ, pⱼ, pₖ))) / min(norm(e1), norm(e2), norm(e3))
            @test ρ ≈ ρℓ
        end
    end
end

@testset "Testing if a point encroaches on an edge" begin
    for _ in 1:500
        pts = rand(2, 500)
        T, adj, adj2v, DG = triangulate_bowyer(pts)
        for e in edges(DG)
            u, v = e
            if DelaunayTriangulation.BoundaryIndex ∉ e
                pu, pv = pts[:, u], pts[:, v]
                mid = (pu + pv) / 2
                rad = norm(pu - pv) / 2
                for _ in 1:500
                    q = rand(2)
                    encroaches = norm(q - mid) < rad
                    @test DT.encroaches_upon(q, e, pts) == encroaches
                end
            end
        end
    end
end

@testset "Testing if a segment is encroached" begin
    pts = rand(2, 30)
    T, adj, adj2v, DG = triangulate_bowyer(pts)
    for e in edges(DG)
        if !DT.is_ghost_edge(e...)
            encroached = DT.segment_encroached(e, adj, pts)
            k, ℓ = DT.get_edge(adj, e[1], e[2]), DT.get_edge(adj, e[2], e[1])
            i, j = e
            pᵢ, pⱼ = pts[:, i], pts[:, j]
            if k ≠ DT.BoundaryIndex && k ≠ DT.DefaultAdjacentValue && ℓ ≠ DT.BoundaryIndex && ℓ ≠ DT.DefaultAdjacentValue
                pₖ, pₗ = pts[:, k], pts[:, ℓ]
                pjk = pₖ - pⱼ
                pki = pᵢ - pₖ
                piℓ = pₗ - pᵢ
                pℓj = pⱼ - pₗ
                θij = acos(dot(pjk, pki) / (norm(pjk) * norm(pki)))
                θji = acos(dot(piℓ, pℓj) / norm(piℓ) * norm(pℓj))
                if θij > π / 2 || θji > π / 2
                    @test encroached
                else
                    @test !encroached
                end
            elseif k ≠ DT.BoundaryIndex && k ≠ DT.DefaultAdjacentValue
                pₖ = pts[:, k]
                pjk = pₖ - pⱼ
                pki = pᵢ - pₖ
                θij = acos(dot(pjk, pki) / (norm(pjk) * norm(pki)))
                if θij > π / 2
                    @test encroached
                else
                    @test !encroached
                end
            elseif ℓ ≠ DT.BoundaryIndex && ℓ ≠ DT.DefaultAdjacentValue
                pₗ = pts[:, ℓ]
                piℓ = pₗ - pᵢ
                pℓj = pⱼ - pₗ
                θji = acos(dot(piℓ, pℓj) / norm(piℓ) * norm(pℓj))
                if θji > π / 2
                    @test encroached
                else
                    @test !encroached
                end
            end
        end
    end
    for _ in 1:500
        i, j = rand(1:30, 2)
        e = (i, j)
        encroached = DT.segment_encroached(e, adj, pts)
        if !DT.edge_exists(i, j, adj)
            @test encroached
        else
            k, ℓ = DT.get_edge(adj, e[1], e[2]), DT.get_edge(adj, e[2], e[1])
            i, j = e
            pᵢ, pⱼ = pts[:, i], pts[:, j]
            if k ≠ DT.BoundaryIndex && k ≠ DT.DefaultAdjacentValue && ℓ ≠ DT.BoundaryIndex && ℓ ≠ DT.DefaultAdjacentValue
                pₖ, pₗ = pts[:, k], pts[:, ℓ]
                pjk = pₖ - pⱼ
                pki = pᵢ - pₖ
                piℓ = pₗ - pᵢ
                pℓj = pⱼ - pₗ
                θij = acos(dot(pjk, pki) / (norm(pjk) * norm(pki)))
                θji = acos(dot(piℓ, pℓj) / norm(piℓ) * norm(pℓj))
                if θij > π / 2 || θji > π / 2
                    @test encroached
                else
                    @test !encroached
                end
            elseif k ≠ DT.BoundaryIndex && k ≠ DT.DefaultAdjacentValue
                pₖ = pts[:, k]
                pjk = pₖ - pⱼ
                pki = pᵢ - pₖ
                θij = acos(dot(pjk, pki) / (norm(pjk) * norm(pki)))
                if θij > π / 2
                    @test encroached
                else
                    @test !encroached
                end
            elseif ℓ ≠ DT.BoundaryIndex && ℓ ≠ DT.DefaultAdjacentValue
                pₗ = pts[:, ℓ]
                piℓ = pₗ - pᵢ
                pℓj = pⱼ - pₗ
                θji = acos(dot(piℓ, pℓj) / norm(piℓ) * norm(pℓj))
                if θji > π / 2
                    @test encroached
                else
                    @test !encroached
                end
            end
        end
    end
end

@testset "Skinny and oversized triangles" begin
    A = [0.727775900868, 2.4317500969586]
    B = [1.181022874267, 2.404192847742]
    C = [1.64, 2.44] # needle
    T = (1, 2, 3)
    pts = [A, B, C]
    ρ = DT.edge_radius_ratio(T, pts)
    @test ρ ≈ 7.2777532478509 rtol = 1e-3
    @test DT.is_skinny(ρ, √2)
    @test DT.is_skinny(ρ, 1.0)
    @test !DT.is_skinny(ρ, 6.2)
    @test DT.is_skinny(ρ, 3.6)
    @test !DT.is_skinny(ρ, 3.7)
    @test DT.is_oversized(ρ, 3.7)
    @test DT.is_oversized(ρ, 7.0)
    pts[2] = [0.7315173512548, 2.4081172794674] # cap
    ρ = DT.edge_radius_ratio(T, pts)
    @test ρ ≈ 19.2059067379818
    @test DT.is_skinny(ρ, √2)
    @test DT.is_skinny(ρ, 8.0)
    @test !DT.is_skinny(ρ, 10.0)
    @test !DT.is_oversized(ρ, 5.0)
    @test DT.is_oversized(ρ, 15.0)
    pts[2] = [0.9849691359043, 2.2425560330431] # bit large
    ρ = DT.edge_radius_ratio(T, pts)
    @test ρ ≈ 1.786156210909
    @test !DT.is_skinny(ρ, √2)
    @test DT.is_oversized(ρ, √2)
    @test !DT.is_oversized(ρ, 2.0)
    pts[2] = [1.1965196174465, 2.054511605612] # standard
    ρ = DT.edge_radius_ratio(T, pts)
    @test ρ ≈ 0.788660948192 rtol = 1e-3
    @test !DT.is_skinny(ρ, √2)
    @test !DT.is_oversized(ρ, √2)
end