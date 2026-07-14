@testset "Table 2 Rows 1--7: 100 seeded defining-operator quadratures per formula" begin
    rng = MersenneTwister(SEED)
    for row in 1:7
        @testset "Row $row" begin
            accepted = 0
            attempts = 0
            while accepted < POINTS_PER_FORMULA && attempts < 20POINTS_PER_FORMULA
                attempts += 1
                s = accepted < POINTS_PER_FORMULA÷2 ?
                    -(rand(rng)*0.07+0.10) : rand(rng)*0.18+0.20
                exterior = accepted >= POINTS_PER_FORMULA÷2 && iseven(accepted)
                xmag = accepted < POINTS_PER_FORMULA÷2 ? rand(rng)*0.58+0.16 :
                    (exterior ? rand(rng)*0.65+1.22 : rand(rng)*0.70+0.14)
                x = isodd(accepted) ? -xmag : xmag
                n = row <= 4 ? rand(rng,0:(exterior ? 1 : 3)) : rand(rng,0:4)
                a = rand(rng)*0.48+0.52
                b = rand(rng)*0.26+0.30
                lambda = rand(rng)*0.42+0.64
                alpha = rand(rng)*0.51+0.18
                f = formula_1d(row,x;s,n,a,b,lambda,alpha)
                safely_nonconfluent(f) || continue
                actual = real(evaluate(f))
                reference = if s < 0
                    riesz_source_quadrature(row,x;t=-s,n,a,b,lambda,alpha)
                else
                    reference_at_cutoff = fractional_laplacian_source_quadrature(
                        row,x;s,n,a,b,lambda,alpha,tail_cutoff=256.0)
                    doubled = fractional_laplacian_source_quadrature(
                        row,x;s,n,a,b,lambda,alpha,tail_cutoff=512.0)
                    if row > 5
                        @test relative_error(reference_at_cutoff,doubled) < 2e-8
                    end
                    doubled
                end
                if isfinite(abs(actual)) && isfinite(abs(reference))
                    @test relative_error(actual, reference) < 8e-6
                    accepted += 1
                end
            end
            @test accepted == POINTS_PER_FORMULA
        end
    end
end

@testset "Table 2 Rows 8--14: seeded fractional-Laplacian singular integrals" begin
    rng = MersenneTwister(SEED + 2)
    for row in 8:14
        @testset "Row $row" begin
            for i in 1:OSCILLATORY_OPERATOR_POINTS
                s = rand(rng)*0.045+0.44
                xmag = rand(rng)*1.18+0.24
                x = isodd(i) ? -xmag : xmag
                a = rand(rng)*0.44+0.53
                nu = rand(rng)*0.28+0.21
                mu = rand(rng)*0.20+0.16
                formula = formula_1d(row,x;s,a,nu,mu)
                actual = real(evaluate(formula))
                reference_at_cutoff = fractional_laplacian_source_quadrature(
                    row,x;s,a,nu,mu,tail_cutoff=256.0)
                reference = fractional_laplacian_source_quadrature(
                    row,x;s,a,nu,mu,tail_cutoff=512.0)
                @test relative_error(reference_at_cutoff,reference) < 5e-4
                @test relative_error(actual,reference) < 8e-4
            end
        end
    end
end
