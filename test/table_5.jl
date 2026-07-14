@testset "Table 5 one-dimensional special cases" begin
    rng = MersenneTwister(SEED + 3)
    for row in ("1*", "1**", "1***", "6*", "6**")
        @testset "Row $row" begin
            for i in 1:POINTS_PER_FORMULA
                xmag = i <= POINTS_PER_FORMULA÷2 ? rand(rng)*0.68+0.14 : rand(rng)*0.72+1.22
                x = isodd(i) ? -xmag : xmag
                a = rand(rng)*0.58+0.42
                s = rand(rng)*0.18+0.19
                n = row in ("1***", "6**") ? rand(rng, 1:4) : rand(rng, 0:3)
                actual = special_1d(row, x; n, a, s)
                reference = if row == "1*"
                    fractional_laplacian_source_quadrature(1,x;s,n,a=s)
                elseif row == "1**"
                    fractional_laplacian_source_quadrature(1,x;s=0.5,n,a)
                elseif row == "1***"
                    logarithmic_potential_row1(x; a, n)
                elseif row == "6*"
                    at_cutoff = fractional_laplacian_source_quadrature(
                        6,x;s=0.5,n,tail_cutoff=20.0)
                    doubled = fractional_laplacian_source_quadrature(
                        6,x;s=0.5,n,tail_cutoff=40.0)
                    @test relative_error(at_cutoff,doubled) < 2e-10
                    doubled
                else
                    logarithmic_potential_hermite(x; n)
                end
                @test relative_error(actual, real(reference)) < 2e-7
            end
        end
    end
end
