@testset "Table 1 classical-function Meijer-G representations" begin
    rng = MersenneTwister(SEED + 5)
    for row in 1:14
        @testset "Row $row" begin
            for i in 1:POINTS_PER_FORMULA
                xmag = row <= 5 ? rand(rng)*0.74+0.12 : rand(rng)*1.75+0.16
                x = isodd(i) ? -xmag : xmag
                n = rand(rng, 0:5)
                a = rand(rng)*0.62+0.34
                b = rand(rng)*0.38+0.21
                lambda = rand(rng)*0.52+0.61
                alpha = rand(rng)*0.72+0.12
                nu = rand(rng)*0.38+0.19
                mu = rand(rng)*0.31+0.16
                actual = evaluate(source_formula_1d(row, x;
                                  n, a, b, lambda, alpha, nu, mu))
                reference = source_classical_1d(row, x;
                            n, a, b, lambda, alpha, nu, mu)
                @test relative_error(real(actual), reference) < 2e-8
            end
        end
    end
end
