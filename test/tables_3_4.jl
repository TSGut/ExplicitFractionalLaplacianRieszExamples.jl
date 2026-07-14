@testset "Tables 3--4 literal hypergeometric formulas, Rows 1--14" begin
    rng = MersenneTwister(SEED + 7)
    for row in 1:14
        @testset "Row $row" begin
            accepted = 0
            while accepted < POINTS_PER_FORMULA
                xmag = row <= 5 && accepted >= POINTS_PER_FORMULA÷2 ?
                    rand(rng)*0.65+1.22 : rand(rng)*0.68+0.14
                x = isodd(accepted) ? -xmag : xmag
                s = row <= 7 && isodd(fld(accepted,25)) ?
                    -(rand(rng)*0.09+0.08) : rand(rng)*0.15+0.19
                n = rand(rng, 0:4)
                a = rand(rng)*0.54+0.44
                b = rand(rng)*0.31+0.22
                lambda = rand(rng)*0.48+0.64
                alpha = rand(rng)*0.62+0.16
                nu = rand(rng)*0.31+0.47
                mu = rand(rng)*0.20+0.17
                safe_trig = minimum(abs.((sin(pi*(nu-2s)/2),
                            sin(pi*(s-nu/2)), cos(pi*(s-nu/2)), sin(pi*nu),
                            sin(pi*nu/2+pi*s), cos(pi*nu/2+pi*s),
                            sin(pi*(mu+nu-2s)/2)))) > 0.10
                safe_trig || continue
                actual = hypergeometric_1d(row, x;
                         s, n, a, b, lambda, alpha, nu, mu)
                reference = if row <= 4 && abs(x) > 1
                    s > 0 ?
                        fractional_laplacian_source_exterior(row,x;s,n,a,b,lambda) :
                        riesz_source_quadrature(row,x;t=-s,n,a,b,lambda,alpha)
                else
                    evaluate(formula_1d(row, x;
                             s, n, a, b, lambda, alpha, nu, mu))
                end
                if isfinite(abs(actual)) && isfinite(abs(reference))
                    tolerance = row <= 4 && abs(x) > 1 ? 5e-6 : 3e-8
                    @test relative_error(actual, real(reference)) < tolerance
                    accepted += 1
                end
            end
        end
    end
end
