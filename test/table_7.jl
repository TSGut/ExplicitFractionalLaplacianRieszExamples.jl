@testset "Table 7 higher-dimensional special cases" begin
    rng = MersenneTwister(SEED + 4)
    for row in ("A*", "A**", "A***")
        @testset "Row $row" begin
            for i in 1:POINTS_PER_FORMULA
                r = i <= POINTS_PER_FORMULA÷2 ? rand(rng)*0.68+0.14 : rand(rng)*0.62+1.24
                t = rand(rng)*0.04+0.18
                s = -t
                d = rand(rng, 2:5)
                ell = rand(rng, 0:1)
                n = row == "A***" ? 0 : rand(rng, 0:3)
                a = rand(rng)*0.55+0.45
                b = (d+2ell-2)/2
                actual = special_highd(row, r; s, d, ell, n, a,
                                       solid_harmonic=r^ell)
                parent_a = row == "A*" ? a : row == "A**" ? s : 0.0
                parent_n = row == "A***" ? 0 : n
                @test parent_a > -1 && b < (d+ell)/2
                reference = riesz_highd_quadrature(:A,r;t,d,ell,n=parent_n,
                                                    a=parent_a,b)
                @test relative_error(actual,reference) < 2e-6
            end
        end
    end
end

@testset "Table 7 Row A** at the numerical-example Riesz order" begin
    rng = MersenneTwister(SEED + 11)
    s, d = -1/3, 2
    for i in 1:20
        r = i <= 10 ? rand(rng)*0.66+0.16 : rand(rng)*0.58+1.26
        ell = isodd(i) ? 3 : rand(rng,0:1)
        n = rand(rng,0:3)
        b = (d+2ell-2)/2
        t, D = -s, d/2+ell
        @test D+n > t+ell/2
        @test 0 < (d+ell)/2
        @test 1+s > 1-2t
        actual = special_highd("A**",r;s,d,ell,n,solid_harmonic=r^ell)
        reference = riesz_highd_quadrature(:A,r;t=-s,d,ell,n,a=s,b)
        @test relative_error(actual,reference) < 2e-6
    end
end

@testset "Table 7 reductions against operator-validated Table 6 formulas" begin
    rng = MersenneTwister(SEED + 6)
    for row in (:A, :B)
        @testset "Row $row" begin
            for i in 1:POINTS_PER_FORMULA
                r = row == :A && i > POINTS_PER_FORMULA÷2 ?
                    rand(rng)*0.62+1.24 : rand(rng)*0.68+0.14
                s = rand(rng)*0.16+0.20
                d = rand(rng, 2:5)
                ell = rand(rng, 0:2)
                n = rand(rng, 0:3)
                a = rand(rng)*0.52+0.46
                b = rand(rng)*0.31+0.17
                alpha = rand(rng)*0.61+0.14
                actual = hypergeometric_highd(row, r;
                         s, d, ell, n, a, b, alpha)
                if row == :A
                    @test a > -1 && b < (d+ell)/2
                else
                    @test alpha > max(ell/2-s-n-1,-n-1)
                end
                reference = evaluate(formula_highd(row, r;
                            s, d, ell, n, a, b, alpha))
                @test relative_error(actual, real(reference)) < 2e-8
            end
        end
    end
end
