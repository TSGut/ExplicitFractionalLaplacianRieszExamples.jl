@testset "Table 6: 100 seeded defining-operator quadratures per formula" begin
    rng = MersenneTwister(SEED + 1)
    for row in (:A, :B, :C, :D)
        @testset "Row $row" begin
            for i in 1:POINTS_PER_FORMULA
                d = rand(rng,2:5)
                ell = row in (:C,:D) ? rand(rng,0:1) : rand(rng,0:2)
                n = row == :C ? 0 : rand(rng,0:3)
                r = row in (:C,:D) ? rand(rng)*0.58+0.18 :
                    (isodd(i) ? rand(rng)*0.62+0.16 : rand(rng)*0.60+1.24)
                b = row in (:C,:D) ? rand(rng)*0.05+0.05 : rand(rng)*0.18+0.08
                use_fractional = row != :B && iseven(i)
                if row in (:C,:D)
                    a = rand(rng)*0.04-0.92
                else
                    a = rand(rng)*0.42+0.48
                end
                s = use_fractional ? rand(rng)*0.07+0.23 :
                    -(rand(rng)*0.04+0.18)
                alpha_lower = row == :B ? max(0.18,ell/2-s-n-0.95) : 0.18
                alpha = alpha_lower+rand(rng)*0.48
                if row == :A
                    @test a > -1 && b < (d+ell)/2
                elseif row == :B
                    @test alpha > max(ell/2-s-n-1,-n-1)
                elseif row == :C
                    @test -1 < a < s-n-ell/2
                    @test b > -n-1+ell/2-s
                else
                    @test a > -1 && max(a,a+b) < s-ell/2
                end
                actual = real(evaluate(formula_highd(row,r;s,d,ell,n,a,b,alpha,
                                       solid_harmonic=r^ell)))
                if use_fractional
                    cutoff = row in (:C,:D) ? 65536.0 : 4096.0
                    reference_at_cutoff = fractional_laplacian_highd_separated_quadrature(
                        row,r;s,d,ell,n,a,b,alpha,tail_cutoff=cutoff)
                    reference = fractional_laplacian_highd_separated_quadrature(
                        row,r;s,d,ell,n,a,b,alpha,tail_cutoff=2cutoff)
                    if row in (:C,:D)
                        @test relative_error(reference_at_cutoff,reference) < 8e-6
                    end
                else
                    cutoff = row == :B ? 8.0 : row in (:C,:D) ? 65536.0 : 4096.0
                    reference_at_cutoff = riesz_highd_quadrature(
                        row,r;t=-s,d,ell,n,a,b,alpha,tail_cutoff=cutoff)
                    reference = riesz_highd_quadrature(
                        row,r;t=-s,d,ell,n,a,b,alpha,tail_cutoff=2cutoff)
                    if row != :A
                        @test relative_error(reference_at_cutoff,reference) < 8e-6
                    end
                end
                @test relative_error(actual,reference) < 1.5e-5
            end
        end
    end
end

@testset "Table 6 previously omitted operator branches" begin
    @testset "Row A fractional interior and Riesz exterior" begin
        d, ell, n, a, b, alpha = 3, 1, 1, 0.72, 0.20, 0.30
        @test a > -1 && b < (d+ell)/2

        r, s = 0.55, 0.27
        actual = real(evaluate(formula_highd(:A,r;s,d,ell,n,a,b,alpha,
                               solid_harmonic=r^ell)))
        reference = fractional_laplacian_highd_quadrature(
            :A,r;s,d,ell,n,a,b,alpha,rtol=3e-5,tail_cutoff=16.0)
        @test relative_error(actual,reference) < 8e-5

        r, t = 1.50, 0.20
        actual = real(evaluate(formula_highd(:A,r;s=-t,d,ell,n,a,b,alpha,
                               solid_harmonic=r^ell)))
        reference = riesz_highd_quadrature(:A,r;t,d,ell,n,a,b,alpha)
        @test relative_error(actual,reference) < 2e-6
    end

    @testset "Row B positive fractional order" begin
        r, s, d, ell, n = 0.65, 0.27, 2, 0, 1
        a, b, alpha = 0.72, 0.20, 0.30
        @test alpha > max(ell/2-s-n-1,-n-1)
        actual = real(evaluate(formula_highd(:B,r;s,d,ell,n,a,b,alpha)))
        reference_at_cutoff = fractional_laplacian_highd_quadrature(
            :B,r;s,d,ell,n,a,b,alpha,rtol=3e-5,tail_cutoff=16.0)
        reference = fractional_laplacian_highd_quadrature(
            :B,r;s,d,ell,n,a,b,alpha,rtol=3e-5,tail_cutoff=32.0)
        @test relative_error(reference_at_cutoff,reference) < 2e-6
        @test relative_error(actual,reference) < 8e-5
    end

    @testset "Rows C--D Riesz exterior" begin
        r, t, d, ell, n = 1.50, 0.20, 2, 0, 0
        a, b, alpha = -0.70, 0.10, 0.30
        for row in (:C,:D)
            if row == :C
                @test -1 < a < -t-n-ell/2
                @test b > -n-1+ell/2+t
            else
                @test a > -1 && max(a,a+b) < -t-ell/2
            end
            actual = real(evaluate(formula_highd(row,r;s=-t,d,ell,n,a,b,alpha)))
            reference_at_cutoff = riesz_highd_quadrature(
                row,r;t,d,ell,n,a,b,alpha,rtol=1e-7,tail_cutoff=65536.0)
            reference = riesz_highd_quadrature(
                row,r;t,d,ell,n,a,b,alpha,rtol=1e-7,tail_cutoff=131072.0)
            @test relative_error(reference_at_cutoff,reference) < 8e-6
            @test relative_error(actual,reference) < 1.5e-5
        end
    end

    @testset "Row C positive order with n > 0" begin
        r, s, d, ell, n = 0.55, 0.80, 2, 0, 1
        a, b, alpha = -0.90, 0.10, 0.30
        @test -1 < a < s-n-ell/2
        @test b > -n-1+ell/2-s
        actual = real(evaluate(formula_highd(:C,r;s,d,ell,n,a,b,alpha)))
        reference_at_cutoff = fractional_laplacian_highd_separated_quadrature(
            :C,r;s,d,ell,n,a,b,alpha,rtol=1e-7,tail_cutoff=4096.0)
        reference = fractional_laplacian_highd_separated_quadrature(
            :C,r;s,d,ell,n,a,b,alpha,rtol=1e-7,tail_cutoff=8192.0)
        @test relative_error(reference_at_cutoff,reference) < 8e-6
        @test relative_error(actual,reference) < 1.5e-5
    end

    @testset "Rows C--D positive-order exterior branch" begin
        cases = ((:C,0,1,0.20,0.10),
                 (:C,1,0,-0.50,0.10),
                 (:D,2,1,0.20,0.10))
        r, s, d, alpha = 1.45, 1.0, 3, 0.30
        for (row,n,ell,a,b) in cases
            if row == :C
                @test -1 < a < s-n-ell/2
                @test b > -n-1+ell/2-s
            else
                @test a > -1 && max(a,a+b) < s-ell/2
            end
            actual = real(evaluate(formula_highd(row,r;s,d,ell,n,a,b,alpha,
                                   solid_harmonic=r^ell)))
            reference = classical_laplacian_highd(row,r;d,ell,n,a,b,alpha)
            @test relative_error(actual,reference) < 2e-7
        end
    end
end
