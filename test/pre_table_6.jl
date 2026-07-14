@testset "Higher-dimensional Jacobi input identities before Table 6" begin
    rng = MersenneTwister(SEED + 8)
    for row in (:A,:Z,:C,:D)
        @testset "Input $row" begin
            for _ in 1:POINTS_PER_FORMULA
                r = row == :Z ? rand(rng)*0.48+0.14 :
                    row == :A ? rand(rng)*0.72+0.14 : rand(rng)*0.75+1.18
                n = row == :Z ? rand(rng,0:3) : rand(rng,0:4)
                a = rand(rng)*0.54+0.43
                b = rand(rng)*0.42+0.16
                D = rand(rng)*1.4+1.1
                formula = higher_input_formula(row,r;n,a,b,D)
                actual = row == :Z ? higher_input_z_reduction(r;n,a,D) : evaluate(formula)
                reference = higher_input_classical(row,r;n,a,b,D)
                @test relative_error(real(actual),reference) < 3e-8
            end
        end
    end
end
