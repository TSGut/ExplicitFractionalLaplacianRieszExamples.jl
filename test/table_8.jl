@testset "Table 8 generic Meijer-G reduction identities" begin
    rng = MersenneTwister(SEED + 9)
    orders = ((2,2,1,1),(2,3,1,2),(2,3,2,1),(2,4,2,1),(3,3,2,1))
    for (p,q,m,nleft) in orders
        @testset "G($p,$q)^($m,$nleft)" begin
            for i in 1:POINTS_PER_FORMULA
                aa = collect(range(1.21,2.03,length=p)) .+ rand(rng,p).*0.09
                bb = collect(range(0.11,0.63,length=q)) .+ rand(rng,q).*0.07
                z = p == q && i > POINTS_PER_FORMULA÷2 ?
                    rand(rng)*0.64+1.21 : rand(rng)*0.68+0.14
                f = PaperGFormula("Table 8",1.0,aa,bb,m,nleft,z)
                @test relative_error(evaluate(f),residue_reference(f)) < 2e-8
            end
        end
    end
end
