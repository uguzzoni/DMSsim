for testrep=1:5
    Random.seed!(24482708+testrep);
    V = 2; T = 3;

    #= I had to simplify these tests due to a bug
    in Distributions.
    TODO: increase system size =#

    S = 30
    p = 0.1rand(S)

    N0 = rand(10^10 : 10^13)
    N,n = dms_display(V, T, N0, p)

    @test all(sum(N; dims=(1,)) .== sum(N[:,1,1]) .== N0 * S)
    @test all(N .≥ 0)
    @test all(0 .≤ n .≤ N[:,:,1:T-1])

    for t = 1:T-1, v = 1:V, s = 1:S
        if iszero(N[s,v,t])
            @test iszero(N[s,v,t+1])
        end
    end

    θ = vec(mean(n ./ max.(N[:,:,1:T-1],1); dims=(2,3)))

    @test cor(p,θ) ≥ 0.9999
end
