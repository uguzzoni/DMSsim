
"""
    dms!(N, n, p)
    Simulates a DMS experiment, where N and n are S×V×T and S×V×(T-1) tensors, 
    where S is the number of variants in the experiment, V the replicates, T the rounds, 
    N[s,v,t] is the number of variants s in the population in the experiment replicate v at round t, 
    and n[s,v,t] the number of selected variants. Finally p[s] is the selectivity of sequence s.

    For each replicate uses the initial condition N[:,v,1]. For each round and replicate sets Ntot[v,t] from sum(N[:,v,t]). 
    Overwrites N[:,:,2:T] and n[:,:,:] with the simulation results.
    
"""
function dms_display!(N::AbstractArray{<:Integer,3},
                        n::AbstractArray{<:Integer,3},
                        p::AbstractVector{<:Real})
    S,V,T = size(N)
    @assert length(p) == S
    @assert size(n) == (S,V,T-1)
    for s = 1:S @assert 0 ≤ p[s] ≤ 1 end
    for v = 1:V, s = 1:S @assert N[s,v,1] ≥ 0 end
    for v = 1:V, t = 1:T-1
        dms_select!(view(n,:,v,t), view(N,:,v,t), p)
        dms_amplify!(view(N,:,v,t+1), view(n,:,v,t))
    end
    nothing
end


"""
    dms_select!(n, N, p)

simulates a selection step
"""
function dms_select!(n::AbstractVector{<:Integer},
                     N::AbstractVector{<:Integer},
                     p::AbstractVector{<:Real})
    @assert length(n) == length(N) == length(p)
    for x in p @assert 0 ≤ x ≤ 1 end
    for x in N @assert x ≥ 0 end
    n .= rand.(Binomial.(N, p))
    nothing
end

"""
    dms_amplify!(N, n)

simulates amplification of selected variants in dms (n)
up to a population of size Ntot = sum(N).
"""
function dms_amplify!(N::AbstractVector{<:Integer},
                        n::AbstractVector{<:Integer})
    @assert length(N) == length(n)
    for x in n @assert x ≥ 0 end
    for x in N @assert x ≥ 0 end
    ntot = sum(n)
    @assert ntot ≥ 0
    if ntot > 0
        Ntot = sum(N)
        @assert Ntot ≥ 0
        N .= rand(Multinomial(Ntot, n ./ ntot))
    else
        N .= 0
    end
    nothing
end


function dms_select_competitive(n::AbstractVector{Int},
                                N::AbstractVector{Int})
end


"""
Simulates a dms display experiment.

    dms_display(V, T, N0, p)

where V is replicates, T rounds, N0 initial counts
of every sequence at first round of every replicate,
and p the vector of binding probabilities.
"""
function dms_display(V::Int, T::Int, N0::AbstractVector{Int},
                       p::AbstractVector{Float64})
    @assert V ≥ 0 && T ≥ 0
    @assert length(p) == length(N0)
    for x in N0 @assert x ≥ 0 end
    for x in p @assert 0 ≤ x ≤ 1 end

    S = length(p)
    N=fill(0,S,V,T)
    for v=1:V
        for t=1:T
            N[:,v,t]=N0
        end
    end
    n = zeros(Int, S, V, T-1)
    dms_display!(N, n, p)
    N, n
end

function dms_display(V::Int, T::Int, N0::Int,
                       p::AbstractVector{Float64})
    @assert N0 ≥ 0 && V ≥ 0 && T ≥ 0
    for x in p @assert 0 ≤ x ≤ 1 end
    S = length(p)
    N = fill(N0, S, V, T)
    n = zeros(Int, S, V, T-1)
    dms_display!(N, n, p)
    N, n
end



"""
    dms_reads(Rt, N)

Simulates taking reads from each sample with dms abundances
N, where Rt is the total number of reads per sample.
"""
function dms_reads(Rt::Integer, N::AbstractArray{<:Real,3})
    @assert Rt ≥ 0
    S, V, T = size(N)
    R = zeros(typeof(Rt), size(N))
    R[1,:,:] .= Rt
    return dms_reads!(R, N)
end


"""
    dms_reads!(R, N)

Simulates taking reads from each sample with dms abundances
N, where Rtot is the total number of reads per sample.
"""
function dms_reads!(R::AbstractArray{<:Integer,3}, N::AbstractArray{<:Real,3})
    @assert size(N) == size(R)
    for x in N @assert 0 ≤ x < Inf end
    for x in R @assert 0 ≤ x < Inf end
    S, V, T = size(N)
    for t = 1:T, v = 1:V
        Nt = sum(view(N,:,v,t))
        Rt = sum(view(R,:,v,t))
        @assert Nt ≥ 0
        @assert Rt ≥ 0
        if Nt > 0 && Rt > 0
            R[:,v,t] .= rand(Multinomial(Rt, N[:,v,t] ./ Nt))
        else
            R[:,v,t] .= 0
        end
    end
    return R
end



"""
    counts_round0(S, Ntot; distribution)

dms abundances at round 0
S number of unique sequences
Ntot=sum(N) total number of dmss in the experiment after amplification.
"""


#constructing round 0 counts
# very ad hoc manner, Gaussian  mean Ntot/S, correcting negative and zero counts
function counts_round0(S,Ntot;distribution::Symbol=:gaussian,sigma::Float64=1)

    if distribution == :gaussian
        mean_N0=Ntot/S
        sigma_n0=sqrt(mean_N0)*sigma
        N0=convert.(Int64,trunc.((randn(S).*sigma_n0) .+ mean_n0))
    elseif distribution == :flat
        N0=convert.(Int64,trunc.(Ntot/S))
    else
        error("distribution not included")
    end

    return N0

end
