#Function mutate! which takes as input a list of sequences with associated populations
#and performs mutations according to a probability of mutations p. Then updates new
#set of unique sequences and populations.

function mutate_dict!(seq::Array{Array{Int64,1},1}, N::AbstractArray{Int64,1}, L::Int64, p::Float64; A::Int64=20)

    S = size(N,1)
    @assert S==length(seq)
    b = Distributions.Binomial(L,p)
    nm = zeros(Int64,sum(N))
    offset::Int64=0

    dict_seq=Dict{Tuple, Int64}()
    sizehint!(dict_seq,round(Int,p*L*sum(N))+length(seq))
    for s=1:S
        dict_seq[tuple(seq[s]...)] = N[s]
    end

    for i=1:length(nm)
        nm[i] = rand(b)
    end

    new_s = zeros(Int64,L)
    @inbounds for s=1:S
        for n=1:N[s]
            #nm = rand(b)
            if nm[offset+n] > 0
                idx_m = sample(1:L,nm[offset+n],replace=false)
                new_s .= seq[s]
                @simd for i in idx_m
                    new_s[i] = mod( new_s[i] - 1 + rand(1:(A-1)), A) + 1 #mutation
                end
                dict_seq[tuple(seq[s]...)]-=1
                ns=tuple(new_s...)
                if haskey(dict_seq,ns)
                    dict_seq[ns] = dict_seq[ns]+1 #idx
                else
                    #println(new_s)
                    dict_seq[ns] = 1
                    push!(seq,collect(ns))
                end
            end
        end
        offset += N[s]
    end

    # ordered_dict = sort(dict_seq);
    # return ordered_dict;

     N_new = Array{Int64,1}(undef,length(dict_seq))

     counter = 0

    # #for key in sort!(collect(keys(dict_seq)))
    # for i in values(dict_seq)
    #      counter += 1
    #      N_new[counter] = i
    #  end

    for s in seq
        counter += 1
        N_new[counter] = dict_seq[tuple(s...)]
    end

    return N_new
end

#cumulative version: slower and wrong(?) number of unique sequences
function mutate_cum!(seq::Array{Array{Int64,1},1}, N::AbstractArray{Int64,1}, L::Int64, p::Float64; A::Int64=20)

    S = size(N,1)
    @assert S==length(seq)

    sumN = sum(N)

    dict_seq=Dict{Array{Int64,1}, Int64}()
    vec_cumAA=zeros(Int64,S)
    cumAA = 0
    for s=1:S
        cumAA += N[s]*L # cumulative number aa
        vec_cumAA[s]=cumAA
        dict_seq[seq[s]] = s
    end

    nmut = rand(Distributions.Binomial(L*sumN,p)) # total number of mutations
    aa_sub= rand(1:(A-1),nmut) #A-1 possible mutations
    rand_idx=sort(randperm(L*sumN)[1:nmut]) # mutation index

    new_s=Array{Int64,1}[]
    new_N=copy(N)

    s=1 # seq index
    i=1 # mut index
    n_news=0 #number new sequence (different from nmut if not only single mutations)
    while i <= nmut

            im=rand_idx[i]
            #select sequence to mutate
            while vec_cumAA[s] < im #cum
                s+=1
            end
            #
            new_s=copy(seq[s])
            new_N[s]-=1 #if zero should be removed

            #
            #more mutation on the same sequence?
            im_next=im
            while ( div(im_next-1,L) == div(im-1, L) && i <=nmut)
                im_next = rand_idx[i]
                pos = mod(im_next - vec_cumAA[s] ,L) + 1 #mutation position
                new_s[pos] = mod( new_s[pos] - 1 + aa_sub[i], A) + 1 #mutation
                i+=1
            end

            if haskey(dict_seq,new_s)
                new_N[dict_seq[new_s]]+=1 #idx
            else
                S+=1
                push!(seq,new_s)
                push!(new_N,1)
                dict_seq[new_s] = S
            end

    end

    vcat(new_N)

end

#Function dms__muts! which simulates a DMS experiment with mutations occuring
#after the amplification process. Updates also the list of unique sequences.

function dms_muts!(seq::Array{Array{Int64,1},1}, N0::AbstractArray{<:Integer,2}, p_m::Float64, T::Int64, p_b::AbstractVector{<:Real}, par, μ::Float64; b::Float64=1.0)

    V,S = size(N0)
    L = length(seq[1])
    S_new = S
    n = Array{Array{Int64,1}}(undef,V,T-1)
    N_new = Array{Array{Int64,1}}(undef,V,T)
    N_mut = Array{Array{Int64,1}}(undef,V)
    #N_amp = Array{Array{Int64,1}}(undef,V)
    for v=1:V
        N_new[v,1] = N0[v,:]
        N_mut[v] = N0[v,:]
    end
    @assert length(p_b) == S == length(seq)
    for s = 1:S @assert 0 ≤ p_b[s] ≤ 1 end
    for v = 1:V, s = 1:S @assert N_new[v,1][s] ≥ 0 end
    for v = 1:V, t = 1:T-1
        N_mut[v] = mutate_dict!(seq,N_new[v,t],L,p_m)
        for i=S_new+1:length(seq)
            push!(p_b,compute_prob_seq(seq[i],par,μ,beta=b))
        end
        S_new = length(seq)
        n[v,t] = zeros(Int64,S_new)
        dms_select!(n[v,t], N_mut[v], p_b)
        N_new[v,t+1] = dms_amplify(N_mut[v], n[v,t])
    end
    return N_new, n
end

#Function which simulates a MEv experiment, starting from a single WT sequence
function mev_sim!(wt_seq::Array{Int64,1}, N0::Array{Int64,1}, p_m::Float64, T::Int64, p_b::Array{Float64,1}, par, μ::Float64; b::Float64=1.0)

    V = length(N0)
    L = length(wt_seq)
    S_new = 1
    seq = [wt_seq]
    n = Array{Array{Int64,1}}(undef,V,T)
    N_new = Array{Array{Int64,1}}(undef,V,T)
    N_start = Array{Array{Int64,1}}(undef,V)
    n_surv = zeros(Int64,V,T)
    idx_surv = []


    for v=1:V
        N_start[v] = [N0[v]]
    end

    @assert length(p_b) == S_new == length(seq)
    @assert 0 ≤ p_b[1] ≤ 1
    for v = 1:V N0[v] ≥ 0 end
    for v = 1:V, t = 1:T
        N_new[v,t] = mutate_dict!(seq,N_start[v],L,p_m)
        idx_surv = collect(S_new+1:length(seq))
        for i=S_new+1:length(seq)
             push!(p_b,compute_prob_seq(seq[i],par,μ,beta=b))
        end
        n[v,t] = zeros(Int64,length(seq))
        dms_select!(n[v,t], N_new[v,t], p_b)
        for i in idx_surv
            if n[v,t][i]>0
                n_surv[v,t] +=1
            end
        end
        N_new[v,t] = dms_amplify(N_new[v,t], n[v,t])
        N_start[v] = N_new[v,t]
        S_new = length(seq)
    end
    return N_new, n, seq, n_surv
end

#Function dms_amplify not directly modifying population structure N at time t+1,
#providing instead updated population given total abudances at previous round.

function dms_amplify(N::AbstractArray{Int64,1}, n::AbstractArray{Int64,1})

    @assert length(N) == length(n)
    N_amp = zeros(Int64,length(N))
    for x in n @assert x ≥ 0 end
    for x in N @assert x ≥ 0 end
    ntot = sum(n)
    @assert ntot ≥ 0
    if ntot > 0
        Ntot = sum(N)
        @assert Ntot ≥ 0
        N_amp .= rand(Multinomial(Ntot, n ./ ntot))
    end

    return N_amp

end

function compute_prob_seq(seq::Array{Int64,1}, par, μ::Float64; beta::Float64=1.0)

    p=0.0
    E = compute_energy_seq(par.J_E,par.h_E,seq)
    φ = μ-E
    p = 1/(1+exp(-beta*φ))
    return p

end

function compute_energy_seq(J::Array{T,4},h::Array{T,2},s::Array{Int,1})  where T<:Real


    N = size(J,3)

    e=0.0
    for i in 1:N-1
        for j in i+1:N
            e-=J[s[i],s[j],i,j];
        end
    end

    for i in 1:N
        e-=h[s[i],i];
    end

    return e
end


function mut_freq(N::Array{Array{Int64,1}}, seq::Array{Array{Int64,1},1}, ref_seq::Array{Int64,1})

    V,T = size(N)
    hd = zeros(Float64,V,T)

    for v=1:V
        for t=1:T
            for i=1:length(N[v,t])
                hd[v,t] += hamm_dist(seq[i],ref_seq)*N[v,t][i]
            end
            hd[v,t] /= sum(N[v,t])
        end
    end

    return hd

end

function hamm_dist(seq1::Array{Int64,1}, seq2::Array{Int64,1})

    hd = 0.0

    @assert length(seq1) == length(seq2)
    for i=1:length(seq1)
        if seq1[i] != seq2[i]
            hd +=1
        end
    end

    return hd
end
