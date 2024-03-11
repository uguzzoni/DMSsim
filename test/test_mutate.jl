
using Distributions
using Random

include("../src/mutations.jl")

A = 20
L = 10
M = 5000
p = 0.01

N=ones(Int64,M)

for i=1:M
    N[i]=1000
end

seq=Array{Int64,1}[]
for i=1:M
       push!(seq,rand(1:A,L))
end

seq_old=deepcopy(seq)
seq1=deepcopy(seq)
seq2=deepcopy(seq)

@time new_N = mutate!(seq1, N, L, p; A=A);
@time new_N_old = mutate_old!(seq_old, N, L, p*L; A=A);
@time new_N2 = mutate_cum!(seq2, N, L, p; A=A);
