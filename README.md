# Simulation of a Deep Mutational Scanning experiment 
# Julia Package


Simulates a DMS experiment, where N and n are S×V×T and S×V×(T-1) tensors, where S is the number of variants in the experiment, V the replicates, T the rounds, N[s,v,t] is the number of variants s in the population in the experiment replicate v at round t, and n[s,v,t] the number of selected variants. Finally p[s] is the selectivity of sequence s.

For each replicate uses the initial condition N[:,v,1]. For each round and replicate sets Ntot[v,t] from sum(N[:,v,t]). Overwrites N[:,:,2:T] and n[:,:,:] with the simulation results.


## Installation

This code is written in Julia language. To add `DMSsim` package on Julia REPL run
```
Import("Pkg");Pkg.add("git@github.com:uguzzoni/DMSsim.git")
```

## Quick start

```julia
 
            using Random,DMSsim

            V = 2; T = 3;
            S = 30
            p = 0.1rand(S)

            N0 = rand(10^10 : 10^13)
            N,n = dms_display(V, T, N0, p)

```

## License
[MIT license](LICENSE)
