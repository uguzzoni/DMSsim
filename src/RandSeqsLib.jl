
#######################################

function modify_seq(s,i,a,A,L)
    ss=convert(NTuple{L,Int}, s)
    s2=collect(ss)
    s2[i]=a
    #return convert(Sequence{A,L},tuple(s2...))
    return tuple(s2...)
end

####################################

#constructing library of one mutations from a sequence
function one_mutation_library(wt,A,L)

    mut_sequences=[wt]
    for i=1:L
        s=deepcopy(wt)
        for a=1:A
            if a!= wt[i]
                s=modify_seq(wt,i,a,A,L)
                push!(mut_sequences,s)
            end
        end
    end
    return mut_sequences;
end


#constructing library of one and two points mutations from a sequence
function one_two_mutations_library(wt,A,L)

    mut_sequences=[wt]
    for i=1:L
        s=deepcopy(wt)
        for a=1:A
            if a!= wt[i]
                s=modify_seq(wt,i,a,A,L)
                push!(mut_sequences,s)
            end
            for j=i+1:L
                for b=1:A
                    if b!= wt[j]
                        s2=modify_seq(s,j,b,A,L)
                        push!(mut_sequences,s2)
                    end
                end
            end
        end
    end
    return mut_sequences;
end

#constructing library sampled approx uniformly at distance 3:L
function distance_mutations_library(wt,A,L;NS=30000,d_inf=3,d_sup=L)

    #rand_sequences=[rand(Sequence{A,L})]
    rand_sequences = [(rand(1:A, L)...,)]
    #mutate d times wt
    for d=d_inf:d_sup
        for n=1:NS
            s=deepcopy(wt)
            rpos=randperm(L)
            for k=1:d
                a=rand(1:A)
                s=modify_seq(s,rpos[k],a,A,L)
            end
            push!(rand_sequences,s)
        end
    end
    return rand_sequences
end
