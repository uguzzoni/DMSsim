module DMSsim

	using Random,Distributions, DataStructures

    export dms_display!, dms_display, dms_reads,  dms_muts!

    include("simulate.jl")
    include("RandSeqsLib.jl")
	include("mutations.jl")

end # module
