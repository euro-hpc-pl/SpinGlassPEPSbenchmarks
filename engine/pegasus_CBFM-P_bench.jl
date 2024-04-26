using LinearAlgebra
using MKL
using SpinGlassEngine
using SpinGlassNetworks
using SpinGlassTensors
#using SpinGlassExhaustive
using Logging
using CSV
using DataFrames
using Memoization
using JSON3

size = 1 #MPI.Comm_size(MPI.COMM_WORLD)
rank = 0 #MPI.Comm_rank(MPI.COMM_WORLD)

M, N, T = 7, 7, 3
INSTANCE_DIR = "$(@__DIR__)/../test/instances/pegasus_random/P8/CBFM-P/SpinGlass/single"
OUTPUT_DIR = "$(@__DIR__)/results/pegasus_random/P8/CBFM-P/droplets/final_bench"

if !Base.Filesystem.isdir(OUTPUT_DIR)
    Base.Filesystem.mkpath(OUTPUT_DIR)
end

BETAS = [0.5,]
LAYOUT = (GaugesEnergy,)
TRANSFORM = all_lattice_transformations

GAUGE =  NoUpdate
STRATEGY = Zipper
SPARSITY = Sparse
graduate_truncation = :graduate_truncate

INDβ = [3,]
MAX_STATES = [1024,]
BOND_DIM = [8,]

MAX_SWEEPS = [0,]
VAR_TOL = 1E-16
TOL_SVD = 1E-16
ITERS_SVD = 2
ITERS_VAR = 1
DTEMP_MULT = 2
METHOD = :psvd_sparse
I = [1,]
eng = [60, ]
hamming_dist = 74

disable_logging(LogLevel(1))
#BLAS.set_num_threads(1)

function pegasus_sim(inst, trans, β, Layout, bd, ms, eng, hamming_dist, mstates)
    δp = 0.0

    cl_h = clustered_hamiltonian(
        ising_graph(INSTANCE_DIR * "/" * inst),
        spectrum=full_spectrum,
        cluster_assignment_rule=pegasus_lattice((M, N, T))
        )

    params = MpsParameters(bd, VAR_TOL, ms, TOL_SVD, ITERS_SVD, ITERS_VAR, DTEMP_MULT, METHOD)
    search_params = SearchParameters(mstates, δp)

    net = PEPSNetwork{SquareCrossDoubleNode{Layout}, SPARSITY}(M, N, cl_h, trans)
    ctr = MpsContractor{STRATEGY, GAUGE}(net, [β/6, β/3, β/2, β], graduate_truncation, params)
    sol1, schmidts = low_energy_spectrum(ctr, search_params, merge_branches(ctr, :nofit, SingleLayerDroplets(eng, hamming_dist, :hamming)))

    sol2 = unpack_droplets(sol1, β)
    ig_states = decode_clustered_hamiltonian_state.(Ref(cl_h), sol2.states)

    ldrop = length(sol2.states)
    cRAM = round(Base.summarysize(Memoization.caches) * 1E-9; sigdigits=2)
    clear_memoize_cache()
    sol1, ctr, cRAM, schmidts, ldrop, sol2, ig_states
end

function run_bench(inst::String, β::Real, t, l, bd, ms, eng, hamming_dist, mstates, i)
    hash_name = hash(string(inst, β, t, l, bd, ms, eng, hamming_dist, mstates, i))
    out_path = string(OUTPUT_DIR, "/", hash_name, ".json")

    if isfile(out_path)
        println("Skipping for $β, $t, $l, $bd, $eng, $hamming_dist, $ms, $mstates.")
    else
        data = try
            tic_toc = @elapsed sol, ctr, cRAM, schmidts, ldrop, droplets, ig_states = pegasus_sim(inst, t, β, l, bd, ms, eng, hamming_dist, mstates)

            data = DataFrame(
                :instance => inst,
                :β => β,
                :Layout => l,
                :transform => t,
                :energy => sol.energies[begin],
                :probabilities => sol.probabilities,
                :discarded_probability => sol.largest_discarded_probability,
                :statistic => minimum(values(ctr.statistics)),
                :max_states => mstates,
                :bond_dim => bd,
                :max_sweeps => ms,
                :eng => eng,
                :hamming_dist => hamming_dist,
                :drop_eng => [droplets.energies],
                :drop_states => [droplets.states],
                :ig_states => [ig_states],
                :drop_prob => [droplets.probabilities],
                :drop_degeneracy => [droplets.degeneracy],
                :drop_ldp => [droplets.largest_discarded_probability],
                :drop_number => ldrop,
                :iters_svd => ITERS_SVD,
                :iters_var => ITERS_VAR,
                :dtemp_mult => DTEMP_MULT,
                :var_tol => VAR_TOL,
                :time => tic_toc,
                :cRAM => cRAM,
                :schmidts => schmidts
            )
        catch err
            data = DataFrame(
                :instance => inst,
                :β => β,
                :Layout => l,
                :transform => t,
                :max_states => mstates,
                :hamming_dist => hamming_dist,
                :iters_svd => ITERS_SVD,
                :iters_var => ITERS_VAR,
                :dtemp_mult => DTEMP_MULT,
                :bond_dim => bd,
                :max_sweeps => ms,
                :var_tol => VAR_TOL,
                :error => err
            )
        end

        json_data = JSON3.write(data)
        open(out_path, "w") do io
            print(io, json_data)
        end

    end
end

all_params = collect(
    Iterators.product(
        readdir(INSTANCE_DIR, join=false), BETAS, TRANSFORM, LAYOUT, BOND_DIM, MAX_SWEEPS, eng, hamming_dist, MAX_STATES, I)
)

for i ∈ (1+rank):size:length(all_params)
    run_bench(all_params[i]...)
    GC.gc()
end
