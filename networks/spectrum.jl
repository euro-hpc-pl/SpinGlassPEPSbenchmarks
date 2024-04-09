
using SpinGlassNetworks

function bench(instance::String, size::NTuple{3,Int}, max_states::Int = 100)
    ig = ising_graph(instance)
    cl = split_into_clusters(ig, super_square_lattice(size))
    @time sp = brute_force(cl[1, 1], num_states = max_states)
    nothing
end

println("Threads: ", Threads.nthreads())
bench("$(@__DIR__)/pegasus_droplets/2_2_3_00.txt", (2, 2, 24))
