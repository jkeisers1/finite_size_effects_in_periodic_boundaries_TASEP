cd(@__DIR__)  # This sets current directory to where main.jl is located

include("src/GillespieCore.jl") 
include("src/GillespieSimulation.jl")  # Assuming this is the file where Gillespie_pbc is defined
include("src/Transient_Basics.jl")
include("src/Clustering.jl")

using .GillespieCore
using .GillespieSimulations  # Assuming this is the module/file where Gillespie_pbc is defined
using .GillespieClusterAnalysis

function main()
    # Define parameters for your simulation
    L = 500                      # lattice length
    ρ = 0.01                      # density of particles
    γ = 1.0                      # hopping rate
    k₋ = 10^-5                   # unpausing rate
    k₊ = 10^-4                    # pausing rate
    number_of_measurements = 100
    steps_to_steady_state = 100
    number_of_steps = 5_000_000
    l_ribosome = 1              # ribosome length
    kymo = false                 # whether to generate kymograph or not
    
    # Run the Gillespie simulation
    current_matrix, density_matrix, time_vec = GillespieSimulations.Gillespie_pbc(
        L, ρ, γ, k₋, k₊, number_of_measurements, steps_to_steady_state, number_of_steps, l_ribosome
    )
    
    # Simple output to confirm run completed
    println("Simulation complete.")
    println("Density matrix size: ", size(density_matrix))
    println("Current matrix size: ", size(current_matrix))
    println("Time vector length: ", length(time_vec))
    
    # Optionally return results for further processing
    return current_matrix, density_matrix, time_vec
end

# Run main if the script is called directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end


J, dens, time = main()

J