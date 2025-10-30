using Random
using Distributions
using Distributed
# using BenchmarkTools

# -- Utility: Create logarithmic ranges --
function create_log10_range_list(starting_point, stopping_point, step_expo)
    list = []
    tmp = collect(range(starting_point, stop=stopping_point, step=step_expo))
    for item in tmp
        push!(list, 10^(item))
    end
    return list
end

# -- Lattice Creation Functions --

# Creates lattice with density ρ and random distribution of particles
function create_lattice(L, ρ)
    lattice = zeros(L)
    while mean(lattice) < ρ
        index = findfirst(x -> x == 0, lattice)
        lattice[index] += 1
    end
    shuffle!(lattice)
    return lattice
end

# Creates lattice with density ρ but particles placed consecutively without shuffling
function create_lattice_stuck(L, ρ)
    lattice = zeros(L)
    while mean(lattice) < ρ
        index = findfirst(x -> x == 0, lattice)
        lattice[index] += 1
    end
    return lattice
end

# Creates lattice for extended particles (ribosomes) of length l_ribosome, with random placement if l_ribosome=1
function create_lattice_l(L, ρ, l_ribosome)
    lattice = zeros(L)
    count = 1
    for i in 1:floor(Int64, L * ρ / l_ribosome)
        lattice[count] = 1
        count += l_ribosome
    end
    if l_ribosome == 1
        shuffle!(lattice)
    end
    return lattice
end

# Creates lattice for extended particles (ribosomes) without shuffling
function create_lattice_l_stuck(L, ρ, l_ribosome)
    lattice = zeros(L)
    count = 1
    for i in 1:floor(Int64, L * ρ / l_ribosome)
        lattice[count] = 1
        count += l_ribosome
    end
    return lattice
end

# -- Rate Initialization --

# Creates an array of hopping rates γ for each site
function create_rates(L, γ)
    rates = ones(L) .* γ
    # rates[1] = α  # Uncomment if initiation rate α to be used
    # rates[end] = β # Uncomment if termination rate β to be used
    return rates
end

# -- Position Extraction --

# Returns a vector of the same length as lattice, with index at positions where lattice == 1, else 0
function positionRibosomes(lattice)::Vector{Int64}
    position = zeros(length(lattice))
    i = 1
    for site in lattice
        if site == 1
            position[i] = i
        end
        i += 1
    end
    return position
end
