# src/GillespieCore.jl

module GillespieCore

"""
GillespieCore module provides core functions to simulate
ribosome dynamics on a lattice with periodic boundary conditions (PBC)
using a Gillespie stochastic simulation algorithm.

Functions:

- get_elongation_vector_pbc(lattice, rates, mobile, jammed, l_ribosome)
    Computes the elongation propensity vector indicating
    where ribosomes can move on the lattice. Also counts mobile and jammed particles.

- elongation_process_pbc(elongation_vector, internal_state_vec, posR, lattice, rates, J, J_w, J_c, k₊, mobile, jammed, paused, l_ribosome, num_clust)
    Executes one elongation event by moving a ribosome on the lattice
    and updating lattice state, internal states, currents, and elongation vector.

- change_internal_state_pbc(internal_state_vec, pßosR, lattice, rates, elongation_vector, k₋, k₊, mobile, jammed, paused, l_ribosome)
    Simulates switching of ribosomes between paused and active states,
    updating corresponding vectors and counters.

- create_lattice_l(L, ρ, l_ribosome)
    Creates a lattice configuration representing ribosomes of size l_ribosome
    randomly distributed according to density ρ on a lattice of length L.

- create_rates(L, γ)
    Creates a hopping rate array for each site on the lattice, defaulting to rate γ.

- positionRibosomes(lattice)
    Returns a vector marking positions of ribosomes (occupied sites) on the lattice.
"""

using LinearAlgebra
using Random
using Distributions
using StatsBase

export get_elongation_vector_pbc, elongation_process_pbc, change_internal_state_pbc, create_lattice_l, create_rates, positionRibosomes, get_internal_state_vec_pbc


function get_elongation_vector_pbc(lattice, rates, mobile, jammed, l_ribosome)

    elongation_vector = zeros(length(lattice)) # create propensity vector for lattice configuration

    for i in eachindex(lattice[1:end-l_ribosome]) # check all sites except last one
        
        if lattice[i] == 1 && sum(lattice[i:i+l_ribosome]) == 1 # mobile particle which means it contributes to possible transition in the elongation vector
            elongation_vector[i] += rates[i] # adding rate allows for weighted sampling
            mobile[1] += 1.0 # particle is mobile
        elseif lattice[i] == 1  && sum(lattice[i:i+l_ribosome]) == 2 # neigbor site is occupied
            jammed[1] += 1.0 # particle is jammed
        end
    end


    for i in length(lattice)-l_ribosome+1:length(lattice)
        distance_to_end = length(lattice)-i
        overlap = l_ribosome - distance_to_end
        if lattice[i] == 1 && sum(lattice[1:overlap]) == 0 # periodic boundary conditions
            elongation_vector[i] += rates[i]
            mobile[1] += 1.0
        elseif lattice[i] == 1 && sum(lattice[1:overlap]) == 1
            jammed[1] += 1.0
        end
    end

    return elongation_vector
end

function elongation_process_pbc(
    elongation_vector, internal_state_vec, posR, lattice, 
    rates, J, J_w, J_c, k₊, mobile, jammed, paused, l_ribosome, num_clust
    )
    
    moving_particle = sample(posR, Weights(elongation_vector)) # which particles moves, the position of the particle is weighted by the rates at the site
    
    if 1+l_ribosome <= moving_particle <= length(lattice)-l_ribosome # bulk hopping
        # particle jump --> change lattice
        
        lattice[moving_particle] = 0
        lattice[moving_particle+1] = 1

        # particle jump --> change position
        posR[moving_particle+1] = posR[moving_particle]+1 # change position first
        posR[moving_particle] = 0
        # only mobile particle can jump --> change internal state propensity
        internal_state_vec[moving_particle+1] = internal_state_vec[moving_particle]
        internal_state_vec[moving_particle] = 0

        # track site current
        # J[moving_particle] += 1
        # change elongation_vector

        elongation_vector[moving_particle] = 0 #particle always leaves current position
        #change_elong_vector_bulk_pbc(elongation_vector, moving_particle, lattice, rates, internal_state_vec,k₊, jammed ,mobile, l_ribosome, num_clust)
        change_elong_vector_bulk_pbc(elongation_vector, moving_particle, lattice, rates, internal_state_vec,k₊, jammed ,mobile, l_ribosome, num_clust)

    elseif moving_particle <= l_ribosome # initiation
        #hopping from site 1
        lattice[moving_particle] = 0
        lattice[moving_particle+1] = 1
        # changing position Vector
        posR[moving_particle+1] = posR[moving_particle]+1
        posR[moving_particle] = 0
        # update internal_state_vec
        internal_state_vec[moving_particle+1] = internal_state_vec[moving_particle]
        internal_state_vec[moving_particle] = 0
        # track site current
        # J[moving_particle] += 1
        # function to change elongation_vector
        elongation_vector[moving_particle] = 0 #particle always leaves current position
        change_elong_vector_init_pbc(elongation_vector, moving_particle, lattice, rates, internal_state_vec,k₊, mobile , jammed, l_ribosome, num_clust)
    
    elseif length(lattice)-l_ribosome+1 <= moving_particle <= length(lattice) #termination p.b.c

        change_elong_vector_term_pbc(elongation_vector, moving_particle, lattice, rates, internal_state_vec,k₊, mobile, jammed, l_ribosome, num_clust)
        
        if moving_particle == length(lattice)
            J[1] += 1
            
            if paused[1] == 0
                J_w[1] += 1
            else
                J_c[1] += 1
            end 
            
            lattice[moving_particle] = 0
            lattice[1] = 1

            posR[moving_particle] = 0
            posR[1] = 1

            internal_state_vec[1] = internal_state_vec[moving_particle]
            internal_state_vec[moving_particle] = 0

        else    
            # change lattice
            lattice[moving_particle] = 0
            lattice[moving_particle+1] = 1
            # change position
            posR[moving_particle+1] = posR[moving_particle] +1
            posR[moving_particle] = 0

            internal_state_vec[moving_particle+1] = internal_state_vec[moving_particle]
            internal_state_vec[moving_particle] = 0
        
        end
        # update elongation vector
        elongation_vector[moving_particle] = 0 # particle leaves position
        
    end

    return moving_particle
end

function change_elong_vector_bulk_pbc(
    elongation_vector, moving_particle, lattice, rates, 
    internal_state_vec, k₊, jammed, mobile, l_ribosome, num_clust
    )
    
    if 1+l_ribosome <= moving_particle < length(lattice)-l_ribosome

        # check site behind particles which moves
        if lattice[moving_particle-l_ribosome] == 1 && internal_state_vec[moving_particle-l_ribosome] == k₊# particle behind must also be active
            elongation_vector[moving_particle-l_ribosome] = rates[moving_particle-l_ribosome] # behind particle can move now
            # a jammed particle becomes mobile
            jammed[1] -= 1.0 
            mobile[1] += 1.0
            # num_clust[1] += 1
        end
        
        # check moving_particle 
        if sum(lattice[moving_particle+1:moving_particle+l_ribosome+1]) == 1 # particle is still mobile
            elongation_vector[moving_particle+1] = rates[moving_particle+1]
        else # we joined a cluster
            jammed[1] += 1
            mobile[1] -= 1
            # num_clust[1] -= 1
        end

        # cluster decoagulation
        if lattice[moving_particle+l_ribosome+1] == 0 && lattice[moving_particle-l_ribosome] == 1
            num_clust[1] += 1
        elseif lattice[moving_particle+l_ribosome+1] == 1 && lattice[moving_particle-l_ribosome] == 0
            num_clust[1] -= 1
        end
        # cluster coagulation
    elseif moving_particle == length(lattice)-l_ribosome

        if lattice[moving_particle-l_ribosome] == 1 && internal_state_vec[moving_particle-l_ribosome] == k₊# particle behind must also be active
            elongation_vector[moving_particle-l_ribosome] = rates[moving_particle-l_ribosome] # behind particle can move now
            # a jammed particle becomes mobile
            jammed[1] -= 1.0 
            mobile[1] += 1.0
            # num_clust[1] += 1
        end

        # check moving_particle 
        if lattice[1] == 0 # particle is still mobile
            elongation_vector[moving_particle+1] = rates[moving_particle+1]
        else # we joined a cluster
            jammed[1] += 1
            mobile[1] -= 1
            # num_clust[1] -= 1
        end

        if lattice[1] == 0 && lattice[moving_particle-l_ribosome] == 1
            num_clust[1] += 1
        elseif lattice[1] == 1 && lattice[moving_particle-l_ribosome] == 0
            num_clust[1] -= 1
        end

    end
    
end

function change_elong_vector_init_pbc(
    elongation_vector, moving_particle, lattice, 
    rates, internal_state_vec,k₊, mobile, 
    jammed, l_ribosome, num_clust
    )
    
    if sum(lattice[moving_particle:moving_particle+l_ribosome+1]) == 1 # particle is still mobile
        elongation_vector[moving_particle+1] = rates[moving_particle+1]
    else # we joined a cluster
        jammed[1] += 1
        mobile[1] -= 1
        # num_clust[1] -= 1
    end
    #check particle behind
    overlap = l_ribosome - moving_particle
    if lattice[end-overlap] == 1 && internal_state_vec[end-overlap] == k₊
        elongation_vector[end-overlap] = rates[end-overlap]
        jammed[1] -= 1.0
        mobile[1] += 1.0
        # num_clust[1] += 1
    end
    # decoagulation
    if lattice[end-overlap] == 1 && lattice[moving_particle+l_ribosome+1] == 0
        num_clust[1] += 1
    end

    if lattice[end-overlap] == 0 && lattice[moving_particle+l_ribosome+1] == 1
        num_clust[1] -= 1
    end

    # if lattice[end-overlap] == 1
    #     num_clust[1] += 1
    # end
    
end

function change_elong_vector_term_pbc(
    elongation_vector, moving_particle, lattice, rates, 
    internal_state_vec,k₊, mobile, jammed, l_ribosome, 
    num_clust
    )
    #check behind
   
    if lattice[moving_particle-l_ribosome] == 1 && internal_state_vec[moving_particle-l_ribosome] == k₊
        elongation_vector[moving_particle-l_ribosome] = rates[moving_particle-l_ribosome]
        mobile[1] += 1.0
        jammed[1] -= 1.0
        # num_clust[1] += 1
    end
    #l_ribosome = 1 -> check L[1]
    #l_ribosome = 2 -> 
    distance_to_end = length(lattice) - moving_particle
    overlap = l_ribosome - distance_to_end - 1
    # if moving_particle==length(lattice)
    #     overlap = l_ribosome
    # else
    #     overlap = l_ribosome - length(lattice[moving_particle+1:end]) 
    # end
    if l_ribosome == 1
        overlap = 0
    end

    if moving_particle == length(lattice) && sum(lattice[1:overlap+2]) == 0
        elongation_vector[1] = rates[1]
    elseif moving_particle == length(lattice) && sum(lattice[1:overlap+2]) == 1  
        jammed[1] += 1.0
        mobile[1] -= 1.0
        # num_clust[1] -= 1
    end

    
    if moving_particle != length(lattice)  && sum(lattice[1:overlap+2]) == 0
        elongation_vector[moving_particle+1] = rates[moving_particle+1]
    elseif moving_particle != length(lattice) && sum(lattice[1:overlap+2]) == 1
        jammed[1] += 1.0
        mobile[1] -= 1.0
        # num_clust[1] -= 1
    end
    # cluster decoagulation
    if sum(lattice[1:overlap+2]) == 0 && lattice[moving_particle-l_ribosome] == 1
        num_clust[1] += 1
    end
    # cluster coagulation
    if  sum(lattice[1:overlap+2]) == 1 && lattice[moving_particle-l_ribosome] == 0
        num_clust[1] -= 1
    end

    # if lattice[moving_particle-l_ribosome] == 1
    #     num_clust[1] += 1
    # end
end

function get_internal_state_vec_pbc(lattice, elongation_vector, k₋,k₊,mobile, jammed, paused) # where k₋ = waiting time paused state,k₊= rate entering active state

    internal_state_vec = zeros(length(lattice)) # initialize internal_state_vec 
    k₋_init = (k₋) / (k₋ + k₊) # makes sure we get close to steady state value for paused fraction

    for i in eachindex(lattice)

        if lattice[i] == 1 # only occupied sites can become paused
            RV = rand()

            if RV >= k₋_init # particle is in a paused state

                internal_state_vec[i] = k₋ # change propensity in internal state vec to rate at which particle is leaving the paused state
                elongation_vector[i] = 0 # can not elongate anymore
                paused[1] += 1.0 # add to paused number
                
                if i != length(lattice) && lattice[i+1] == 1 # jammed particle is chosen
                    jammed[1] -= 1.0 # reduce amount of jammed particles
                elseif i != length(lattice) && lattice[i+1] == 0  # mobile particle is chosen
                    mobile[1] -= 1.0 # reduce number of mobile particles
                end

                if i == length(lattice) && lattice[1] == 1
                    jammed[1] -= 1
                elseif i == length(lattice) && lattice[1] == 0
                    mobile[1] -= 1
                end

            else # particle is set to stay mobile state
                internal_state_vec[i] = k₊
            end

        end

    end
    return internal_state_vec
end


function change_internal_state_pbc(
    internal_state_vec, posR, lattice, rates, 
    elongation_vector, k₋,k₊, mobile, jammed, paused, l_ribosome)

    switching_pos = sample(posR, Weights(internal_state_vec))

    if internal_state_vec[switching_pos] == k₋ # switch from paused to active
        
        internal_state_vec[switching_pos] =k₊# propensity changes; now the particle can become paused with a rate f
        paused[1] -= 1 # substract from paused population

        if switching_pos <= length(lattice)-l_ribosome # all sites except last site
            
            if lattice[switching_pos+l_ribosome] == 0 # if site on the right is empty
                elongation_vector[switching_pos] = rates[switching_pos] # particle contributes to elongation
                mobile[1] += 1 # from paused population to mobile one
            elseif lattice[switching_pos+l_ribosome] == 1 # if right site is occupied
                jammed[1] += 1 # add to jammed population
            end

        elseif length(lattice)-l_ribosome+1 <= switching_pos <= length(lattice)-1 # from site L to site 1
            distance_to_end = length(lattice) - switching_pos
            overlap = l_ribosome - distance_to_end-1  
            if sum(lattice[1:1+overlap]) == 0 # empty
                elongation_vector[switching_pos] = rates[switching_pos] #last particle contributes
                mobile[1] += 1 # add to mobile
            elseif sum(lattice[1:1+overlap]) == 1 # first site occupied
                jammed[1] += 1 # add to jammed
            end

        elseif switching_pos == length(lattice)
            if sum(lattice[1:l_ribosome]) == 0 # empty
                elongation_vector[end] = rates[end] #last particle contributes
                mobile[1] += 1 # add to mobile
            elseif sum(lattice[1:l_ribosome]) == 1 # first site occupied
                jammed[1] += 1 # add to jammed
            end

         end
         

    else # means an non-paused particle was chosen which has to switch to a paused state
        internal_state_vec[switching_pos] = k₋ # non-paused particle can become paused
        elongation_vector[switching_pos] = 0 # paused particle does not contribute to elongation vector
        paused[1] += 1 # add to paused population
        
        if switching_pos <= length(lattice)-l_ribosome #

            if lattice[switching_pos+l_ribosome] == 0 # if particle choosen was mobile
                mobile[1] -= 1 # substract mobile particle
            elseif lattice[switching_pos+l_ribosome] == 1 # if particle choosen was jammed
                jammed[1] -= 1 # substract jammed particle
            end
        
        elseif length(lattice)-l_ribosome+1 <= switching_pos <= length(lattice)-1 # from site L to site 1
            distance_to_end = length(lattice) - switching_pos
            overlap = l_ribosome - distance_to_end-1  
            if sum(lattice[1:1+overlap]) == 0 # empty
                mobile[1] -= 1 # substract from mobile
            elseif sum(lattice[1:1+overlap]) == 1 # first site occupied
                jammed[1] -= 1 # substract from jammed
            end

        elseif switching_pos == length(lattice) # last site was choosen

            if sum(lattice[1:l_ribosome]) == 0 # empty
                mobile[1] -= 1 # substract from mobile
            elseif sum(lattice[1:l_ribosome]) == 1 # first site occupied
                jammed[1] -= 1 # substract from jammed
            end

        end
    end

    return switching_pos
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


# Creates an array of hopping rates γ for each site
function create_rates(L, γ)
    rates = ones(L) .* γ
    # rates[1] = α  # Uncomment if initiation rate α to be used
    # rates[end] = β # Uncomment if termination rate β to be used
    return rates
end


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


end # module
