# src/GillespieSimulations.jl

module GillespieSimulations

using ..GillespieCore
using ..GillespieClusterAnalysis   # relative import of Clustering module (one level up)
using Random
using Distributions

export Gillespie_pbc, standard_pbc


function Gillespie_pbc(
    L, ρ, γ, k₋,k₊, number_of_measurements, steps_to_steady_state, number_of_steps, l_ribosome, kymo=false
    )

    #number_of_measurements = Int64(round((run_time-starting_t)/delta_t))

    # initialize particle classes
    mobile = [0.0]
    jammed = [0.0]
    paused = [0.0]

    # initialize vectors
    ρ_vec = Vector{Float64}(undef,number_of_measurements) # density vector
    ρ_weighted = [0.0] # density after weighted by time step the configuration

    mobile_vec = Vector{Float64}(undef, number_of_measurements) # mobile vector
    mobile_weighted = [0.0] # mobile fraction of particles

    jammed_vec= Vector{Float64}(undef,number_of_measurements)# jammed vector
    jammed_weighted = [0.0]# jammed fraction of particles

    paused_vec = Vector{Float64}(undef,number_of_measurements)# paused vector
    paused_weighted = [0.0]# paused fraction of particles

    cluster_number_dist = initiate_cluster_buckets(ρ, L, l_ribosome)
    cluster_length_dist = initiate_cluster_buckets(ρ, L, l_ribosome)

    current = Vector{Float64}(undef,number_of_measurements) # current vector
    
    J = [0] # current counter at last site
    J_w = [0] #current when there are no paused particles
    J_c = [0] #current when there is at least 1 paused particle
    
    t_c_vec = Vector{Float64}(undef, number_of_measurements)
    t_w_vec = Vector{Float64}(undef, number_of_measurements)  
    t_vec = Vector{Float64}(undef, number_of_measurements)
    
    J_c_vec = Vector{Float64}(undef, number_of_measurements)
    J_w_vec = Vector{Float64}(undef, number_of_measurements)
    # rates can not be the same
    if k₋ == k₊
        k₊*= 0.99999
    end

    lattice = create_lattice_l(L, ρ, l_ribosome)
    rates = create_rates(L, γ)
    posR = positionRibosomes(lattice)
    elongation_vector= get_elongation_vector_pbc(lattice, rates, mobile, jammed,l_ribosome)
    internal_state_vec = get_internal_state_vec_pbc(lattice, elongation_vector, k₋, k₊, mobile, jammed, paused)
    
    m_start = sum(elongation_vector .== γ)
    p_start = sum(internal_state_vec .== k₋)
    j_start = sum(lattice) - m_start - p_start
    mobile, paused, jammed = [m_start], [p_start], [j_start]

    if l_ribosome == 1
        num_clust = [length(cluster(lattice))]
    else
        num_clust = [1.0]
    end

    t_c = 0.0
    t_w = 0.0

    t = 0.0 # setup time counter
    j = 0

    #while t < starting_t # if time is below starting_t, no measurements are taken
    while j <= steps_to_steady_state 
        elong_sum = sum(elongation_vector) # sum of all the rates in the elongation vector
        state_sum = sum(internal_state_vec) # sum of all the rates in the internal state vector
        total_sum = elong_sum + state_sum 
        RV = rand() # Random variable to choose between elongation and internal state switching

        if RV <= elong_sum/total_sum # elongation is chosen
            elongation_process_pbc(
                elongation_vector, internal_state_vec, posR, lattice, rates, J, J_w, J_c, k₊, mobile, jammed, paused, l_ribosome, num_clust
                )
        else # switching happens
            change_internal_state_pbc(
                internal_state_vec, posR, lattice, rates, elongation_vector, k₋,k₊, mobile, jammed, paused, l_ribosome
                )
        end
        t += rand(Exponential(1/total_sum)) # add time depending on the lattice configuration
        j += 1
    end
    
    J[1] = 0 # reset current
    J_c[1] = 0
    J_w[1] = 0
    #t_vec = []
    if kymo == true
        kymo_matr = lattice
        intern_matr = internal_state_vec
        t_vec = [0.0]
    end

    for i in 1:number_of_measurements
        j = 0
        t = 0.0 # reset time
    
        if kymo == true
            kymo_matr = lattice
            intern_matr = internal_state_vec
            t_vec = [0.0]
            first_t = 0
            t_step = 1/minimum([k₊,k₋])
        end

        #while t <= delta_t # if time is below delta_t (time interval)
        while j <= number_of_steps
            elong_sum = sum(elongation_vector)
            state_sum = sum(internal_state_vec)
            total_sum = elong_sum + state_sum
            RV = rand()
            time_added = rand(Exponential(1/total_sum)) # observeables are weigthed by their configuration
            t += time_added # add to time

            if paused[1] == 0
                t_w += time_added
            else
                t_c += time_added
            end

            mobile_weighted[1] += (mobile[1]/L) * time_added # density of mobile particles weighted by time from previous configuration
            paused_weighted[1] += (paused[1]/L) * time_added # density of paused particles weighted by time from previous configuration
            jammed_weighted[1] += (jammed[1]/L) * time_added # density of jammed particles weighted by time from previous configuration 
            
            #find_all_cluster_lengths(lattice, cluster_length_dist, cluster_number_dist, time_added, num_clust)
            #normalization_cluster_length += num_clust[1] time_added
            #print("\ntest_num: $tmp1\n num_clust:$num_clust")
    
            if kymo == true
                if i == number_of_measurements # for the last measurement use kymo
                    if t >= first_t + t_step
                        push!(t_vec, t)
                        kymo_matr = hcat(kymo_matr, lattice)
                        intern_matr = hcat(intern_matr, internal_state_vec)
                        first_t += t_step
                    end
                end
            end

            if RV <= elong_sum/total_sum # evolve system
                elongation_process_pbc(
                elongation_vector, internal_state_vec, posR, lattice, rates, J, J_w, J_c, k₊, mobile, jammed, paused, l_ribosome, num_clust
                )
            else
                change_internal_state_pbc(
                    internal_state_vec, posR, lattice, rates, elongation_vector, k₋,k₊, mobile, jammed, paused, l_ribosome
                    )
            end

            j += 1 
        
        end
        
        # divide the measurements by delta_t/batch size
        current[i]= J[1]/t
        J[1] = 0
        J_c_vec[i] = J_c[1]/t
        J_c[1] = 0
        J_w_vec[i] = J_w[1]/t
        J_w[1]=0

        t_vec[i] = t
        t_w_vec[i] = t_w 
        t_w = 0
        t_c_vec[i] = t_c 
        t_c = 0

        mobile_vec[i] = mobile_weighted[1]/t
        paused_vec[i] = paused_weighted[1]/t
        jammed_vec[i] = jammed_weighted[1]/t

        mobile_weighted[1] = 0.0    # reset mobile counter
        paused_weighted[1] = 0.0    # reset paused counter
        jammed_weighted[1] = 0.0    # reset jammed counter
        #push!(t_vec, t)
    end
    # J_mean::Vector{Float64} = current
    J_mean = vcat(mean(current), mean(J_w_vec), mean(J_c_vec)) 
    # ρ_matrix= hcat(mobile_vec, paused_vec,jammed_vec)
    ρ_matrix = vcat(mean(mobile_vec), mean(paused_vec), mean(jammed_vec), ρ)
    # cluster_matrix = hcat(h_cluster_number, p_cluster_number)
    #(number_of_measurements*mean(t_vec))
    cluster_number_dist = Dict(key => value / 1  for (key,value) in cluster_number_dist)
    cluster_length_dist = Dict(key => value / 1 for (key,value) in cluster_length_dist)
    cluster_dict = Dict{String, Dict{Float64,Float64}}("NumberOfClusters" => cluster_number_dist, "ClusterLengthDistribution" =>cluster_length_dist) 
    
    time = vcat(mean(t_vec), mean(t_w_vec), mean(t_c_vec))

    stdev::Tuple{Float64, Float64} = (std(current) / mean(current), std(t_vec)/mean(t_vec))
    
    # if kymo == true
    #     return J_mean, ρ_matrix, t_vec, kymo_matr, intern_matr
    # end

    return J_mean, ρ_matrix, cluster_dict, stdev, time
end



function standard_pbc(L, ρ_list, init, elong, term, run_time)
    tmp = []
    J= [0]
    for item in ρ_list
        t = 0
        J[1] = 0
        lattice = create_lattice(L, item)
        rates = create_rates(L, elong, α = init, β = term)
        posR = positionRibosomes(lattice)
        elongation_vector= get_elongation_vector_pbc(lattice, rates, mobile, jammed)
        internal_state_vec = get_internal_state_vec_pbc(lattice, elongation_vector, k₋,k₊, mobile, jammed, paused)

        while t <= run_time
            elongation_process_pbc(elongation_vector, internal_state_vec, posR, lattice, rates, J,k₊, mobile, jammed)
            t += rand(Exponential(1/sum(elongation_vector)))
        end
        push!(tmp, J[1]/t)
    end
    return tmp
end



function kymo_Gillespie(L, ρ_list, init, term, elong, k₋, k₊_list, run_time, starting_t, delta_t)
    data_matrix = Array{Any}(undef, length(ρ_list), length(k₊_list))
    i= 1
    for ρ in ρ_list
        j = 1
        for k₊ in k₊_list

            if k₊ == k₋
                k₊ *=0.999999
            end
            tmp_J,tmp_ρ, tmp_t, tmp_kymo, tmp_intern =  Gillespie_pbc(L, ρ, init, term, elong, k₋, k₊, run_time, starting_t, delta_t;kymo=true)
            data = [tmp_J, tmp_ρ, tmp_t, tmp_kymo, tmp_intern]
            data_matrix[i,j] = data
            j += 1
        end

        i += 1
    end

    return data_matrix
end
end # module
