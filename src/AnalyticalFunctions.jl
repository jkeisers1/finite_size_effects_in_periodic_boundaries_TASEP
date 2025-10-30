# === Core analytical functions for the TASEP with pausing dynamics ===
# Parameters:
# α = initiation rate
# β = termination rate
# γ = hopping rate
# k₊ = rate of entering paused state
# k₋ = rate of leaving paused state
# L = system length (e.g., number of sites)
# l_ribosome = size of extended particle (ribosome length)

# Pausing fraction
Φ(k₋, k₊) = k₊ / (k₋ + k₊)

# Max density (closed form)
function ρ_max_fct(k₋, k₊, γ)
    x = (1 + k₊ / k₋) / (γ * Φ(k₋, k₊) / k₋)
    return -x + sqrt(x^2 + x)
end

# Bulk densities for boundary-driven phases
ρ_HD(β, γ) = 1 - β / γ
ρ_LD(α, k₋, k₊, γ) = (α / γ) * (1 + k₊ / k₋)

# Wang current
Wang(ρ, k₋, k₊, γ) = γ * ρ * (1 - ρ) * k₋ / (k₋ + k₊ + γ * ρ * Φ(k₋, k₊))

# Boundary currents via Wang
Wang_HD(β, k₋, k₊, γ) = Wang(ρ_HD(β, γ), k₋, k₊, γ)
Wang_LD(α, k₋, k₊, γ) = Wang(ρ_LD(α, k₋, k₊, γ), k₋, k₊, γ)

# Critical rates
α_crit_fct(k₋, k₊, γ) = γ * ρ_max_fct(k₋, k₊, γ) / (1 + k₊ / k₋)
β_crit_fct(k₋, k₊, γ) = (1 - ρ_max_fct(k₋, k₊, γ)) * γ
# (LD-only heuristic from your notes)
k₊_crit_fct(k₋, α, γ) = (-2α + γ) / ((α / k₋) * (α / k₋ + 2))

# Particle fractions
active_fct(ρ, k₋, k₊) = ρ / (1 + k₊ / k₋)
paused_fct(ρ, k₋, k₊) = ρ * (k₊ / k₋) / (1 + k₊ / k₋)

# Extended-particle currents
J_ext(γ, ρ, l_ribosome) = (γ * ρ * (1 - ρ * l_ribosome)) / (1 - ρ * (l_ribosome - 1))
J_ext_k₊(γ, l_ribosome, k₋, k₊, ρ) = (k₋ * γ * ρ * (1 - ρ * l_ribosome)) /
                                     ((1 - ρ * (l_ribosome - 1)) * (k₋ + k₊ + γ * ρ * Φ(k₋, k₊)))

# === "Woosh" mixture (cluster vs. free) ===
pw_fct(ρ, L, k₋, k₊) = (1 - Φ(k₋, k₊))^(L * ρ)  # cluster probability surrogate
J_whoosh(γ, ρ, L, k₋, k₊) = γ * ρ * (1 - ρ) * pw_fct(ρ, L, k₋, k₊) +
                            Wang(ρ, k₋, k₊, γ) * (1 - pw_fct(ρ, L, k₋, k₊))
Woosh_fraction(γ, ρ, L, k₋, k₊) = (γ * ρ * (1 - ρ) * pw_fct(ρ, L, k₋, k₊)) / J_whoosh(γ, ρ, L, k₋, k₊)
Cluster_Wang_fraction(γ, ρ, L, k₋, k₊) = (Wang(ρ, k₋, k₊, γ) * (1 - pw_fct(ρ, L, k₋, k₊))) / J_whoosh(γ, ρ, L, k₋, k₊)

# Particle classes (decomposition using Wang)
ψ(ρ, k₋, k₊) = ρ * (k₊) / (k₋ + k₊)                 # paused fraction
μ(ρ, k₋, k₊, γ) = Wang(ρ, k₋, k₊, γ) / γ            # mobile fraction (flow/γ)
ν(ρ, k₋, k₊, γ) = ρ - μ(ρ, k₋, k₊, γ) - ψ(ρ, k₋, k₊) # jammed/clustered fraction

# Alternative mean-field current (Greulich-style)
greulich(ρ, k₋, k₊, γ) = γ * ρ * (1 - ρ) / (1 + (γ / k₋) * ρ * (k₊ / (k₊ + k₋)))

# === Transient current pieces ===
function av_current_half(ρ, L, γ, t)
    if t < ρ * L
        return γ / L * (t / 6)
    elseif t < L / (4 * ρ)
        return (γ * ρ^2 * L / 6 + (γ * ρ / 3) * (ρ * L + 3t - 4sqrt(ρ * L * t))) / t
    else
        return (ρ^2 * L / 6 + (ρ / 3) * (ρ * L + 3L / (4ρ) - 2L) + ρ * (1 - ρ) * t -
                L / 4 * (1 - ρ) - (L^2 / 48) * (4ρ / L - 1 / t)) * (γ / t)
    end
end

av_trans_cur(ρ, L, γ, t) = ρ <= 0.5 ? av_current_half(ρ, L, γ, t) : av_current_half(1 - ρ, L, γ, t)

# Full transient current and its split (woosh vs. cluster) with kp≡k₊, km≡k₋
function J_transient(ρ, kp, km, L)
    fa = 1 - kp / (km + kp)
    D = 1 / (1 - fa)
    P₀ = fa^(ρ * L)
    return ρ == 0 ? 0 :
           P₀ * av_trans_cur(ρ, L, 1, 1 / (ρ * L * kp)) +
           (1 - P₀) * km * (1 - ρ) * ρ * L * D / (ρ * L + D - 1)
end

J_transient_w(ρ, kp, km, L) = ρ == 0 ? 0 : (fa = 1 - kp / (km + kp); P₀ = fa^(ρ * L);
                                            P₀ * av_trans_cur(ρ, L, 1, 1 / (ρ * L * kp)))
J_transient_c(ρ, kp, km, L) = ρ == 0 ? 0 : (fa = 1 - kp / (km + kp); D = 1 / (1 - fa); P₀ = fa^(ρ * L);
                                            (1 - P₀) * km * (1 - ρ) * ρ * L * D / (ρ * L + D - 1))

# === Error / comparison utilities ===
absolute_error(numerical, analytical) = abs.(numerical .- analytical)

function relative_max_error(numerical, analytical)
    max_num = maximum(numerical)
    max_ana = maximum(analytical)
    return abs(max_num - max_ana) / abs(max_num)
end

function density_at_max(numerical_current, density_grid)
    idx = argmax(numerical_current)
    return density_grid[idx]
end

function relative_max_error_rel_ana(numerical, analytical)
    max_ana = maximum(analytical)
    return abs(maximum(numerical) - max_ana) / abs(max_ana)
end

function positional_peak_error(J_num, J_ana, ρ_grid; normalised=true)
    @assert length(J_num) == length(J_ana) == length(ρ_grid)
    Δρ = abs(ρ_grid[argmax(J_num)] - ρ_grid[argmax(J_ana)])
    return normalised ? Δρ / (maximum(ρ_grid) - minimum(ρ_grid)) : Δρ
end

absolute_bounded_peak_error(J_num, J_ana) = abs(maximum(J_num) - maximum(J_ana)) /
                                             (abs(maximum(J_num)) + abs(maximum(J_ana)))
absolute_peak_error(J_num, J_ana) = abs(maximum(J_num) - maximum(J_ana)) / abs(maximum(J_ana))

rmse(J_num, J_ana) = sqrt(mean((J_num .- J_ana).^2))
relative_error(numerical, analytical) = abs.(numerical .- analytical) ./ abs.(analytical)

mean_squared_error(numerical, analytical) = mean((numerical .- analytical).^2)

function norm_mean_squared_error(numerical, analytical)
    mse = mean((numerical .- analytical).^2)
    range_ana = maximum(analytical) - minimum(analytical)
    return mse / (range_ana^2)
end

squared_difference(numerical, analytical) = (numerical .- analytical).^2
