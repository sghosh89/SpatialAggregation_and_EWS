# ============================================================
# ONE complete replicate with observation noise
#
# Purpose:
#   Simulate one N-patch ecological trajectory, add independent patch-level observation noise, aggregate the observed
#   trajectories, and calculate EWS Kendall's tau for all m.
#
# Observation model:
#   X_obs,i(t) = X_i(t) + eta_i(t)
#
# where:
#   eta_i(t) ~ N(0, sigma_obs^2)
#
# Input:
#   seed       = replicate seed
#   pars       = model parameter list
#   m_values   = aggregation scales
#   sigma_obs  = SD of observation noise
#   stationary = TRUE for stationary null, FALSE for EWS trajectory (non-stationary) 
#
# Output:
#   Data frame containing Kendall's tau for every m and  sampling strategy.
# ============================================================

source(here("R/simulate_ring_trajectory.R"))
source(here("R/get_ews_tau_all_m.R"))
source(here("R/aggregate_consecutive_ring.R"))
source(here("R/aggregate_random_ring.R"))

run_one_trajectory_on_ring_network_obsnoise <- function(seed, pars, m_values, sigma_obs, stationary = FALSE){
    
    
    # Simulate the TRUE N-patch ecological trajectory
    sim <- simulate_ring_trajectory(pars = pars, stationary = stationary, seed = seed)
    
    # ----------------------------------------------------------
    # Add independent patch-level observation noise
    # X_obs = X + eta
    # eta ~ N(0, sigma_obs^2)
    
    observation_error <- matrix(rnorm(length(sim$X), mean = 0, sd = sigma_obs), nrow = nrow(sim$X), ncol = ncol(sim$X))
    X_obs <- sim$X + observation_error
    
    # ----------------------------------------------------------
    #  Consecutive sampling for ALL m using OBSERVED trajectories
    Y_con <- aggregate_consecutive_ring(X = X_obs,  m_values = m_values)
    
    # ----------------------------------------------------------
    # Random sampling for ALL m using OBSERVED trajectories
    # Reset the SAME replicate seed specifically for spatial sampling
    set.seed(seed)
    Y_rand_obj <- aggregate_random_ring(X = X_obs, m_values = m_values)
    Y_rand <- Y_rand_obj$Y
    
    # ----------------------------------------------------------
    # EWS tau: consecutive
    tau_con <- get_ews_tau_all_m(Y_con, window_fraction = pars$window_fraction)
    tau_con$strategy <- "Consecutive"
    
    # ----------------------------------------------------------
    # EWS tau: random
    tau_rand <- get_ews_tau_all_m(Y_rand, window_fraction = pars$window_fraction)
    tau_rand$strategy <- "Random"
    
    # ----------------------------------------------------------
    #  Combine
    
    result <- bind_rows(tau_con, tau_rand)
    
    result$seed <- seed
    result$stationary <- stationary
    result$D <- pars$D
    result$m_over_N <- result$m / pars$N
    result$sigma_obs <- sigma_obs
    
    return(result)
}