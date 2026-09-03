# ============================================================
# Run one ring trajectory with spatially correlated observation error
#
# PURPOSE
# -------
# Simulate the ecological ring system, add spatially correlated observation error at the PATCH level,
# aggregate patches, and calculate Kendall's tau of rolling variance for all aggregation scales.
#
# Observation model:
#
# X_obs_i(t) = X_i(t) + eta_i(t)
#
# eta_i(t) = sigma_obs * [sqrt(rho_obs) * Z_common(t) + sqrt(1-rho_obs) * Z_i(t)]
#
# where: Z_common(t) ~ N(0,1) shared among all patches at time t
#        Z_i(t) ~ N(0,1) independent among patches
#
# Therefore: Var(eta_i) = sigma_obs^2, Corr(eta_i, eta_j) = rho_obs
#
# for i != j.
# Observation errors are independent THROUGH TIME.
# Only their spatial correlation is varied.
#
# ============================================================

library(tidyverse)
library(here)

source(here("R/simulate_ring_trajectory.R"))
source(here("R/get_ews_tau_all_m.R"))
source(here("R/aggregate_consecutive_ring.R"))
source(here("R/aggregate_random_ring.R"))


# ============================================================
# Worker
# ============================================================

run_one_trajectory_on_ring_network_correlated_obsnoise <-  function(seed, pars, m_values, sigma_obs, rho_obs, stationary = FALSE){
    
    # -------------------------------
    # Check rho_obs
    if(rho_obs < 0 || rho_obs > 1){
      stop("rho_obs must lie between 0 and 1.")
    }
    
    
    #------------------------------------------------
    # 1. Simulate TRUE ecological trajectory
    # simulate_ring_trajectory() internally calls:
    # set.seed(seed)
    
    sim <- simulate_ring_trajectory(pars = pars, stationary = stationary, seed = seed)
    
    #---------------------------------------------------------
    # 2. Generate spatially correlated observation error
    #
    # IMPORTANT:
    # We do NOT call set.seed() again here.
    #
    # Therefore observation-error draws continue from the same RNG stream initialized by 'seed' inside the
    # ecological simulation.
    #
    # Because EWS and Null simulations use the same seed and consume the same number of random numbers during
    # the ecological simulation, their observation-error innovations are matched.
    
    n_time  <- nrow(sim$X) # time points
    N       <- ncol(sim$X) #patch number
    # --------------------------------------------------------
    # Common observation-error component    #
    # One random number at each time point, shared across all patches.
    
    Z_common <- rnorm(n_time, mean = 0, sd = 1)
    Z_common_matrix <- matrix(Z_common, nrow = n_time, ncol = N, byrow=F)# each column gets same sequence of random numbers
    # --------------------------------------------------------
    # Independent patch-specific component
    Z_local <- matrix( rnorm(n_time * N, mean = 0, sd = 1), nrow = n_time, ncol = N, byrow=F)
    
    # --------------------------------------------------------
    # Construct correlated observation error
    observation_error <-  sigma_obs * (sqrt(rho_obs) * Z_common_matrix + sqrt(1 - rho_obs) * Z_local)
    
    # --------------------------------------------------------
    # Observed ecological state
    X_obs <- sim$X + observation_error
    
    # --------------------------------------------------------
    # 3. Consecutive aggregation
    Y_con <- aggregate_consecutive_ring(X = X_obs, m_values = m_values)
    
    #----------------------------------------------------------
    # 4. Random aggregation
    #
    # Reset SAME replicate seed before choosing random spatial locations.
    set.seed(seed)
    
    Y_rand_obj <- aggregate_random_ring(X = X_obs, m_values = m_values)
    Y_rand <- Y_rand_obj$Y
    
    # ----------------------------------------------------------
    # 5. Calculate EWS Kendall tau
    tau_con <- get_ews_tau_all_m(Y_con, window_fraction = pars$window_fraction)
    tau_con$strategy <- "Consecutive"
    
    tau_rand <- get_ews_tau_all_m(Y_rand, window_fraction = pars$window_fraction)
    tau_rand$strategy <- "Random"
    
    #----------------------------------------------------------
    # 6. Combine
    result <- bind_rows(tau_con,  tau_rand)
    
    #-----------------------------------------------------
    # 7. Add metadata
    
    result$seed       <- seed
    result$stationary <- stationary
    result$D          <- pars$D
    result$m_over_N   <- result$m / pars$N
    result$sigma_obs  <- sigma_obs
    result$rho_obs    <- rho_obs
    
    return(result)
  }