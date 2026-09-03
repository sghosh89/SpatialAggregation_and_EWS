# ============================================================
# ONE complete replicate
#
# THIS WILL LATER BE THE PARALLEL WORKER
# ============================================================
source(here("R/simulate_ring_trajectory.R"))
source(here("R/get_ews_tau_all_m.R"))
source(here("R/aggregate_consecutive_ring.R"))
source(here("R/aggregate_random_ring.R"))

run_one_trajectory_on_ring_network <- function(seed, pars, m_values, stationary = FALSE){
  
  # ----------------------------------------------------------
  # 1. Simulate one N-patch ecological trajectory
  
  sim <- simulate_ring_trajectory(pars = pars, stationary = stationary, seed = seed)
  
  # ----------------------------------------------------------
  # 2. Consecutive sampling for ALL m
  
  Y_con <- aggregate_consecutive_ring(X = sim$X, m_values = m_values)
  
  # ----------------------------------------------------------
  # 3. Reset replicate seed before random spatial sampling
  #
  # This ensures that a given replicate seed always produces
  # the same nested random patch ordering, independently of
  # the random numbers consumed during ecological simulation.
  
  set.seed(seed)
  Y_rand_obj <- aggregate_random_ring(X = sim$X,  m_values = m_values)
  Y_rand <- Y_rand_obj$Y
  
  # ----------------------------------------------------------
  # 4. EWS tau: consecutive
  
  tau_con <- get_ews_tau_all_m(Y_con, window_fraction = pars$window_fraction)
  tau_con$strategy <- "Consecutive"
  
  # ----------------------------------------------------------
  # 5. EWS tau: random
 
  tau_rand <- get_ews_tau_all_m(Y_rand, window_fraction = pars$window_fraction)
  tau_rand$strategy <- "Random"
  
  # ----------------------------------------------------------
  # 6. Combine
  
  result <- bind_rows(tau_con,  tau_rand)
  result$seed <- seed
  result$stationary <- stationary
  result$D <- pars$D
  result$m_over_N <-  result$m / pars$N
  
  return(result)
}