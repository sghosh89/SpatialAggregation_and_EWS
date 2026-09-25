library(here)
library(tidyverse)
library(future)
library(future.apply)
source(here("R/run_one_trajectory_on_ring_network_obsnoise.R"))

# ============================================================
# Fixed model parameters
# ============================================================

pars <- list(
  
  # Number of patches
  N = 100,
  
  # Dispersal strength
  D = 0.1,
  
  # Recovery rate
  lambda_start = -0.5,
  lambda_end   = -0.05,
  
  # Process noise
  sigma_local = 1,
  
  # gamma = sigma_shared / sigma_local = 0.5
  sigma_shared = 0.5,
  
  # Integration time step
  dt = 0.05,
  
  # Burn-in and transition periods
  burnin_time = 50,
  transition_time = 200,
  
  # Save every 5 integration steps
  sample_every = 5,
  
  # Rolling variance window
  window_fraction = 0.5
)

m_values <- c(1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

sigma_obs_values <- c(0, 0.25, 0.5, 1, 2)

n_rep <- 1000
seeds <- seq_len(n_rep)

# ============================================================
# Run one matched EWS-null replicate
#
# The SAME seed is used for:
#
#   - EWS trajectory
#   - stationary null trajectory
#
# Within the worker:
#
#   - ecological dynamics start from set.seed(seed)
#   - observation noise continues from that RNG stream
#   - set.seed(seed) is reset before random spatial sampling
#
# ============================================================

run_ews_null_obsnoise_once <- function(seed, sigma_obs, pars, m_values){
  
  # ----------------------------------------------------------
  # Non-stationary EWS trajectory
  out_ews <-  run_one_trajectory_on_ring_network_obsnoise(seed = seed,  pars = pars,  m_values = m_values,  sigma_obs = sigma_obs,  stationary = FALSE) %>%
    mutate(type = "EWS")
  
  # ----------------------------------------------------------
  # Matched stationary null trajectory
  
  out_null <- run_one_trajectory_on_ring_network_obsnoise(seed = seed,  pars = pars,  m_values = m_values,  sigma_obs = sigma_obs,  stationary = TRUE) %>%
    mutate(type = "Null")
  
  # ----------------------------------------------------------
  # Combine
  
  out<-bind_rows(out_ews,out_null) %>%
    mutate(sigma_ratio =  pars$sigma_shared / pars$sigma_local)
  
  return(out)
}

job_grid <- expand_grid(sigma_obs = sigma_obs_values, seed = seeds)

n_cores <- future::availableCores()
n_workers <- max(1, n_cores - 2)
n_workers

plan(multisession,  workers = n_workers)

# ============================================================
# PARALLEL RUN ACROSS all parameters
# ============================================================

system.time({
  
  results_list <- future_lapply(X = seq_len(nrow(job_grid)), 
                                FUN = function(kk){
                                  run_ews_null_obsnoise_once(seed = job_grid$seed[kk],
                                                             sigma_obs=job_grid$sigma_obs[kk],
                                                             pars = pars,
                                                             m_values = m_values)}, future.seed = TRUE)
})

results_all <- bind_rows(results_list)

saveRDS(results_all, here("Results/ring_EWS_tau_with_obsnoise_1000rep.RDS"))




