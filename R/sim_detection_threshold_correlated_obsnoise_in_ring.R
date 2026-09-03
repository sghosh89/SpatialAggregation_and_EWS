library(here)
library(tidyverse)
library(future)
library(future.apply)

source(here("R/run_one_trajectory_on_ring_network_correlated_obsnoise.R"))

# -----------------------------------------------
# Fixed ecological model parameters
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


#---------------------------------------
# Spatial aggregation scales
m_values <- c(1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

#-------------------------------------
# Observation-noise parameters
# Fix observation-noise magnitude at the strong-noise case, where independent observation error showed a clear aggregation benefit.
sigma_obs <- 2

# Spatial correlation of observation error
# rho_obs = 0: completely independent among patches
# rho_obs = 1:  completely shared among patches
rho_obs_values <- c(0, 0.25, 0.5, 0.75, 1)

# -------------------------------------------
# Replicates
n_rep <- 1000
seeds <- seq_len(n_rep)

#--------------------------------
# Run one matched EWS-null replicate
#
# SAME seed is used for:
#
#   - EWS trajectory
#   - stationary Null trajectory
#
# For a given seed:
#
#   - ecological process innovations are matched
#   - underlying observation-error innovations are matched
#   - random spatial sampling configuration is matched
#
# rho_obs only changes how common and local observation
# errors are combined.
#
run_ews_null_correlated_obsnoise_once <- function(seed,  rho_obs,  sigma_obs,  pars, m_values){
    
    
    # Non-stationary EWS trajectory
    out_ews <-  run_one_trajectory_on_ring_network_correlated_obsnoise(seed = seed, pars = pars, 
                                                                       m_values = m_values,
                                                                       sigma_obs = sigma_obs,
                                                                       rho_obs = rho_obs,
                                                                       stationary = FALSE) %>%mutate(type = "EWS")
    
    
    # Matched stationary null trajectory
    out_null <- run_one_trajectory_on_ring_network_correlated_obsnoise(seed = seed, pars = pars, 
                                                                       m_values = m_values,
                                                                       sigma_obs = sigma_obs,
                                                                       rho_obs = rho_obs,
                                                                       stationary = TRUE) %>%mutate(type = "Null")
    
    # --------------------------------------------------------
    # Combine EWS and Null
    # --------------------------------------------------------
    
    out <- bind_rows(out_ews, out_null) %>% mutate(
        # Shared/local ecological process-noise ratio
        sigma_ratio =  pars$sigma_shared / pars$sigma_local)
    
    return(out)
  }


# --------------------------------------------------
# Parameter grid
job_grid <- expand_grid(rho_obs = rho_obs_values, seed = seeds)

# ---------------------------------------------------
# Parallel processing
n_cores <- future::availableCores()
n_workers <- max(1, n_cores - 2)
plan(multisession, workers = n_workers)

#---------------------------------------------
# PARALLEL RUN
#
# Number of jobs:
# 5 rho_obs values x 1000 seeds = 5000 jobs
#
# Each job produces:
# 2 trajectory types
# x 2 strategies
# x 12 m values
#
# = 48 rows
#
# Expected total:
#
# 5000 x 48 = 240,000 rows
#
system.time({results_list <- future_lapply(X = seq_len(nrow(job_grid)), FUN = function(kk){
  run_ews_null_correlated_obsnoise_once(seed =job_grid$seed[kk],
                                        rho_obs = job_grid$rho_obs[kk],
                                        sigma_obs =  sigma_obs,
                                        pars = pars,
                                        m_values =  m_values)}, future.seed = TRUE)})


#----------------------------------------------
# Combine results
results_all <- bind_rows(results_list)

#----------------------------------------------
# Basic checks
nrow(results_all)# Expected:# 240000

# --------------------------------------------
# Save
saveRDS(results_all,here("Results/ring_EWS_tau_with_correlated_obsnoise_1000rep.RDS"))





