library(here)
library(tidyverse)
library(future)
library(future.apply)
source(here("R/run_one_trajectory_on_ring_network.R"))

# parameter list
pars <- list(
  
  # Spatial system
  N = 100,
  D = 0.1,
  
  # Stability
  lambda_start = -0.5,
  lambda_end   = -0.05,
  
  # Stochastic forcing
  sigma_local  = 1,
  sigma_shared = 0.5,
  
  # Numerical integration
  dt = 0.05,
  burnin_time = 50,
  transition_time = 200,
  
  # Save only every kth numerical step
  # This considerably reduces memory and EWS computation
  sample_every = 5,
  
  # Rolling EWS
  window_fraction = 0.5
)


m_values <- c(1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

D_values <- c(0, 0.02, 0.1, 0.5)

n_rep <- 1000
seeds <- seq_len(n_rep)


sigma_ratio_values <- c(0, 0.25, 0.5, 1, 2)# sigma_s/sigma_l

run_ews_null_once <- function(seed,D,sigma_ratio,pars,m_values){
  
  # Local copy
  pars_run <- pars
  
  # Set dispersal
  pars_run$D <- D
  
  # Keep sigma_local fixed
  pars_run$sigma_local <- 1
  
  # sigma_ratio = sigma_shared / sigma_local
  pars_run$sigma_shared <-  sigma_ratio * pars_run$sigma_local
  
  
  # ----------------------------------------------------------
  # EWS trajectory
  # ----------------------------------------------------------
  
  ews <- run_one_trajectory_on_ring_network(seed = seed, pars = pars_run, m_values = m_values, stationary = FALSE)
  ews$type <- "EWS"
  
  # ----------------------------------------------------------
  # Matched stationary null
  # ----------------------------------------------------------
  
  null <- run_one_trajectory_on_ring_network(seed = seed, pars = pars_run, m_values = m_values, stationary = TRUE)
  null$type <- "Null"
  
  # ----------------------------------------------------------
  # Combine
  # ----------------------------------------------------------
  
  out <- bind_rows(ews, null)
  
  out$D <- D
  out$sigma_ratio <- sigma_ratio
  
  return(out)
}

job_grid <- expand.grid(D = D_values,  sigma_ratio = sigma_ratio_values,  seed = seeds)
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
                                  run_ews_null_once(seed = job_grid$seed[kk],
                                                    D = job_grid$D[kk],
                                                    sigma_ratio = job_grid$sigma_ratio[kk],
                                                    pars = pars,
                                                    m_values = m_values)},
                                future.seed = TRUE)
})

results_all <- bind_rows(results_list)

saveRDS(results_all, here("Results/ring_EWS_tau_D_sigmaRatio_1000rep.RDS"))
