# ============================================================
# Simulate ONE N-patch ring trajectory
# 
# PURPOSE: Simulate the stochastic dynamics of N spatial patches arranged
# on a ring, with nearest-neighbour dispersal and a combination of local and shared environmental noise.
#
# The model for patch i is approximately
#
#   dX_i = [lambda(t) * X_i + D * (X_{i-1} + X_{i+1} - 2 X_i)] dt + sigma_local  dW_i + sigma_shared dW_shared
#
# where
#
# X_i           = deviation of patch i from equilibrium
# lambda(t)     = local recovery/stability parameter
# D             = dispersal strength
# sigma_local   = strength of independent local noise
# sigma_shared  = strength of noise shared by all patches
#
# ----------------------------------------------------
# INPUTS:
#
# pars : a list containing
#   pars$N
#       Number of patches in the ring.
#   pars$D
#       Nearest-neighbour dispersal strength.
#   pars$lambda_start
#       Initial value of lambda.
#       Used during burn-in and at the beginning of the observation period.
#   pars$lambda_end
#       Final value of lambda during a non-stationary trajectory.
#   pars$sigma_local
#       Standard deviation of independent local stochastic forcing.
#   pars$sigma_shared
#       Standard deviation of the stochastic forcing shared by all patches.
#   pars$dt
#       Numerical integration time step for Euler-Maruyama.
#   pars$burnin_time
#       Duration of the initial burn-in period.
#   pars$transition_time
#       Duration of the observation / forcing period.
#   pars$sample_every
#       Store the state only every this many integration steps.
#
#
# stationary : logical
#   FALSE: lambda changes linearly from lambda_start to lambda_end.
#   TRUE: lambda remains fixed at lambda_start.
#
#
# seed : optional integer
#   Random-number seed used to make a trajectory reproducible.
#
# ----------------------------------------
# OUTPUT
#
# A list containing:
#
#   X =  Matrix of stored patch states (each patch timeseries along each column; Rows = observation times, Columns = patches)
#   lambda
#       Value of lambda corresponding to each stored observation.
#
#   time
#       Time corresponding to each stored observation.
#
# ============================================================
simulate_ring_trajectory <- function(pars, stationary = FALSE, seed = NULL){
  
  if (!is.null(seed)){set.seed(seed)}
  
  # ----------------------------------------------------------
  # Extract parameters
  # ----------------------------------------------------------
  
  N <- pars$N
  D <- pars$D
  
  lambda_start <- pars$lambda_start
  lambda_end   <- pars$lambda_end
  
  sigma_local  <- pars$sigma_local
  sigma_shared <- pars$sigma_shared
  
  dt <- pars$dt
  
  burnin_time <- pars$burnin_time
  transition_time <- pars$transition_time
  
  sample_every <- pars$sample_every
  
  
  # ----------------------------------------------------------
  # Numerical quantities
  # ----------------------------------------------------------
  
  sqrt_dt <- sqrt(dt) # Brownian increments scale as sqrt(dt)
  n_burn <- round(burnin_time / dt) # Number of Euler-Maruyama steps during burn-in
  n_steps <- round(transition_time / dt) # Number of Euler-Maruyama steps during the observation period
  n_save <- floor(n_steps / sample_every)# Number of observations that will actually be stored
  
  # ----------------------------------------------------------
  # Initial state
  # ----------------------------------------------------------
  
  X <- numeric(N) # Start every patch at the deterministic equilibrium X = 0.
  
  # ==========================================================
  # Burn-in
  #
  # lambda is kept fixed at lambda_start during burn-in.
  # ==========================================================
  
  for (tt in seq_len(n_burn)){
    
    # Periodic neighbours
    X_left <- c(X[N], X[1:(N - 1)])
    X_right <- c(X[2:N], X[1])
    # (left neighbour,focal patch,right neighbour)
    #example: 
    # Patch:        1    2    3    4    5
    # X_left:       50   10   20   30   40
    #               ↓    ↓    ↓    ↓    ↓
    # X:           10   20   30   40   50
    #               ↑    ↑    ↑    ↑    ↑
    # X_right:      20   30   40   50   10
    
    # Deterministic drift
    drift <- lambda_start * X +  D * (X_left +  X_right -  2 * X)
    
    
    # Independent local stochastic increments
    local_noise <- sigma_local * sqrt_dt * rnorm(N, mean=0, sd=1)
    
    
    # One common (or shared) environmental increment
    shared_noise <- sigma_shared * sqrt_dt * rnorm(1, mean=0, sd=1)
    
    
    # Euler-Maruyama update
    X <-  X +  (drift * dt) + local_noise + shared_noise
  }
  
  
  # ==========================================================
  # Observation or transition period
  # ==========================================================
  
  X_store <- matrix(NA_real_, nrow = n_save, ncol = N)# Matrix storing the observed state of all patches. (Rows = sampled times, Columns = patches)
  
  lambda_store <- numeric(n_save) # Store lambda at each sampled observation
  time_store <- numeric(n_save) # Store observation time
  
  save_counter <- 1
  
  for (tt in seq_len(n_steps)){
    
    # --------------------------------------------------------
    # Recovery rate
    # --------------------------------------------------------
    if (stationary){
      lambda_t <- lambda_start# Stationary control trajectory: lambda is held fixed throughout the observation period.
    }else{
      lambda_t <- lambda_start + (lambda_end - lambda_start) * (tt - 1) / (n_steps - 1)# Nonstationary: lambda changes linearly from lambda_start
                                                                                                            # to lambda_end over n_steps.
    }
    # --------------------------------------------------------
    # Periodic neighbours
    # --------------------------------------------------------
    X_left <- c(X[N], X[1:(N - 1)])
    X_right <- c(X[2:N], X[1])
    # --------------------------------------------------------
    # Deterministic drift
    # --------------------------------------------------------
    drift <- lambda_t * X +  D * (X_left +  X_right -  2 * X)
    # --------------------------------------------------------
    # Noise
    # --------------------------------------------------------
    local_noise <-  sigma_local * sqrt_dt * rnorm(N, mean=0, sd=1)
    shared_noise <- sigma_shared * sqrt_dt * rnorm(1, mean=0, sd=1)
    # --------------------------------------------------------
    # Euler-Maruyama step
    # --------------------------------------------------------
    X <- X +  (drift * dt) + local_noise + shared_noise
    # --------------------------------------------------------
    # Store only every sample_every steps
    # --------------------------------------------------------
    if (tt %% sample_every == 0){
      X_store[save_counter, ] <- X
      lambda_store[save_counter] <- lambda_t
      time_store[save_counter] <- tt * dt
      save_counter <- save_counter + 1
    }
  }
  
  out<- list(X = X_store,
             lambda = lambda_store,
             time = time_store)
  
  return(out)
}

# call the function
#sim <- simulate_ring_trajectory(pars = pars, stationary = FALSE, seed = 123)
#dim(sim$X)
#head(sim$lambda)
#tail(sim$lambda)












