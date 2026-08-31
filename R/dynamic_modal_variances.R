source(here("R/dynamic_ring_retention.R"))
# ===================================================================================================================================================
# Dynamic modal variances for a nonstationary N-patch ring
#
# Purpose: 
#   Compute the time-dependent variance V_k(t) of each spatial Fourier mode k when the local recovery rate lambda(t) changes through time.
#
# For each Fourier mode, the variance obeys dV_k/dt = -2 * alpha_k(t) * V_k + q_k^2
#                                                where, alpha_k(t) = |lambda(t)| + 4D sin^2(pi k / N)
# is the effective recovery rate of mode k.
#
# The function uses exactly the same linear lambda(t) trajectory as the nonstationary EWS simulations.
# ====================================================================================================================================================
# Input:
#N = Number of patches in the ring
#D =  Dispersal/diffusion strength
#lambda_start =  Initial local recovery rate
#lambda_end =  Final recovery rate, closer to zero
#transition_time =  Duration over which lambda changes
#dt =  Numerical integration time step
#sigma_local =  SD of independent/local noise
#sigma_shared =  SD of spatially shared noise


dynamic_modal_variances <- function(N = 100, D = 0.1, lambda_start = -0.5, lambda_end = -0.05, 
                                    transition_time = 200, dt = 0.05, sigma_local = 1, sigma_shared = 0.5){
  
  # ----------------------------------------------------------
  # Time grid
  # ----------------------------------------------------------
  n_steps <- round(transition_time / dt)# Number of integration steps during the transition.
  time <- seq(0, transition_time, length.out = n_steps + 1)
  
  # Linear nonstationary trajectory for lambda: lambda(t) = lambda_start + (lambda_end - lambda_start) * t / transition_time
  # Thus lambda gradually approaches zero during the simulation. This is the same lambda(t) schedule used in the EWS simulations.
  lambda_t <- lambda_start + ((lambda_end - lambda_start) * time / transition_time)
  
  
  # ----------------------------------------------------------
  # Fourier modes
  # ----------------------------------------------------------
  
  k <- 0:(N - 1)
  alpha_start <-  abs(lambda_start) +  (4 * D * sin(pi * k / N)^2) # Initial modal recovery rates
  
  # ----------------------------------------------------------
  # Noise intensity q_k^2
  #
  # Local noise contributes to all modes.
  # Shared forcing contributes only to k = 0 (uniform mode).
  # ----------------------------------------------------------
  q2 <- rep(sigma_local^2, N)
  q2[1] <- sigma_local^2 + (N * sigma_shared^2)
  
  # ----------------------------------------------------------
  # Initial modal variances:
  # At t = 0, assume the system has already equilibrated under lambda = lambda_start.
  # For a stationary OU mode, dY_k = -alpha_k Y_k dt + q_k dW_k, the stationary variance is V_k* = q_k^2 / (2 alpha_k).
  # Therefore the dynamic calculation starts from the stationary variance associated with lambda_start.
  # ----------------------------------------------------------
  V <- matrix(NA_real_, nrow = n_steps + 1, ncol = N)
  V[1, ] <- q2 /(2 * alpha_start)
  
  # ----------------------------------------------------------
  # Integrate modal variance ODEs
  # Euler method, using same dt as ecological simulation
  # ----------------------------------------------------------
  
  for (tt in seq_len(n_steps)){
    alpha_t <- abs(lambda_t[tt]) +  (4 * D * sin(pi * k / N)^2)
    V[tt + 1, ] <-  V[tt, ] + (-2 * alpha_t * V[tt, ] + q2) * dt
  }
  
  out<- list(time = time,
            lambda = lambda_t,
            V = V, # V is a matrix with:   rows    = time points,   columns = Fourier modes k = 0,...,N-1
                        # Therefore V[t, k+1] gives the variance of Fourier mode k at the corresponding time point.
            k = k)
    
  return(out)
}
#=============================================





