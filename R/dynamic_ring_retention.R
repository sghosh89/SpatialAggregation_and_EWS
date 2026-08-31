# ============================================================
# Reconstruct C_r(t) from modal variances
# Purpose:
#   Reconstruct the time-dependent spatial covariance C_r(t)
#   between two patches separated by ring distance r, using the
#   dynamic Fourier-mode variances V_k(t).
#
# Inputs:
#   modal_result : Output from dynamic_modal_variances();
#                  contains V(t), Fourier modes k, time, lambda(t)
#   r            : Ring distance between two patches
#
# Output:
#   Numeric vector giving C_r(t) at every time point
# ============================================================

dynamic_ring_covariance <- function(modal_result,  r){
  
  V <- modal_result$V
  k <- modal_result$k
  N <- ncol(V)
  
  weights <- cos(2 * pi * k * r / N)
  
  # C_r(t) for every time point
  Cr <- as.vector(V %*% weights) / N
  
  return(Cr)
}

# ============================================================
# Dynamic R_m(t): consecutive sampling
# Purpose:
#   Calculate the time-dependent variance-retention factor R_m(t)   when m consecutive patches on the ring are spatially averaged.
#
#   R_m(t) measures how much of the single-patch variance C_0(t) remains after averaging over m neighbouring patches:
#
#       Var(mean of m patches) = C_0(t) * R_m(t)
#
# Inputs:
#   modal_result : Output from dynamic_modal_variances()
#   m            : Number of consecutive patches being averaged
#
# Output:
#   Numeric vector giving R_m(t) at every time point
# ============================================================

dynamic_retention_consecutive <- function(modal_result, m){
  
  N <- ncol(modal_result$V)
  
  # Local variance C_0(t)
  C0 <- dynamic_ring_covariance(modal_result, r = 0)
  
  if (m == 1){
    return(rep(1,length(C0)))
  }
  
  s_values <- 1:(m - 1)
  r_values <- pmin(s_values, N - s_values)
  pair_counts <-  m - s_values
  
  
  # Start with diagonal contribution
  Rm <-  rep(1 / m, length(C0))
  
  
  # Add covariance contributions
  for (ii in seq_along(r_values)){
    Cr <- dynamic_ring_covariance(modal_result, r = r_values[ii])
    rho_r <- Cr / C0
    Rm <- Rm + (2 / m^2) * pair_counts[ii] * rho_r
  }
  
  return(Rm)
}

# ============================================================
# Dynamic expected R_m(t): random sampling
# Purpose:
#   Calculate the expected time-dependent variance-retention
#   factor R_m(t) when m patches are selected randomly from the
#   N-patch ring.
#
#   Unlike consecutive sampling, random sampling does not have a
#   fixed spatial configuration. Therefore the expected R_m(t)
#   is calculated using the mean pairwise correlation across all
#   possible distinct patch pairs on the ring.
#
# Inputs:
#   modal_result : Output from dynamic_modal_variances()
#   m            : Number of randomly sampled patches
#
# Output:
#   Numeric vector giving the expected R_m(t) at every time point
#   under random sampling
# ============================================================

dynamic_retention_random <- function(modal_result,  m){
  
  N <- ncol(modal_result$V)
  
  C0 <- dynamic_ring_covariance(modal_result, r = 0)
  
  if (m == 1){
    return(rep(1, length(C0)))
  }
  
  
  # All distinct patch pairs on the full ring
  pairs <- combn(1:N, 2)
  separation <- abs(pairs[1, ] - pairs[2, ])
  ring_distance <- pmin(separation, N - separation)
  
  
  # Compute rho_r(t) only once for each distinct distance
  unique_r <- sort(unique(ring_distance))
  
  
  # Matrix: rows = time, cols = ring distance
  rho_mat <- sapply(unique_r, function(r){dynamic_ring_covariance(modal_result, r = r) / C0})
  
  
  # Number of pairs at each ring distance
  pair_counts <- table(factor(ring_distance, levels = unique_r))
  
  # Mean pairwise correlation at each time point
  rho_bar <- as.vector(rho_mat %*% as.numeric(pair_counts)) / choose(N, 2)
  
  # Expected random retention
  Rm_random <-  1 / m + ((m - 1) / m) * rho_bar
  
  return(Rm_random)
}









