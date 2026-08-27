# ============================================================
# Retention for m CONSECUTIVE observed patches in a N patch ring
#
# The observed patches are a contiguous block embedded in the full N-patch ring.
# N= patch numbers, must be >2
# m = CONSECUTIVE observed patches
# lambda < 0 = local recovery eigenvalue
# D = nearest-neighbour dispersal strength
# sigma_local = sigma_l
# sigma_shared = sigma_s
# ============================================================
#source(here("R/get_rho_ring.R"))
ring_retention_consecutive <- function(N, m, lambda, D, sigma_local = 1, sigma_shared = 0){
  
  if (m < 1 || m > N)
    stop("m must satisfy 1 <= m <= N.")
  
  # No aggregation when only one patch is observed
  if (m == 1)
    return(1)
  
  # Separations within the observed block
  s_values <- 1:(m - 1)
  
  # Convert index separation s to shortest ring distance r
  r_values <- pmin(s_values, N - s_values) # pmin used for each pair index
  
  # Exact correlation corresponding to each ring distance
  rho_values <- sapply(r_values, function(r){get_rho_ring(N=N, r=r, lambda=lambda, D=D, sigma_local = sigma_local, sigma_shared = sigma_shared)})
  
  # Number of observed pairs having index separation s
  pair_counts <- m - s_values
  
  # Exact retention factor
  Rm <- (1 / m)  +   ((2 / m^2) * sum(pair_counts * rho_values))
  
  return(Rm)
}