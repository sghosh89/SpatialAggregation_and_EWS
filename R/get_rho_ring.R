# ============================================================
# Exact stationary covariance for the N-patch ring
# N= patch numbers, must be >2
# r = shortest network distance between two patches
# lambda < 0 = local recovery eigenvalue
# D = nearest-neighbour dispersal strength
# sigma_local = sigma_l
# sigma_shared = sigma_s
# ============================================================

get_cov_ring <- function(N, r, lambda, D, sigma_local = 1, sigma_shared = 0){
  
  if (N < 2)
    stop("N must be at least 2.")
  
  if (lambda >= 0)
    stop("lambda must be negative for a stable equilibrium.")
  
  if (D < 0)
    stop("D must be non-negative.")
  
  # epsilon = |lambda| = -lambda
  epsilon <- abs(lambda)
  
  # Fourier mode indices
  k <- 0:(N - 1)
  
  # Recovery strength of each spatial mode
  alpha_k <- epsilon + (4 * D * sin(pi * k / N)^2)
  
  # Shared environmental forcing contributes equally to covariance at every spatial separation.
  shared_component <- sigma_shared^2 / (2 * epsilon)
  
  # Patch-specific stochastic forcing produces distance-dependent covariance through dispersal.
  local_component <- sigma_local^2 / (2 * N) *  sum(cos(2 * pi * k * r / N) / alpha_k)
  
  cov_ring<-shared_component + local_component
  
  return(cov_ring)
}


# ============================================================
# Exact spatial correlation rho_r = C_r / C_0
# ============================================================
# N= patch numbers, must be >2
# r = shortest network distance between two patches
# lambda < 0 = local recovery eigenvalue
# D = nearest-neighbour dispersal strength
# sigma_local = 
# sigma_shared = sigma_s
get_rho_ring <- function(N, r, lambda, D, sigma_local = 1, sigma_shared = 0){
  
  C0 <- get_cov_ring(N = N, r = 0, lambda = lambda, D = D, sigma_local = sigma_local, sigma_shared = sigma_shared)
  
  Cr <- get_cov_ring(N = N, r = r, lambda = lambda, D = D, sigma_local = sigma_local, sigma_shared = sigma_shared)
  
  rho_ring<- Cr/C0
  
  return(rho_ring)
}
