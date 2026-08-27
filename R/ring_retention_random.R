# ============================================================
# Mean pairwise correlation over ALL distinct pairs
# in the N-patch ring
#
# This gives rho_bar_N needed for expected random sampling.
# N= patch numbers, must be >2
# m = CONSECUTIVE observed patches
# lambda < 0 = local recovery eigenvalue
# D = nearest-neighbour dispersal strength
# sigma_local = sigma_l
# sigma_shared = sigma_s
# ============================================================
#source(here("R/get_rho_ring.R"))
ring_mean_pairwise_rho <- function(N, m, lambda, D, sigma_local = 1, sigma_shared = 0) {
  
  # Generate every distinct pair of patch indices.
  # R uses 1,...,N here; this does not affect distances.
  pairs <- combn(1:N, 2)
  
  # Ordinary index separation for every pair
  s <- abs(pairs[1, ] - pairs[2, ])
  
  # Shortest distance around the ring
  ring_distance <- pmin(s, N - s) #this is a vectors -> minimum at EACH combo of s, N-s 
  
  # Compute rho_r only once for each distinct distance
  unique_distances <- sort(unique(ring_distance))
  
  rho_by_distance <- sapply(unique_distances, function(r){get_rho_ring(N = N, r = r, lambda = lambda, D = D, sigma_local = sigma_local, sigma_shared = sigma_shared)})
  
  
  # match(x, table) asks:
  #
  # "For each value in x, at which POSITION does that value
  #  occur in table?"
  #
  # In our N = 4 example:
  #
  # ring_distance    = 1, 2, 1, 1, 2, 1
  # unique_distances = 1, 2
  #
  # Therefore:
  #
  # distance 1 occurs at position 1 in unique_distances
  # distance 2 occurs at position 2 in unique_distances
  #
  # So match() returns: distance_index = 1, 2, 1, 1, 2, 1
  distance_index<- match(ring_distance, unique_distances)
  rho_all_pairs <- rho_by_distance[distance_index]
  
  # Mean correlation across all N choose 2 patch pairs
  mean(rho_all_pairs)
}

# ============================================================
# Expected retention for m patches selected uniformly
# at random without replacement
#
# This is E[R_m^random], not retention of one particular
# randomly drawn set.
# ============================================================

ring_retention_random <- function(N, m, lambda, D, sigma_local = 1, sigma_shared = 0){
  
  if (m < 1 || m > N)
    stop("m must satisfy 1 <= m <= N.")
  
  if (m == 1)
    return(1)
  
  rho_bar_N <- ring_mean_pairwise_rho(N = N, lambda = lambda, D = D, sigma_local = sigma_local, sigma_shared = sigma_shared)
  
  Rm_expected <- (1 / m) + (((m - 1) / m) *  rho_bar_N)
  
  return(Rm_expected)
}