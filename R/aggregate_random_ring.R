# ============================================================
# Aggregate ONE trajectory for all m: RANDOM spatial sampling
# ============================================================
# Starting from one simulated N-patch trajectory X, construct an aggregated time series for each requested number of randomly sampled patches, m.

# IMPORTANT TRICK:
# The function first generates ONE random ordering of all N patches. For each m, it then uses the first m patches from this same random ordering.
# Therefore the random samples are NESTED.
#
# Example:
#
# random_order = 4, 1, 7, 2, 5, 3, 6
#
# m = 1  -> patch 4
# m = 2  -> patches 4, 1
# m = 3  -> patches 4, 1, 7
# m = 4  -> patches 4, 1, 7, 2
#
# Thus increasing m adds patches to the existing sample rather than drawing a completely new random sample.
#
#-------------------------------------------------------
# INPUT

# X: Matrix containing ONE simulated N-patch trajectory.  (Rows = observation times, Columns = patches)
# m_values: Vector containing the aggregation scales to calculate.
#--------------------------------------------------------------
#
# OUTPUT
#
# A list containing:
#   Y : A list of aggregated time series, one for each m.
#   order : The random ordering of patches used to construct the nested random samples.
#
# ============================================================
aggregate_random_ring <- function(X, m_values){
  
  
  N <- ncol(X) #total number of patches
  
  if (any(m_values < 1 | m_values > N)){    # Check if every m must satisfy 1 <= m <= N
    stop("All m must satisfy 1 <= m <= N.")
  }
  
  # ----------------------------------------------------------
  # Generate ONE random ordering of all N patches
  # ----------------------------------------------------------
  
  random_order <- sample(seq_len(N), size = N, replace = FALSE) # draw N patches without replacement from 1:N 
  
  # Example: for N = 5, random_order might be: 4  2  5  1  3
  # This random ordering will be used for ALL values of m within this trajectory.
  
  # ----------------------------------------------------------
  # Construct the aggregated time series for every m
  # ----------------------------------------------------------
 
  Y_list <- lapply(m_values, function(m){
      
      
      # ------------------------------------------------------
      # Select the m randomly ordered patches to be observed
      # ------------------------------------------------------
      
      observed <- random_order[seq_len(m)]# Select positions 1 through m from random_order.
        
      # Example: Suppose      
      # random_order = 4, 2, 5, 1, 3      
      # If m = 1: observed = 4      
      # If m = 2: observed = 4, 2     
      # If m = 3: observed = 4, 2, 5
      #
      # Thus the samples are nested as m increases.
      
      
      # ------------------------------------------------------
      # Calculate the spatially aggregated time series
      # ------------------------------------------------------
      # Select only the columns corresponding to the randomly selected patches.
      X_sub<-X[ ,observed, drop = FALSE]# Keep as a matrix even when m = 1.Without drop = FALSE, selecting only one column could cause R to simplify the matrix into a vector.
      rowMeans(X_sub)
    }
  )
  
  
  names(Y_list) <- m_values
  
  out<-list(Y = Y_list,# List containing one aggregated time series for each m
       order = random_order# Random ordering of patches used to construct all of the nested random samples
  )
  
  return(out)
}

