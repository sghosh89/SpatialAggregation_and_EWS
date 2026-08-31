# ===========================================================================================================
# Aggregate ONE N-patch trajectory over consecutive patches
# ===========================================================================================================
# Starting from one simulated N-patch trajectory, construct an aggregated time series for each requested aggregation scale m.
#
# For a given m, the function averages the first m consecutive patches:
#
#              1
#   Y_m(t) =  ---  sum X_i(t)
#              m   i=1,...,m
#
# Thus:
#
#   m = 1  -> use patch 1 only
#   m = 2  -> average patches 1 and 2
#   m = 3  -> average patches 1, 2 and 3
#   ...
#   m = N  -> average all N patches
#
# Because the patches lie on a homogeneous ring, choosing patches 1,...,m is equivalent to choosing any other block of m consecutive patches.
#-----------------------------
# INPUTS
# X: Matrix containing ONE simulated N-patch trajectory.  (Rows = observation times, Columns = patches)
#
# m_values: Vector containing the aggregation scales to calculate.
#---------------------------
# OUTPUT

# A list called Y_list: Each element contains one aggregated time series.
#
# For example:
#
#   Y_list[["1"]] = time series for patch 1
#
#   Y_list[["2"]] = mean time series of patches 1:2

# Each time series has the same number of observations as the number of rows in X.
#
# ============================================================
aggregate_consecutive_ring <- function(X, m_values){
  
  N <- ncol(X) #total number of patches
  
  if (any(m_values < 1 | m_values > N)){    # Check if every m must satisfy 1 <= m <= N
    stop("All m must satisfy 1 <= m <= N.")
  }
  
  # ----------------------------------------------------------
  # Construct aggregated time series for every m
  # ----------------------------------------------------------
  
  # lapply() goes through each value contained in m_values.
  #
  # For each m:
  #
  #   1. select the first m consecutive columns (patches) of X
  #
  #   2. calculate the mean across those m patches separately
  #      for every time point
  #
  #   3. return the resulting aggregated time series
  #
  # The resulting time series are collected into a list.
  
  Y_list <- lapply(m_values, function(m){
      
      
      # ------------------------------------------------------
      # Select the first m consecutive patches
      # ------------------------------------------------------
      
      # seq_len(m) generates:
      #
      # m = 1  -> 1
      # m = 2  -> 1, 2
      # m = 3  -> 1, 2, 3
      # ...
      #
      # Therefore:
      #
      # X[, seq_len(m)]
      #
      # selects all observation times (all rows) but only patches 1,...,m (the first m columns).
      
      
      # ------------------------------------------------------
      # Average across the selected patches
      # ------------------------------------------------------
      
      # rowMeans() calculates one mean for every row of X.
      #
      # Since each row represents one observation time, this gives the spatial mean across the selected m consecutive patches at each time point.
      X_sub<-X[ ,seq_len(m), drop = FALSE]# Keep as a matrix even when m = 1.Without drop = FALSE, selecting only one column could cause R to simplify the matrix into a vector.
      rowMeans(X_sub)
    }
  )
  
  names(Y_list) <- m_values# Give each list: timeseries, the corresponding m value as its name
  
  return(Y_list)
}

#X <- matrix(
#  c(
#    10, 20, 30, 40, 50,
#    12, 22, 32, 42, 52,
#    14, 24, 34, 44, 54
#  ),
#  nrow = 3,
#  byrow = TRUE
#)

#X

#aggregate_consecutive(X, m_values = c(1,2,3))









