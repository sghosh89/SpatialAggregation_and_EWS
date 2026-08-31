# ============================================================
# Calculate the EWS trend statistic for ONE aggregated time series
# ============================================================
#
# Calculate a variance-based early-warning signal (EWS) for one aggregated time series Y.
#
# The function:
#
#   1. calculates rolling variance through time,
#   2. removes positions where rolling variance is not available,
#   3. calculates Kendall's tau between time and rolling variance.
#
# A positive tau means that rolling variance tends to increase through time.
#
# -------------------------------------------------------------------------------
# INPUTS
#
# Y :Numeric vector containing ONE aggregated time series obtained from either consecutive or random spatial sampling.
# window_fraction :  Fraction of the complete time series used as the rolling window length.
#--------------
# OUTPUT

# tau : A single numeric value containing Kendall's rank correlation between time and rolling variance.
#
# ============================================================
source(here("R/rolling_variance.R"))

get_ews_tau <- function(Y,  window_fraction = 0.5){
  
  V_roll <- rolling_variance(Y, window_fraction)
  valid <- is.finite(V_roll)
  
  # ----------------------------------------------------------
  # STEP 3: Calculate Kendall's tau
  # ----------------------------------------------------------
  
  tau <- suppressWarnings(
    
    cor(which(valid),# which(valid) returns the position indices for which valid == TRUE : acts as the time ordering of the rolling-variance estimates.
            V_roll[valid], method = "kendall")
  )
  
  
  # suppressWarnings() prevents warning messages from being
  # printed in unusual cases, for example when all rolling
  # variance values are identical and the correlation cannot
  # be estimated normally.
  

  
  return(tau)
}