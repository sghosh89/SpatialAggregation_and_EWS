# ============================================================
# # Calculate rolling variance for ONE aggregated time series
# ============================================================
# INPUTS

# Y :   Numeric vector containing ONE time series.
#       In our analysis, Y will usually be an aggregated time series obtained from either: aggregate_consecutive_ring(), or aggregate_random_ring()
#
# window_fraction :  Fraction of the complete time series used as the  rolling-window length.
#-------------------------------------
# OUTPUT
#
# V_roll : Numeric vector containing the rolling variance.  It has the SAME length as Y. The initial entries are NA because there are not yet enough
#   observations to construct a complete rolling window.
#
library(RcppRoll)
rolling_variance <- function(Y, window_fraction = 0.5){
  
  n <- length(Y)
  
  w <- floor(window_fraction * n)
  
  if (w < 3){
    stop("Rolling window is too short.")
  }
  
  V_roll <- RcppRoll::roll_var(Y,  n = w, align = "right", fill = NA)
  
  return(V_roll)
}