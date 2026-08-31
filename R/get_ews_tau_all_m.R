# ============================================================
# Calculate Kendall's tau for ALL aggregation scales m
# ============================================================
# Take a list of aggregated time series, Y_list, and calculate the variance-based EWS statistic (Kendall's tau) separately for every m.
#
# Each element of Y_list corresponds to one value of m.
#
# For example: Y_list[["2"]]   -> aggregated time series for m = 2 patches
#--------------------------------------------------------------------
# INPUTS
#
# Y_list : A named list of aggregated time series. The names of the list must correspond to the aggregation scales m.
# window_fraction : Fraction of the time series used for the rolling-variance window.  Default = 0.5.
#
# ----------------------------------------------------------------------------
# OUTPUT
#
# A data.frame with two columns:
#
#   m : aggregation scale, tau : Kendall's tau measuring the temporal trend in rolling variance for that aggregation scale
#
# ============================================================
source("R/get_ews_tau.R")
get_ews_tau_all_m <- function(Y_list,  window_fraction = 0.5){
  
  m_values <- as.integer(names(Y_list))

  tau_values <- vapply(Y_list, 
                       FUN = get_ews_tau, 
                       FUN.VALUE = numeric(1),# Tell R that calculate_ews_tau() must return exactly ONE numeric value for each time series. This makes vapply() stricter and safer than sapply().
                       window_fraction = window_fraction)
 
  
  df<-data.frame(m = m_values,  tau=unname(tau_values))
  rownames(df) <- NULL
  return(df)
}