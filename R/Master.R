library(here)
# =======================================
source(here("R/check_general_aggregation_framework.R"))# General aggregation framework

source(here("R/ring_retention.R"))# Stationary analytical theory for the ring
source(here("R/plot_ring_Rm_vs_mbyN.R"))# R_m versus spatial sampling fraction m/N for consecutive and random sampling
source(here("R/get_Rm_vs_lambda.R"))# State dependence of retention: stationary R_m increases as lambda approaches zero

# Non-equilibrium validation
source(here("R/dynamic_modal_variances.R"))# Solve time-dependent modal variance equations
source(here("R/compare_retention_qs_vs_dynamic.R"))# Compare actual dynamic R_m(lambda) against quasi-stationary R_m(lambda)
                                                      # QS versus dynamic R_m for D = 0.1, gamma = 0.5, m = 20
source(here("R/summarise_qs_dynamic_lag.R"))# Generalize QS-dynamic discrepancy across parameter space (put in suppmat, maybe)


# Baseline stochastic EWS detection: Ecological process noise only
# Does aggregation-induced variance dilution translate into reduced EWS detectability?
source(here("R/sim_detection_threshold_in_ring.R"))# Simulate non-stationary EWS trajectories and matched stationary null trajectories across parameter space
source(here("R/plot_ring_Pdet_vs_mbyN.R"))# Plot P_det versus m/N this is plotting results only for sigmaRatio=0.5, # maybe delete this
source(here("R/plot_Pdet_heatmap.R"))# this is heatmap across parameter space

#==============================================
# start with observation noise model
#source(here("R/sim_detection_threshold_obsnoise_in_ring.R"))
#source(here("R/plot_ring_Pdet_vs_mbyN_with_obsnoise.R"))

#=================================================
# Observation-error extension
# Spatially correlated observation error: When does observation structure alter the baseline detectability result?
source(here("R/sim_detection_threshold_correlated_obsnoise_in_ring.R"))
source(here("R/plot_ring_Pdet_vs_mbyN_with_correlated_obsnoise.R"))# main text fig
source(here("R/plot_tau_histograms_correlated_obs_noise.R"))# suppmat figs.

#================================================






