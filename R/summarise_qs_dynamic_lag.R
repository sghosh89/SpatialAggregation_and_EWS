# ============================================================
# Quasi-stationary versus dynamic retention across parameter space
# ============================================================
#
# PURPOSE:
# Quantify the error introduced by the quasi-stationary approximation when lambda changes through time.
#
# The previous script: compare_retention_qs_vs_dynamic.R
# compares quasi-stationary and dynamic R_m(lambda) for one representative parameter combination: D = 0.1, sigma_shared / sigma_local = 0.5, m = 20
#
# Here, that comparison is extended across the same parameter space used in the quasi-stationary R_m(lambda) analysis.
#
# For each combination of D, gamma and m, calculate: lag = R_m^QS - R_m^dynamic, relative lag = (R_m^QS - R_m^dynamic) / R_m^QS
# Positive values indicate that the quasi-stationary approximation overestimates the dynamically realized variance retention.
#
#--------------------------------------------------------------------------------------
# INPUT:
# Functions already loaded through these files:
source(here("R/ring_retention.R"))
source(here("R/dynamic_modal_variances.R"))

# Parameter values are chosen to match get_Rm_vs_lambda.R:
#
#   N = 100
#   D = 0, 0.02, 0.1, 0.5
#   gamma_val = sigma_shared / sigma_local = 0, 0.5, 2
#   m = 5, 20, 50, 100
#
# The lambda trajectory matches the EWS simulations:
#
#   lambda_start = -0.5
#   lambda_end   = -0.05
#   transition_time = 200
#   dt = 0.05
#
#---------------------------------------
# OUTPUT:
# data_qs_dynamic_lag:
#   Full lambda-resolved dataset containing dynamic and quasi-stationary R_m and their discrepancy.

# lag_summary:
#   Summary for each D x gamma_val x m x sampling strategy, including final and maximum absolute/relative lag.

# ============================================================
# Parameter values
# ============================================================
#
# These match the parameter values used in get_Rm_vs_lambda.R
# and the dynamic EWS simulations.
# ============================================================

N <- 100
sigma_local <- 1
lambda_start <- -0.5
lambda_end   <- -0.05
transition_time <- 200
dt <- 0.05
D_values_lag <- c(0, 0.02, 0.1, 0.5)
gamma_values_lag <- c(0, 0.5, 2)
m_values_lag <- c(5, 20, 50, 100)

# ============================================================
# Parameter grid for dynamic covariance calculation
# ============================================================
#
# Modal dynamics depend on D and gamma, but NOT on m.
# Therefore dynamic_modal_variances() only needs to be solved once for each D x gamma combination.
# ============================================================

parameter_grid_lag <- expand_grid(D = D_values_lag,
                                  gamma_val = gamma_values_lag)


# ============================================================
# Calculate dynamic and quasi-stationary R_m
# ============================================================

lag_results <- vector("list", length = nrow(parameter_grid_lag))


for (ii in seq_len(nrow(parameter_grid_lag))) {
  
  # Current dispersal strength
  D_now <- parameter_grid_lag$D[ii]
  
  # Current shared/local forcing ratio
  gamma_now <- parameter_grid_lag$gamma_val[ii]
  
  # Convert gamma into shared-noise SD
  sigma_shared_now <- gamma_now * sigma_local
  
  
  # ----------------------------------------------------------
  # Solve dynamic modal variances
  #
  # This is done once for the current D x gamma combination.
  # ----------------------------------------------------------
  
  modal_now <- dynamic_modal_variances(
    N = N,
    D = D_now,
    lambda_start = lambda_start,
    lambda_end = lambda_end,
    transition_time = transition_time,
    dt = dt,
    sigma_local = sigma_local,
    sigma_shared = sigma_shared_now
  )
  
  
  # ----------------------------------------------------------
  # Calculate retention for each aggregation scale m
  # ----------------------------------------------------------
  
  results_m <- vector("list", length = length(m_values_lag))
  
  
  for (jj in seq_along(m_values_lag)){
    
    m_now <- m_values_lag[jj]
    
    # Dynamic retention
    Rm_dyn_con <- dynamic_retention_consecutive(modal_now, m = m_now)
    Rm_dyn_rand <- dynamic_retention_random(modal_now, m = m_now)
    
    # Quasi-stationary retention
    # Evaluate the stationary analytical expression at each instantaneous lambda along the dynamic trajectory.
     Rm_QS_con <- sapply(modal_now$lambda, function(lambda_now){
        ring_retention(N = N,
                       m = m_now,
                       lambda = lambda_now,
                       D = D_now,
                       sigma_local = sigma_local,
                       sigma_shared = sigma_shared_now,
                       strategy = "consecutive")})
    
    
    Rm_QS_rand <- sapply(modal_now$lambda, function(lambda_now){
      ring_retention(N = N,
                     m = m_now,
                     lambda = lambda_now,
                     D = D_now,
                     sigma_local = sigma_local,
                     sigma_shared = sigma_shared_now,
                     strategy = "random")})
    
    
    # ========================================================
    # Store both sampling strategies in long format
    # ========================================================
    
    results_m[[jj]] <- bind_rows(
      
      tibble(
        time = modal_now$time,
        lambda = modal_now$lambda,
        D = D_now,
        gamma_val = gamma_now,
        m = m_now,
        m_over_N = m_now / N,
        strategy = "Consecutive",
        Rm_dynamic = Rm_dyn_con,
        Rm_QS = Rm_QS_con
      ),
      
      tibble(
        time = modal_now$time,
        lambda = modal_now$lambda,
        D = D_now,
        gamma_val = gamma_now,
        m = m_now,
        m_over_N = m_now / N,
        strategy = "Random",
        Rm_dynamic = Rm_dyn_rand,
        Rm_QS = Rm_QS_rand
      )
    )
  }
  
  
  # Combine all m for the current D x gamma combination
  lag_results[[ii]] <- bind_rows(results_m)
}


# Combine all parameter combinations
data_qs_dynamic_lag <- bind_rows(lag_results)


# ============================================================
# Calculate QS-dynamic discrepancy
# ============================================================
#
# Absolute lag:
#
#   lag = R_m^QS - R_m^dynamic
#
# Relative lag:
#
#   relative_lag =
#       (R_m^QS - R_m^dynamic) / R_m^QS
#
# ============================================================

data_qs_dynamic_lag <- data_qs_dynamic_lag %>%
  mutate(lag = Rm_QS - Rm_dynamic,
    relative_lag =lag / Rm_QS)


# ============================================================
# Consistency check at the beginning of the transition
# The dynamic modal variances are initialized at their stationary values for lambda_start.
# Therefore: R_m^dynamic(t = 0) = R_m^QS(lambda_start) for every parameter combination.

initial_lag_check <- data_qs_dynamic_lag %>%
  group_by(
    D,
    gamma_val,
    m,
    strategy
  ) %>%
  slice_head(n = 1) %>%
  ungroup()


max_initial_error <- max(
  abs(initial_lag_check$lag),
  na.rm = TRUE
)

max_initial_error


# ============================================================
# Summarize QS-dynamic discrepancy
# ============================================================
#
# final_lag:
#   discrepancy at lambda_end
#
# max_lag:
#   largest QS - dynamic difference during the transition
#
# final_relative_lag:
#   proportional discrepancy at lambda_end
#
# max_relative_lag:
#   largest proportional discrepancy during the transition
# ============================================================

lag_summary <- data_qs_dynamic_lag %>% group_by(D,gamma_val,m,m_over_N,strategy) %>%
  summarise(final_lag = last(lag),
            max_lag =  max(lag, na.rm = TRUE),
            final_relative_lag = last(relative_lag),
            max_relative_lag = max(relative_lag, na.rm = TRUE), .groups = "drop")

lag_summary %>%
  arrange(desc(max_relative_lag))


# Save results
saveRDS(data_qs_dynamic_lag, here("Results/data_qs_dynamic_lag.RDS"))
saveRDS(lag_summary, here("Results/summary_qs_dynamic_lag.RDS"))


lag_summary %>%
  mutate(
    difference =
      max_relative_lag - final_relative_lag
  ) %>%
  summarise(
    max_difference = max(abs(difference)),
    mean_difference = mean(abs(difference))
  )

#lag_summary # therefore in this table you can keep only the final_relative_lag
lag_summary<-readRDS(here("Results/summary_qs_dynamic_lag.RDS"))
g1 <- ggplot(lag_summary, aes(x = m_over_N, y = final_relative_lag, colour = factor(D), linetype = strategy, group = interaction(D, strategy))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ gamma_val, nrow = 1,labeller = label_both)+ #, labeller = labeller(gamma_val = function(x) paste0("\u03B3 = ", x))) +
  labs(x = expression("Relative aggregation scale, " * m/N),
       y = "Final relative discrepancy",
       colour = expression(D),
       linetype = "Sampling strategy") +
  theme_bw(base_size = 14, base_family = "sans") +
  theme(legend.position = "top")

pdf(here("Results/Final_relative_lag_vs_mbyN_qs_and_dynamic_approach.pdf"),  width = 10, height = 5)
print(g1)
dev.off()

#Figure Sx. Finite-rate effects on spatial variance retention. Relative discrepancy between quasi-stationary and dynamically calculated 
#variance retention, \(R_m\), at the end of the non-stationary transition (\(\lambda=-0.05\)), shown as a function of relative aggregation scale
#\(m/N\). Relative lag was calculated as \((R_m^{QS}-R_m^{dyn})/R_m^{QS}\), where \(R_m^{QS}\) is obtained from the stationary covariance 
#evaluated at the instantaneous recovery rate and \(R_m^{dyn}\) is obtained from the time-dependent modal variance equations. 
#Panels correspond to different ratios of shared to local stochastic forcing, \(\gamma=\sigma_{\mathrm{shared}}/\sigma_{\mathrm{local}}\); 
#colours indicate dispersal strength \(D\), and line types distinguish consecutive and random sampling. The quasi-stationary approximation 
#increasingly overestimates dynamically realized retention under weak shared forcing, whereas the discrepancy is strongly reduced as 
#shared forcing increases and is negligible under strong shared forcing. For \(D=0\), the relative lag is zero.











