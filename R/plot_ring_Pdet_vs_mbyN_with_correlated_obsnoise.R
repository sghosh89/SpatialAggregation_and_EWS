library(tidyverse)
library(gridExtra)
library(here)

# ============================================================
# Read observation-noise simulation results
results_all <- readRDS(here("Results/ring_EWS_tau_with_correlated_obsnoise_1000rep.RDS"))
data_ews <- results_all %>% filter(type == "EWS")
data_null <- results_all %>% filter(type == "Null")

# ============================================================
# get detection threshold
#
# IMPORTANT:
# threshold is recalculated separately for every
# rho_obs x m x strategy combination
null_thresholds <- data_null %>%  group_by(rho_obs, m, m_over_N, strategy) %>%
  summarise(tau_crit = quantile(tau, probs = 0.95, na.rm = TRUE),
    .groups = "drop")

# ============================================================
# Detection probability: Pdet = P(tau_EWS > tau_crit)
data_Pdet <- data_ews %>%
  left_join(null_thresholds,
    by = c("rho_obs", "m", "m_over_N", "strategy")) %>%
  mutate(detected = tau > tau_crit) %>%
  group_by(rho_obs, m, m_over_N, strategy) %>%
  summarise(Pdet = mean(detected, na.rm = TRUE),
            n_detected = sum(detected, na.rm = TRUE),
            n_rep = sum(is.finite(tau)),
            .groups = "drop")

# ============================================================
# Probability of a positive temporal variance trend
# Ppositive = P(tau_EWS > 0)
data_Ppositive <- data_ews %>% group_by(rho_obs, m, m_over_N, strategy) %>%
      summarise(Ppositive = mean(tau > 0, na.rm = TRUE),
                median_tau_ews = median(tau, na.rm = TRUE),
                mean_tau_ews = mean(tau, na.rm = TRUE), .groups = "drop")

# ============================================================
# Combine
data_detection_summary <- data_Pdet %>% left_join(data_Ppositive,  by = c("rho_obs", "m", "m_over_N", "strategy"))
# Inspect results
print(data_detection_summary, n = Inf)

# Plot Pdet and Ppositive together
data_plot <- data_detection_summary %>% select(rho_obs,  m_over_N, strategy, Pdet, Ppositive) %>%
  pivot_longer(cols = c(Pdet, Ppositive),
    names_to = "metric", values_to = "probability") %>%
  mutate(metric = recode(metric, Pdet = "Detection probability",
      Ppositive = "Positive variance trend"))


g_corr_obs <- ggplot(data_plot, aes(x = m_over_N,  y = probability, 
                                    colour = factor(rho_obs), 
                                    linetype = metric, 
                                    group = interaction(rho_obs,  metric))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ strategy) +
  labs(x = expression("Relative aggregation scale, " * m/N),
    y = "Probability",
    colour = expression(rho[obs]),
    linetype = NULL) +
  scale_x_continuous(limits = c(0, 1),breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1),breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 15,base_family = "sans") +
  theme(legend.position = "top")

pdf(here("Results/plot_ring_Pdet_Ppositive_vs_mbyN.pdf"),  width = 12, height = 6)
print(g_corr_obs)
dev.off()

#Figure X. Effects of spatially correlated observation error on the detectability of variance-based early-warning signals 
#under spatial aggregation. Detection probability (solid lines), \(P_{\mathrm{det}}=P(\tau_{\mathrm{EWS}}>\tau_{\mathrm{crit}})\), 
#and the probability of a positive temporal trend in variance (dashed lines), \(P(\tau_{\mathrm{EWS}}>0)\), are shown as 
#functions of the relative aggregation scale \(m/N\) for consecutive (A) and random (B) spatial sampling. Colours indicate
#the spatial correlation of observation error, \(\rho_{\mathrm{obs}}\), ranging from independent (\(\rho_{\mathrm{obs}}=0\)) 
#to fully correlated (\(\rho_{\mathrm{obs}}=1\)) errors. Observation-error magnitude was fixed at \(\sigma_{\mathrm{obs}}=2\), 
#with \(D=0.1\) and \(\gamma=\sigma_{\mathrm{shared}}/\sigma_{\mathrm{local}}=0.5\). Detection thresholds \(\tau_{\mathrm{crit}}\) 
#were defined separately for each parameter combination as the 95th percentile of Kendall's \(\tau\) obtained from the 
#corresponding stationary null simulations. Each estimate was based on 1000 stochastic replicates. Increasing spatial 
#correlation in observation error progressively reduced both the prevalence of positive variance trends and their formal 
#detection, while the improvement in \(P_{\mathrm{det}}\) with spatial aggregation observed for independent errors 
#disappeared as observation errors became increasingly correlated.



