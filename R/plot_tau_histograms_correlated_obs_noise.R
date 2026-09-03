# ============================================================
# Histograms of Kendall's tau under spatially correlated
# observation error
#
# Shows:
#   1. Null tau distribution
#   2. EWS tau distribution
#   3. 95% Null threshold
#   4. Pdet = P(tau_EWS > tau_crit)
#   5. P+   = P(tau_EWS > 0)
# ============================================================
library(tidyverse)
library(here)
# ============================================================
# Read results
results_all <- readRDS(here("Results/ring_EWS_tau_with_correlated_obsnoise_1000rep.RDS"))
data_ews <- results_all %>% filter(type == "EWS")
data_null <- results_all %>% filter(type == "Null")
# ============================================================
# Calculate 95% Null threshold
#
# IMPORTANT:
# Separate threshold for each:
#
# rho_obs x m x strategy
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
# Calculate probability of positive EWS trend
# Ppositive = P(tau_EWS > 0)
# ============================================================
data_Ppositive <- data_ews %>% group_by(rho_obs, m, m_over_N, strategy) %>%
  summarise(Ppositive = mean(tau > 0, na.rm = TRUE),
            median_tau_ews = median(tau, na.rm = TRUE),
            mean_tau_ews = mean(tau, na.rm = TRUE), .groups = "drop")

# ============================================================
# Combine threshold, Pdet and Ppositive
# This data frame will be used for vertical lines and text annotations.
# ============================================================

hist_annotations <- null_thresholds %>%  left_join(data_Pdet,  by = c("rho_obs", "m", "m_over_N", "strategy")) %>%
  left_join(data_Ppositive,   by = c("rho_obs", "m", "m_over_N", "strategy"))


# ============================================================
# Factor ordering
rho_obs_values <- sort(unique(results_all$rho_obs))
m_values <- sort(unique(results_all$m))


results_all <- results_all %>%  mutate(rho_obs = factor(rho_obs, levels = rho_obs_values),
                                       m = factor(m, levels = m_values),
                                       type = factor(type,  levels = c("Null",  "EWS")))


hist_annotations <- hist_annotations %>%  mutate(rho_obs = factor(rho_obs,  levels = rho_obs_values),
                                                 m = factor(m, levels = m_values))


# ============================================================
# CONSECUTIVE SAMPLING
hist_annotations_con <- hist_annotations %>%filter(strategy == "Consecutive")


g_tau_hist_consecutive <-  results_all %>% filter(strategy == "Consecutive") %>%
  ggplot(aes(x = tau, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.45, bins = 30) +
  geom_vline(xintercept = 0,  linetype = "dotted",  linewidth = 0.4) +
  geom_vline(data = hist_annotations_con,
  aes(xintercept = tau_crit), inherit.aes = FALSE, linetype = "dashed", linewidth = 0.7) +
  geom_text(data = hist_annotations_con,  aes(x = -1, y = Inf,
                                              label = paste0("Pdet = ", sprintf("%.2f", Pdet),
                                                             "\n",
                                                             "P+ = ", sprintf("%.2f", Ppositive))),
            inherit.aes = FALSE, hjust = 0, vjust = 1.15, size = 3) +
  facet_grid(rho_obs ~ m,
  labeller = labeller(rho_obs = function(x){paste0("rho_obs = ", x)},
    m = function(x){paste0("m = ", x)})) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  labs(title ="Consecutive sampling", x = expression("Kendall's " * tau), y = "Frequency", fill = NULL) +
  theme_bw(base_size = 14,  base_family = "sans") +
  theme(legend.position = "top", strip.text =  element_text(face = "bold"))


pdf(here("Results/plot_tauHist_correlated_obsnoise_consecutive.pdf"), width = 18, height = 10)
print(g_tau_hist_consecutive)
dev.off()

# ============================================================
# RANDOM SAMPLING

hist_annotations_rand <- hist_annotations %>% filter(strategy == "Random")
g_tau_hist_random <-  results_all %>% filter(strategy == "Random") %>%
  ggplot(aes(x = tau, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.45, bins = 30) +
  geom_vline(xintercept = 0,  linetype = "dotted",  linewidth = 0.4) +
  geom_vline(data = hist_annotations_rand,
             aes(xintercept = tau_crit), inherit.aes = FALSE, linetype = "dashed", linewidth = 0.7) +
  geom_text(data = hist_annotations_rand,  aes(x = -1, y = Inf,
                                              label = paste0("Pdet = ", sprintf("%.2f", Pdet),
                                                             "\n",
                                                             "P+ = ", sprintf("%.2f", Ppositive))),
            inherit.aes = FALSE, hjust = 0, vjust = 1.15, size = 3) +
  facet_grid(rho_obs ~ m,
             labeller = labeller(rho_obs = function(x){paste0("rho_obs = ", x)},
                                 m = function(x){paste0("m = ", x)})) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  labs(title ="Random sampling", x = expression("Kendall's " * tau), y = "Frequency", fill = NULL) +
  theme_bw(base_size = 14,  base_family = "sans") +
  theme(legend.position = "top", strip.text =  element_text(face = "bold"))
  

pdf(here("Results/plot_tauHist_correlated_obsnoise_random.pdf"), width = 18, height = 10)
print(g_tau_hist_random)
dev.off()

