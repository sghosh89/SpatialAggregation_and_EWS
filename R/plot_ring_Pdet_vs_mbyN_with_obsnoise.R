#rm(list=ls())

library(tidyverse)
library(gridExtra)
library(here)

# ============================================================
# Read observation-noise simulation results
results_all <- readRDS(here("Results/ring_EWS_tau_with_obsnoise_1000rep.RDS"))
data_ews_obs <- results_all %>% filter(type == "EWS")
data_null_obs <- results_all %>% filter(type == "Null")

# ============================================================
# Null detection threshold
#
# IMPORTANT:
# Calculate a separate 95% null threshold for every:
#
#   sigma_obs x m x strategy
#
# because observation noise may alter the null distribution
# of Kendall's tau.
# ============================================================

null_thresholds_obs <- data_null_obs %>% group_by(sigma_obs, m, m_over_N, strategy) %>%
  summarise(tau_crit = quantile(tau, probs = 0.95, na.rm = TRUE), .groups = "drop")

# ============================================================
# Detection probability
# EWS is detected when: tau_EWS > tau_crit

data_Pdet_obs <- data_ews_obs %>%  left_join(null_thresholds_obs,  by = c("sigma_obs","m","m_over_N","strategy")) %>%
  mutate(detected = tau > tau_crit) %>%
  group_by(sigma_obs, m, m_over_N, strategy) %>%
  summarise(Pdet = mean(detected, na.rm = TRUE),
    n_detected = sum(detected, na.rm = TRUE),
    n_rep = n(),
    .groups = "drop"
  )

# ============================================================
# Make sigma_obs a factor for plotting

sigma_obs_values <- sort(unique(data_Pdet_obs$sigma_obs))
data_Pdet_obs$sigma_obs <- factor(data_Pdet_obs$sigma_obs,levels = sigma_obs_values)

# ============================================================
# Put Consecutive and Random in separate columns

data_Pdet_difference <- data_Pdet_obs %>%
  select(sigma_obs, m, m_over_N, strategy, Pdet) %>%
  pivot_wider(names_from = strategy, values_from = Pdet) %>%
  mutate(difference = Consecutive - Random)

# ============================================================
# PANEL A: Detection probability — Consecutive sampling

gA <- ggplot(data_Pdet_difference,aes(x = m_over_N, y = Consecutive,colour = sigma_obs, group = sigma_obs)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  annotate("text",x = 0.04, y = Inf,label = "A",fontface = "bold",size = 6, vjust = 1.5) +
  labs(title = "Consecutive",
    x = expression(
      "Relative aggregation scale, " * m/N
    ),
    y = expression(P[det]),
    colour = expression(sigma[obs])
  ) +
  scale_x_continuous(limits = c(0, 1),breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1),breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 15,base_family = "sans") +
  theme(plot.title = element_text(face = "bold",hjust = 0.5),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.75))


# ============================================================
# PANEL B
# Detection probability — Random sampling
# ============================================================

gB <- ggplot(data_Pdet_difference,aes(x = m_over_N, y = Random, colour = sigma_obs, group = sigma_obs)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  annotate("text",x = 0.04, y = Inf,label = "B",fontface = "bold",size = 6,vjust = 1.5) +
  labs(title = "Random",
    x = expression("Relative aggregation scale, " * m/N),
    y = expression(P[det]),
    colour = expression(sigma[obs])) +
  scale_x_continuous(limits = c(0, 1),breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1),breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 15,base_family = "sans") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.75))


# ============================================================
# Combine plots

pdf(here("Results/indep_obs_noise/plot_ring_Pdet_vs_mbyN_obsnoise.pdf"),width = 10,height = 5)
grid.arrange(gA,gB, nrow = 1, widths = c(1, 1))
dev.off()
#===========================================================================================================================================


# ============================================================
# Combine EWS and Null Kendall-tau values for supplementary distributions
# ============================================================

# ============================================================
# Null thresholds for D = 0.1
# ============================================================

null_thresholds <- results_all %>%
  filter(
    type == "Null",
    D == 0.1
  ) %>%
  group_by(m, strategy, sigma_obs) %>%
  summarise(
    tau_crit = quantile(tau, probs = 0.95, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# Prepare data
# ============================================================

tau_distribution_data <- results_all %>%
  filter(D == 0.1) %>%
  mutate(
    sigma_obs = factor(
      sigma_obs,
      levels = sigma_obs_values
    ),
    m = factor(
      m,
      levels = sort(unique(m))
    ),
    type = factor(
      type,
      levels = c("Null", "EWS")
    )
  )

# ============================================================
# Detection probability
# D = 0.1, Consecutive & random sampling
# ============================================================

Pdet_labels_con <- results_all %>%filter(D == 0.1,strategy == "Consecutive",type == "EWS") %>%
  # Attach the corresponding 95% null threshold
  left_join(null_thresholds %>% filter(strategy == "Consecutive"),
    by = c("m", "strategy", "sigma_obs")) %>%
  # Detection probability:
  # fraction of EWS replicates exceeding the 95% null threshold
  group_by(m, sigma_obs) %>%  summarise(Pdet = mean(tau > tau_crit, na.rm = TRUE),
    .groups = "drop") %>%
  # Match facet-variable formatting used in the histogram
  mutate(m = factor(m,levels = sort(unique(results_all$m))),
    sigma_obs = factor(sigma_obs,levels = sigma_obs_values),
    # Plotmath label: P_det = value
    label = paste0("P[det] == ",sprintf("%.2f", Pdet)))

Pdet_labels_rand <- results_all %>%filter(D == 0.1,strategy == "Random",type == "EWS") %>%
  # Attach the corresponding 95% null threshold
  left_join(null_thresholds %>% filter(strategy == "Random"),
            by = c("m", "strategy", "sigma_obs")) %>%
  # Detection probability:
  # fraction of EWS replicates exceeding the 95% null threshold
  group_by(m, sigma_obs) %>%  summarise(Pdet = mean(tau > tau_crit, na.rm = TRUE),
                                        .groups = "drop") %>%
  # Match facet-variable formatting used in the histogram
  mutate(m = factor(m,levels = sort(unique(results_all$m))),
         sigma_obs = factor(sigma_obs,levels = sigma_obs_values),
         # Plotmath label: P_det = value
         label = paste0("P[det] == ",sprintf("%.2f", Pdet)))

# ============================================================
# Histogram of Kendall's tau
# Consecutive sampling, D = 0.1
# ============================================================

g_tau_hist_consecutive <-  tau_distribution_data %>% filter(D == 0.1, strategy == "Consecutive") %>%
  ggplot(aes(x = tau, fill = type)) +
  # Null and EWS distributions
  geom_histogram(position = "identity", alpha = 0.45,  bins = 30) +
  # 95% null threshold for each m × sigma_obs combination
  geom_vline(data = null_thresholds %>%
      filter(strategy == "Consecutive"),
    aes(xintercept = tau_crit),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.6) +
  # Detection probability in the top-left corner
  geom_text(data = Pdet_labels_con,  aes(x = -Inf,y = Inf,label = label),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = -0.08,
    vjust = 1.25,
    size = 3.2
  ) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  labs(title = "Consecutive sampling",
    subtitle = expression(D == 0.1),
    x = expression("Kendall's " * tau),
    y = "Frequency",
    fill = NULL) +
  facet_grid(sigma_obs ~ m,labeller = label_both) +
  theme_bw(base_size = 15,base_family = "sans") +  
  theme(legend.position = "top")

pdf(here("Results/indep_obs_noise/plot_ring_tauHist_consecutivesampling_D0.1_obsnoise.pdf"),width = 15,height = 9)
print(g_tau_hist_consecutive)
dev.off()

#g_tau_hist_consecutive
#Figure Sx. Distribution of variance-based early-warning trends under observational noise and spatial aggregation. 
#Distributions of Kendall’s \(\tau\) for the rolling variance under approaching-transition (EWS; blue) and stationary 
#null (red) simulations for consecutive spatial sampling at \(D=0.1\). Columns show increasing aggregation scale \(m\), 
#and rows show increasing observational-noise intensity \(\sigma_{\mathrm{obs}}\). Dashed vertical lines indicate 
#the 95th percentile of the corresponding null distribution, used as the detection threshold. The detection 
#probability \(P_{\mathrm{det}}\), shown in each panel, is the proportion of EWS realizations with \(\tau\) 
#exceeding this threshold. Although EWS distributions are shifted toward positive \(\tau\), substantial overlap 
#with the upper tail of the null distribution results in moderate detection probabilities (\(P_{\mathrm{det}}<0.5\)). 
#Detection probability nevertheless remains broadly stable across aggregation scales and observational-noise levels, 
#except for reduced detectability at small \(m\) under strong observational noise.


g_tau_hist_random <-  tau_distribution_data %>% filter(D == 0.1, strategy == "Random") %>%
  ggplot(aes(x = tau, fill = type)) +
  # Null and EWS distributions
  geom_histogram(position = "identity", alpha = 0.45,  bins = 30) +
  # 95% null threshold for each m × sigma_obs combination
  geom_vline(data = null_thresholds %>%
               filter(strategy == "Random"),
             aes(xintercept = tau_crit),
             inherit.aes = FALSE,
             linetype = "dashed",
             linewidth = 0.6) +
  # Detection probability in the top-left corner
  geom_text(data = Pdet_labels_rand,  aes(x = -Inf,y = Inf,label = label),
            inherit.aes = FALSE,
            parse = TRUE,
            hjust = -0.08,
            vjust = 1.25,
            size = 3.2
  ) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  labs(title = "Random sampling",
       subtitle = expression(D == 0.1),
       x = expression("Kendall's " * tau),
       y = "Frequency",
       fill = NULL) +
  facet_grid(sigma_obs ~ m,labeller = label_both) +
  theme_bw(base_size = 15,base_family = "sans") +  
  theme(legend.position = "top")

pdf(here("Results/indep_obs_noise/plot_ring_tauHist_randomsampling_D0.1_obsnoise.pdf"),width = 15,height = 9)
print(g_tau_hist_random)
dev.off()

