#rm(list=ls())
library(tidyverse)
library(gridExtra)
library(here)
#========================
results_all<-readRDS(here("Results/ring_EWS_tau_D_sigmaRatio_1000rep.RDS"))

data_null_all <- results_all %>% filter(type == "Null")
data_ews_all <- results_all %>% filter(type == "EWS")


# Null threshold specific to every parameter combination
null_thresholds_all <- data_null_all %>%
  group_by(D, sigma_ratio, m, m_over_N, strategy) %>%
  summarise(tau_crit = quantile(tau, probs = 0.95, na.rm = TRUE), .groups = "drop")


# Detection probability
data_Pdet_all <- data_ews_all %>% left_join(null_thresholds_all,  by = c("D", "sigma_ratio", "m", "m_over_N", "strategy")) %>%
  mutate(detected = tau > tau_crit) %>%
  group_by(D, sigma_ratio, m, m_over_N, strategy) %>%
  summarise(Pdet = mean(detected, na.rm = TRUE),
    n_detected = sum(detected, na.rm = TRUE),
    n_rep = n(),
    .groups = "drop")

ggplot(data_Pdet_all, aes(x = m_over_N, y = Pdet, colour = factor(D), linetype = strategy, group = interaction(D, strategy))) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  facet_wrap(~ sigma_ratio, nrow = 1,
    labeller = labeller(sigma_ratio = function(x) paste0("\u03B3 = ", x))) +
  labs(x = expression("Relative aggregation scale, " * m/N),
    y = expression(P[det]),
    colour = expression(D),
    linetype = "Sampling strategy") +
  scale_y_continuous(limits = c(0, 1),breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 14, base_family = "sans") +
  theme(legend.position = "top")


data_Pdet_all <- data_Pdet_all %>%
  mutate(
    gamma_plot = factor(
      sigma_ratio,
      levels = c(0, 0.25, 0.5, 1, 2)
    )
  )

library(RColorBrewer)

g_Pdet_heat <- ggplot(data_Pdet_all, aes(x = factor(m_over_N), y = gamma_plot, fill = Pdet)) +
  geom_tile() +
  facet_grid(strategy ~ D, labeller = labeller(D = function(x) paste0("D = ", x))) +
  scale_fill_gradientn(
    colours = brewer.pal(9, "YlGnBu"),
    limits = c(0.35, 0.55),
    breaks = c(0.35, 0.40, 0.45, 0.50, 0.55),
    name = expression(P[det])
  ) +
  scale_x_discrete(breaks = c("0.01", "0.2", "0.4", "0.6", "0.8", "1")) +
  labs(x = expression("Relative aggregation scale, " * m/N),
    y = expression(gamma)) +
  theme_bw(base_size = 14, base_family = "sans") +
  theme(legend.position = "top", panel.grid = element_blank())
g_Pdet_heat

pdf(here("Results/plot_Pdet_heatmap.pdf"),  width = 12, height = 6)
print(g_Pdet_heat)
dev.off()
#Figure Sx. Early-warning signal detection probability across the spatial parameter space. 
#Detection probability (\(P_{\rm det}\)) as a function of relative aggregation scale \(m/N\), the ratio of shared to local 
#stochastic forcing \(\gamma=\sigma_{\mathrm{shared}}/\sigma_{\mathrm{local}}\), dispersal strength \(D\), and spatial 
#sampling strategy. Detection was defined by comparing Kendall's \(\tau\) of the rolling variance from non-stationary 
#trajectories with the 95th percentile of the corresponding stationary-null distribution for each parameter combination. 
#each estimate was based on 1000 replicate trajectories. Across the parameter space examined, \(P_{\rm det}\) remained 
#within a relatively narrow range (0.371–0.532), indicating that variance-based EWS detection was comparatively 
#insensitive to spatial aggregation and sampling geometry under the ideal-observation model.