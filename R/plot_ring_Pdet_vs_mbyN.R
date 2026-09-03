#rm(list=ls())
library(tidyverse)
library(gridExtra)
library(here)
#========================
results_all<-readRDS(here("Results/ring_EWS_tau_D_sigmaRatio_1000rep.RDS"))
results_all<-results_all%>%filter(sigma_ratio==0.5)

data_ews_D <- results_all %>% filter(type == "EWS")
data_null_D <- results_all %>% filter(type == "Null")

null_thresholds_D <- data_null_D %>% group_by(D, m, m_over_N, strategy) %>% summarise(tau_crit = quantile(tau, probs = 0.95, na.rm = TRUE), .groups = "drop")

# ============================================================
# Detection probability
# ============================================================

data_Pdet_D <- data_ews_D %>% left_join(null_thresholds_D, by = c("D","m","m_over_N","strategy")) %>%
  mutate(detected = tau > tau_crit) %>% group_by(D, m, m_over_N, strategy) %>%
  summarise(Pdet = mean(detected, na.rm = TRUE), n_detected = sum(detected, na.rm = TRUE), n_rep = n(), .groups = "drop")

D_values<-unique(data_Pdet_D$D)
data_Pdet_D$D <- factor(data_Pdet_D$D,levels = D_values)

data_Pdet_difference <- data_Pdet_D %>% select(D,m,m_over_N,strategy,Pdet) %>%
  pivot_wider(names_from = strategy, values_from = Pdet) %>%
  mutate(difference = Consecutive - Random)

head(data_Pdet_difference)

# ============================================================
# PANEL D: Pdet for consecutive sampling
# ============================================================

gD <- ggplot(data_Pdet_difference, aes(x = m_over_N,  y = Consecutive, colour = D,  group = D)) + 
  geom_line(linewidth = 1) +  
  geom_point(size = 2) +
  annotate("text", x = 0.04, y = Inf, label = "D", fontface = "bold", size = 6, vjust = 1.5) +
  labs(title = "Consecutive",
    x = expression("Relative aggregation scale, " * m/N),
    y = expression(P[det]),
    colour = expression(D)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 15, base_family = "sans") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.8))


# ============================================================
# PANEL E: Pdet for random sampling
# ============================================================

gE <- ggplot(data_Pdet_difference, aes(x = m_over_N,  y = Random, colour = D,  group = D)) + 
  geom_line(linewidth = 1) +  
  geom_point(size = 2) +
  annotate("text", x = 0.04, y = Inf, label = "E", fontface = "bold", size = 6, vjust = 1.5) +
  labs(title = "Random",
       x = expression("Relative aggregation scale, " * m/N),
       y = expression(P[det]),
       colour = expression(D)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 15, base_family = "sans") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.8))


# ============================================================
# PANEL F: Median EWS Kendall tau
# ============================================================
# Summary of EWS Kendall tau across replicate trajectories

tau_summary_all_D <- data_ews_D %>% group_by(D, m, m_over_N, strategy) %>% 
  summarise( mean_tau_ews = mean(tau, na.rm = TRUE),
             median_tau_ews = median(tau, na.rm = TRUE),
             sd_tau_ews = sd(tau, na.rm = TRUE),
             n_rep = sum(is.finite(tau)),
             .groups = "drop")


# Keep D ordering consistent with the Pdet plots
tau_summary_all_D$D <- factor(tau_summary_all_D$D, levels = D_values)

gF <- ggplot(tau_summary_all_D, aes(x = m_over_N, y = median_tau_ews, 
                                    colour = D, 
                                    linetype = strategy, 
                                    group = interaction(D, strategy))) +
  
  # Horizontal reference showing no temporal trend
  geom_hline(yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  # Median EWS tau across replicate trajectories
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  # Panel label
  annotate("text", x = 0.04, y = Inf, label = "F", fontface = "bold",  size = 6, vjust = 1.5) +
  labs(title = "EWS trend strength",
    x = expression("Relative aggregation scale, " * m/N),
    y = expression("Median " * tau[EWS]),
    colour = expression(D),
    linetype = "Sampling strategy") +
  scale_x_continuous(limits = c(0, 1),  breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0.5, 1),  breaks = seq(-1, 1, 0.25)) +
  theme_bw(base_size = 15, base_family = "sans") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "inside",
    legend.position.inside = c(0.6, 0.5))


pdf(here("Results/plot_ring_Pdet_vs_mbyN_sigmaRatio0.5.pdf"),  width = 16, height = 5)
grid.arrange(gD, gE, gF, nrow = 1, widths = c(1, 1, 1))
dev.off()


# ============================================================
# Combine EWS and Null Kendall-tau values for supplementary distributions
# ============================================================

tau_distribution_data <- results_all
tau_distribution_data$D <- factor(tau_distribution_data$D,levels = D_values)
tau_distribution_data$m <- factor(tau_distribution_data$m,levels = sort(unique(tau_distribution_data$m)))
tau_distribution_data$type <- factor(tau_distribution_data$type,levels = c("Null", "EWS"))


# ============================================================
# Prepare threshold + Pdet information for each histogram facet
# ============================================================
#
# Each facet corresponds to:
#
#   D x m x strategy
#
# For each facet we need:
#
#   tau_crit = 95th percentile of Null tau
#   Pdet     = fraction of EWS tau > tau_crit
#
null_thresholds_D$D <- factor(null_thresholds_D$D,levels = D_values)
hist_annotations <- null_thresholds_D %>% left_join(data_Pdet_D %>% 
                                                      select(D,m,m_over_N,strategy,Pdet),
                                                    by = c("D","m","m_over_N","strategy"))


# Keep factor ordering identical to histogram data

hist_annotations$m <- factor(hist_annotations$m,levels = sort(unique(results_all$m)))
# ============================================================
# Supplement:
# Histogram of Kendall's tau — Consecutive sampling
# ============================================================

hist_annotations_con <- hist_annotations %>%filter(strategy == "Consecutive")

g_tau_hist_consecutive <-
  tau_distribution_data %>%
  filter(strategy == "Consecutive") %>%
  ggplot(aes(x = tau, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.45,bins = 30) +
  geom_vline(data = hist_annotations_con, aes(xintercept = tau_crit),
  inherit.aes = FALSE,
  linetype = "dashed",
  linewidth = 0.7) +
  geom_text(data = hist_annotations_con,
  aes(x = -Inf,y = Inf,
      label = paste0("P[det] == ",sprintf("%.2f", Pdet))),
  inherit.aes = FALSE,
  parse = TRUE,
  hjust = -0.08,
  vjust = 1.25,
  size = 3.2) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  labs(title = "Consecutive sampling",
  x = expression("Kendall's " * tau), y = "Frequency",fill = NULL) +
  facet_grid(D ~ m,labeller = label_both) +
  theme_bw(base_size = 15,base_family = "sans") +
  theme(legend.position = "top")

pdf(here("Results/plot_ring_tauHist_consecutivesampling_sigmaRatio0.5.pdf"),  width = 12,  height = 6)
print(g_tau_hist_consecutive)
dev.off()

# ============================================================
# Supplement:
# Histogram of Kendall's tau — Random sampling
# ============================================================
hist_annotations_rand <- hist_annotations %>%filter(strategy == "Random")

g_tau_hist_random <-
  tau_distribution_data %>%
  filter(strategy == "Random") %>%
  ggplot(aes(x = tau, fill = type)) +
  geom_histogram(position = "identity", alpha = 0.45,bins = 30) +
  geom_vline(data = hist_annotations_rand, aes(xintercept = tau_crit),
             inherit.aes = FALSE,
             linetype = "dashed",
             linewidth = 0.7) +
  geom_text(data = hist_annotations_rand,
            aes(x = -Inf,y = Inf,
                label = paste0("P[det] == ",sprintf("%.2f", Pdet))),
            inherit.aes = FALSE,
            parse = TRUE,
            hjust = -0.08,
            vjust = 1.25,
            size = 3.2) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1)) +
  labs(title = "Random sampling",
       x = expression("Kendall's " * tau), y = "Frequency",fill = NULL) +
  facet_grid(D ~ m,labeller = label_both) +
  theme_bw(base_size = 15,base_family = "sans") +
  theme(legend.position = "top")

pdf(here("Results/plot_ring_tauHist_randomsampling_sigmaRatio0.5.pdf"),  width = 12,  height = 6)
print(g_tau_hist_random)
dev.off()
