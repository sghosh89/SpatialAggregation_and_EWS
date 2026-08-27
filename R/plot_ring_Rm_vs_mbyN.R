# ============================================================
# FIGURE 3A-B
# Analytical variance retention for the ring vs. m/N
#
# A: m consecutive observed patches
# B: m randomly chosen patches (expected retention)
# ============================================================
library(here)
library(tidyverse)
library(gridExtra)
# -----------------------------
# Baseline parameters
# -----------------------------

N <- 100
lambda <- -0.2
sigma_local <- 1
# Start with moderate shared forcing
sigma_shared <- 0.5
#sigmaRatio<- sigma_shared/sigma_local

# Different dispersal strengths
D_values <- c(0, 0.02, 0.1, 0.5)

# Aggregation scales
m_values <- 1:N

# -----------------------------
# Store results
# -----------------------------

data_retention <- data.frame()

for (D in D_values){
  
  # Consecutive sampling
  R_consecutive <- sapply(
    m_values,
    function(m){
      ring_retention(N = N, m = m, lambda = lambda, D = D, sigma_local = sigma_local, sigma_shared = sigma_shared,strategy = "consecutive")
    }
  )
  
  
  # Random sampling
  R_random <- sapply(
    m_values,
    function(m){
      ring_retention(N = N, m = m, lambda = lambda, D = D, sigma_local = sigma_local, sigma_shared = sigma_shared,strategy = "random")
    }
  )
  
  
  # Save in long format
  temp_consecutive <- data.frame(m = m_values,
                                 m_over_N = m_values / N,
                                 D = D,
                                 strategy = "Consecutive",
                                 retention = R_consecutive)
  
  temp_random <- data.frame(m = m_values,
                            m_over_N = m_values / N,
                            D = D,
                            strategy = "Random",
                            retention = R_random)
  
  data_retention <- rbind(data_retention, temp_consecutive,  temp_random)
  
}



# Treat D as a factor so ggplot gives one distinct curve
# for each dispersal strength
data_retention$D <- factor(data_retention$D)
data_retention$strategy <- factor(data_retention$strategy, levels = c("Consecutive","Random"))

# ---------------------------------------------------------
# Panel A: Consecutive sampling
# ---------------------------------------------------------

gA <- data_retention %>%
  filter(strategy == "Consecutive") %>%
  ggplot(aes(x = m_over_N, y = retention, colour = D)) +
  geom_line(linewidth = 1) +
  labs(title = "Consecutive",
       x = expression("Relative aggregation scale, " * m/N),
    y = expression("Variance retention, " * R[m]),
    colour = expression(D)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")


# ---------------------------------------------------------
# Panel B: Random sampling
# ---------------------------------------------------------

gB <- data_retention %>%
  filter(strategy == "Random") %>%
  ggplot(aes(x = m_over_N, y = retention, colour = D)) +
  geom_line(linewidth = 1) +
  labs(title = "Random",
       x = expression("Relative aggregation scale, " * m/N),
    y = expression("Expected variance retention, " * E(R[m])),
    colour = expression(D)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")


# ---------------------------------------------------------
# Panel C: Difference
# ---------------------------------------------------------

data_difference <- data_retention %>%
  select(m, m_over_N, D, strategy, retention) %>%
  pivot_wider(names_from = strategy,
    values_from = retention) %>%
  mutate(difference = Consecutive - Random)


gC <- ggplot(data_difference,  aes(x = m_over_N, y = difference, colour = D)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_line(linewidth = 1) +
  labs(title = "Sampling-strategy difference",
       x = expression("Relative aggregation scale, " * m/N),
    y = expression(R[m]^Consecutive - E(R[m]^Random)),
    colour = expression(D)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw(base_size = 15)+
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.98),
    legend.justification = c(1, 1)
  )


# ---------------------------------------------------------
# Combine A, B, C in one row
# ---------------------------------------------------------

pdf(here("Results/plot_ring_Rm_vs_mbyN_sigmaRatio0.5.pdf"),  width = 15, height = 5)
gridExtra::grid.arrange(gA, gB, gC, nrow = 1, widths = c(1, 1, 1))
dev.off()
