# This script will check: Rm(t) itself change as λ(t)→0 ?
#==================================================================
library(here)
library(tidyverse)
library(gridExtra)
source(here("R/ring_retention.R"))
#==================================================================
N <- 100
sigma_local <- 1
# Same lambda range used in the EWS simulations
lambda_values <- seq(-0.5,-0.05,length.out = 100)

# Start with a few representative aggregation scales
m_values_check <- c(5,  20,  50,  100)

# Representative dispersal strengths
D_values_check <- c(0, 0.02, 0.1, 0.5)

# Representative shared/local forcing ratios
gamma_values_check <- c(0, 0.5, 2)# sigma_shared = sigma_local*gamma_value

#parameter grid
data_Rm_lambda <- expand_grid(lambda = lambda_values,
                              D = D_values_check,
                              gamma_val = gamma_values_check,
                              m = m_values_check,
                              strategy = c("Consecutive", "Random"))


data_Rm_lambda <- data_Rm_lambda %>% rowwise() %>%  mutate(sigma_shared =  gamma_val * sigma_local,  
                                                           Rm = ring_retention(N = N, m = m, lambda = lambda, D = D, 
                                                                               sigma_local = sigma_local, 
                                                                               sigma_shared = sigma_shared,
                                                                               strategy = ifelse(strategy == "Consecutive", "consecutive", "random"))) %>% 
  ungroup() %>%
  mutate(m = factor(m, levels = m_values_check),
    gamma_val = factor(gamma_val, levels = gamma_values_check))

saveRDS(data_Rm_lambda,here("Results/data_Rm_lambda.RDS"))

data_Rm_lambda_con<-data_Rm_lambda%>%filter(strategy=="Consecutive")
data_Rm_lambda_ran<-data_Rm_lambda%>%filter(strategy=="Random")

ggplot(data_Rm_lambda_con, aes(x = lambda, y = Rm, colour = factor(D), group = D)) +
  geom_line(linewidth = 0.9) +
  facet_grid(gamma_val ~ m, labeller = label_both) +
  labs(x = expression(lambda),
    y = expression(R[m]),
    colour = expression(D)) +
  theme_bw(base_size = 15,base_family = "sans") +
  theme(legend.position = "top")


ggplot(data_Rm_lambda_ran, aes(x = lambda, y = Rm, colour = factor(D), group = D)) +
  geom_line(linewidth = 0.9) +
  facet_grid(gamma_val ~ m, labeller = label_both) +
  labs(x = expression(lambda),
       y = expression(R[m]),
       colour = expression(D)) +
  theme_bw(base_size = 15,base_family = "sans") +
  theme(legend.position = "top")


# ============================================================
# Difference in retention between observation strategies
#
# Delta R_m(lambda)
#   = R_m^Consecutive(lambda) - R_m^Random(lambda)
# ============================================================

data_delta_Rm_lambda <- data_Rm_lambda %>%select(lambda, D, gamma_val, m, strategy, Rm)%>%
  pivot_wider(names_from = strategy, values_from = Rm) %>%
  mutate(delta_Rm = Consecutive - Random)


# ============================================================
# Plot Delta R_m against lambda
# ============================================================

g_delta_Rm_lambda <- ggplot(data_delta_Rm_lambda, aes(x = lambda, y = delta_Rm, colour = factor(D), group = D)) +
  # Reference line: no difference between strategies
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  # Same layout as the R_m vs lambda figures
  facet_grid(gamma_val ~ m,
    labeller = label_both) +
  labs(x = expression(lambda),
    y = expression(Delta * R[m] ==  R[m]^Consecutive -  E(R[m]^Random)), colour = expression(D)) +
  theme_bw(base_size = 15,base_family = "sans") +
  theme(legend.position = "top")

g_delta_Rm_lambda

pdf(here("Results/plot_delta_Rm_vs_lambda.pdf"),  width = 12, height = 6)
print(g_delta_Rm_lambda)
dev.off()

# interpretation of the plot: 
#The key result is that ΔRm(λ)=Rm_Consecutive(λ)−E[Rm_Random(λ)] is positive whenever dispersal creates appreciable short-range covariance, 
#but its magnitude depends strongly on m, γ, and proximity to criticality.
#At γ=0, the sampling-strategy effect is strongest. For m=5, strong dispersal (D=0.5) gives a substantial advantage to consecutive sampling, 
#and that advantage rises sharply as λ→0−. This is exactly what we would expect: dispersal creates distance-dependent covariance, 
#and consecutive sampling preferentially retains neighboring, highly correlated patches.

#Notice the especially interesting result: dΔRm/dλ>0as λ→0− for D>0, particularly at small m. So the benefit of spatially 
# clustered sampling (i.e., consecutive ones) can increase as the system approaches criticality.

#At γ=0.5, ΔRm is still positive, but considerably smaller. Shared forcing contributes covariance between distant as well as nearby patches, 
#so random sampling loses less relative to consecutive sampling.

#At γ=2, the difference is almost zero.
#Here shared forcing dominates the spatial covariance, so sampling geometry becomes almost irrelevant.

#And at m=100=N, ΔR_100 =0 exactly, regardless of D,γ,λ, because both strategies observe the entire ring. 
#That's an excellent consistency check.

#There is another important pattern:
#The strategy advantage is largest at small/intermediate m and disappears as m→N.
#So we now have three interacting controls on whether monitoring geometry matters: sampling geometry matters most when m/N is small, D is appreciable, γ is small.
#And approaching criticality can amplify this effect.

#This gives you a very clean ecological interpretation: 
#clustered monitoring is most beneficial when spatial synchrony is generated locally through dispersal rather than 
#globally through shared environmental forcing.


#==========================
#One caution, though: our current Rm(λ) calculation is still a quasi-stationary calculation — we are inserting each instantaneous λ 
#into the stationary covariance solution. Since your actual simulation changes λ(t) from −0.5 to −0.05 over finite time, before calling this Rm(t) 
#in the paper, I think we should check whether the quasi-stationary approximation is sufficiently accurate for your forcing rate.
