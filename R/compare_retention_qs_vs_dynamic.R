# comparison between quasi-stationary calculation of Rm vs its dynamic analog


# call the function
source(here("R/dynamic_modal_variances.R"))
modal_dyn <- dynamic_modal_variances(
  N = 100,
  D = 0.1,
  lambda_start = -0.5,
  lambda_end = -0.05,
  transition_time = 200,
  dt = 0.05,
  sigma_local = 1,
  sigma_shared = 0.5
)

Rm_dyn_con <- dynamic_retention_consecutive(modal_dyn, m = 20)
Rm_dyn_rand <- dynamic_retention_random(modal_dyn, m = 20)

data_dyn <- data.frame(time = modal_dyn$time, lambda = modal_dyn$lambda, Consecutive = Rm_dyn_con,  Random = Rm_dyn_rand)
saveRDS(data_dyn, here("Results/data_dyn_retention.RDS"))
#=================================================================



data_dyn<-readRDS(here("Results/data_dyn_retention.RDS"))
# ============================================================
# Quasi-stationary retention

source(here("R/ring_retention.R"))
Rm_QS_con <- sapply(data_dyn$lambda, function(lambda_now){
    ring_retention(
      N = 100,
      m = 20,
      lambda = lambda_now,
      D = 0.1,
      sigma_local = 1,
      sigma_shared = 0.5,
      strategy = "consecutive"
    )
  }
)

Rm_QS_rand <- sapply(data_dyn$lambda, function(lambda_now){
    ring_retention(
      N = 100,
      m = 20,
      lambda = lambda_now,
      D = 0.1,
      sigma_local = 1,
      sigma_shared = 0.5,
      strategy = "random"
    )
  }
)

# ============================================================
# Combine dynamic and quasi-stationary results
# ============================================================

data_compare_Rm <- tibble(
  
  time = data_dyn$time,
  lambda = data_dyn$lambda,
  
  dyn_con = data_dyn$Consecutive,
  QS_con = Rm_QS_con,
  
  dyn_rand = data_dyn$Random,
  QS_rand = Rm_QS_rand
)

data_compare_long <- data_compare_Rm %>%
  
  pivot_longer(cols = c(dyn_con,QS_con,dyn_rand,QS_rand),
    names_to = "curve",
    values_to = "Rm"
  ) %>%
  mutate(strategy = case_when(
      grepl("con", curve) ~ "Consecutive",
      grepl("rand", curve) ~ "Random"
    ),
    method = case_when(
      grepl("dyn", curve) ~ "Dynamic",
      grepl("QS", curve) ~ "Quasi-stationary"
    )
  )


g1<-ggplot(data_compare_long, aes(x = lambda, y = Rm, color = strategy, linetype = method)) +
  geom_line(linewidth = 1) +
  labs(x = expression(lambda),
    y = expression(R[m]),
    color = "Sampling strategy",
    linetype = "Method") +
  theme_bw(base_size = 14, base_family = "sans") + 
  theme(legend.position = "top")

pdf(here("Results/compare_Rm_vs_lambda_qs_and_dynamic_approach.pdf"),  width = 9, height = 6)
print(g1)
dev.off()
#Figure X. Dynamic changes in spatial variance retention during approach to criticality. 
#Variance retention \(R_m\) as the local recovery rate \(\lambda\) increases from \(-0.5\) to \(-0.05\), 
#shown for a representative case with \(N=100\), \(m=20\), \(D=0.1\), and \(\gamma=\sigma_{\mathrm{shared}}/\sigma_{\mathrm{local}}=0.5\). 
#Solid lines show retention calculated from the dynamically evolving modal variances, while dotted lines show the corresponding 
#quasi-stationary prediction obtained from the stationary covariance at each instantaneous value of \(\lambda\). 
#Colours distinguish consecutive and random spatial sampling. Retention increases as the system approaches criticality 
#under both sampling strategies, indicating increasing preservation of local variance under spatial aggregation. 
#The quasi-stationary approximation increasingly overestimates dynamically realized retention near criticality because 
#the covariance structure requires finite time to adjust to the changing recovery rate. The dependence of this discrepancy 
#on dispersal, shared forcing, and aggregation scale is shown in Fig. Sx (from "R/summarise_qs_dynamic_lag.R").

data_compare_Rm <- data_compare_Rm %>%
  mutate(rel_lag_con = (QS_con - dyn_con) / QS_con,
         rel_lag_rand = (QS_rand - dyn_rand) / QS_rand)

ggplot(data_compare_Rm, aes(x = lambda)) +
  geom_line(aes(y = rel_lag_con, linetype = "Consecutive"),linewidth = 1) +
  geom_line(aes(y = rel_lag_rand, linetype = "Random"),linewidth = 1) +
  labs(x = expression(lambda),
    y = "Relative lag") +
  theme_bw(base_size = 14, base_family = "sans") +
  theme(legend.position = "top")
