# ============================================================
# Main analytical retention function
#
# strategy = "consecutive"
#          or "random"
#
# For random sampling, the returned value is the EXPECTED
# retention over equally likely m-patch subsets.

# N= patch numbers, must be >2
# m = CONSECUTIVE observed patches
# lambda < 0 = local recovery eigenvalue
# D = nearest-neighbour dispersal strength
# sigma_local = sigma_l
# sigma_shared = sigma_s
# ============================================================
library(here)
source(here("R/get_rho_ring.R"))
source(here("R/ring_retention_consecutive.R"))
source(here("R/ring_retention_random.R"))

ring_retention <- function(N, m, lambda, D, sigma_local = 1, sigma_shared = 0, strategy = c("consecutive", "random")){
  
  strategy <- match.arg(strategy)
  
  if (strategy == "consecutive"){
    out_c<-ring_retention_consecutive(N = N, m = m, lambda = lambda, D = D, sigma_local = sigma_local, sigma_shared = sigma_shared)
    return(out_c)
  }
  
  if (strategy == "random"){
    out_r<-ring_retention_random(N = N, m = m, lambda = lambda, D = D, sigma_local = sigma_local, sigma_shared = sigma_shared)
    return(out_r)
  }
}


# now call the function

# first, sanity check:
#ring_retention(N = 100, m = 1, lambda = -0.2, D = 0.1, strategy = "consecutive")# should be 1
#ring_retention(N = 100, m = 1, lambda = -0.2, D = 0.1, strategy = "random")# should be 1

# For D=0, sigma_s = 0, Rm= 1/m
#m_test <- 20
#R_con <- ring_retention(N = 100,  m = m_test,  lambda = -0.2,  D = 0,  sigma_local = 1,  sigma_shared = 0,  strategy = "consecutive")
#R_random <- ring_retention(N = 100,  m = m_test,  lambda = -0.2,  D = 0,  sigma_local = 1,  sigma_shared = 0,  strategy = "random")

#c(consecutive = R_con,
#  random = R_random,
#  theoretical_independent = 1 / m_test)

# For sigma_s >> sigma_l; Rm~1
#R_con <- ring_retention(N = 100,  m = m_test,  lambda = -0.2,  D = 0,  sigma_local = 1,  sigma_shared = 20,  strategy = "consecutive")
#R_random <- ring_retention(N = 100,  m = m_test,  lambda = -0.2,  D = 0,  sigma_local = 1,  sigma_shared = 20,  strategy = "random")
#R_con
#R_random

#------------------------------------------------------------------------------------------------------------------------------------------------
# here R_consecutive>R_random: when dispersal creates spatially localized synchrony, clustered/consecutive monitoring preserves more variance 
# than spatially random monitoring at the same sampling effort m.
#ring_retention(N = 100,  m = 20,  lambda = -0.2,  D = 0.1,  sigma_local = 1,  sigma_shared = 0.5,  strategy = "consecutive")
#ring_retention(N = 100,  m = 20,  lambda = -0.2,  D = 0.1,  sigma_local = 1,  sigma_shared = 0.5,  strategy = "random")

# here both are same
#ring_retention(N = 100,  m = 20,  lambda = -0.2,  D = 0,  sigma_local = 1,  sigma_shared = 0.5,  strategy = "consecutive")
#ring_retention(N = 100,  m = 20,  lambda = -0.2,  D = 0,  sigma_local = 1,  sigma_shared = 0.5,  strategy = "random")
#--------------------------------------------------------------------------------------------------------------------------------------






