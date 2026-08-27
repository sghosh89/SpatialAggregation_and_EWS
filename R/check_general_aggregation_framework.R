library(here)
library(tidyverse)
library(showtext)


# ------------------------------------------------------------
# Theoretical retention factor under spatial aggregation
# ------------------------------------------------------------

# m   = number of local patches included in the spatial aggregate
# rho = pairwise correlation (spatial synchrony) among patches
#
# R_m measures the fraction of local-scale variance retained
# after averaging across m patches.
#
# When rho = 0, R_m = 1/m:
# independent local fluctuations are increasingly averaged out.
#
# When rho = 1, R_m = 1:
# perfectly synchronous fluctuations are fully retained.
#
# For intermediate rho, R_m declines with aggregation
# and approaches rho as m becomes large.

get_retention_factor <- function(m, rho) {
  out <- (1 / m) + (((m - 1) / m) * rho)
  return(out)
}


# ------------------------------------------------------------
# Define spatial system and parameter values
# ------------------------------------------------------------

# Total number of patches in the full spatial system
N <- 100

# Aggregation scales, from one patch to the entire system
m_values <- 1:N

# Pairwise spatial correlation values to compare
rho_values <- c(0, 0.25, 0.50, 0.75, 1)

# ------------------------------------------------------------
# Create dataframe for plotting
# ------------------------------------------------------------

# Generate all combinations of aggregation scale m and spatial correlation rho
retention_df <- expand.grid(m   = m_values, 
                            rho = rho_values)

# Calculate relative aggregation scale.
#
# m/N = 1/N corresponds to a single local patch,
# whereas m/N = 1 corresponds to averaging across the entire spatial system.
retention_df$m_relative <- retention_df$m / N

# Calculate the theoretical retention factor for every combination of m and rho
retention_df$Rm <- get_retention_factor(m   = retention_df$m, rho = retention_df$rho)

# Keep rho as a numeric variable in the dataframe,
# but create a separate factor for plotting and legend labels
retention_df$rho_plot <- factor(retention_df$rho, levels = rho_values, labels = paste0("\u03c1 = ", rho_values))


# ------------------------------------------------------------
# Plot theoretical retention curves
# ------------------------------------------------------------

retention_plot <- ggplot(retention_df, aes(x = m_relative, y = Rm, colour = rho_plot, linetype = rho_plot)) +
  
  # Plot one theoretical retention curve for each rho
  geom_line(linewidth = 1) +
  
  # Axis and legend labels
  labs(
    x = "Relative aggregation scale, m/N",
    y = expression(paste("Retention factor, ", R[m])),
    colour = expression(rho),
    linetype = expression(rho)
  ) +
  
  # Display the full range of aggregation scales
  coord_cartesian(
    xlim = c(1 / N, 1),
    ylim = c(0, 1)
  ) +
  
  # Clean manuscript-style theme
  theme_bw(base_size = 15) +
  
  # Place the single combined legend above the figure
  theme(legend.title = element_blank(),
    legend.position = "top",
    #legend.position.inside = c(0.80, 0.85),
    legend.background = element_rect(
      fill = "white",
      colour = NA#"black"
    )
  )

# Display figure
retention_plot

pdf(here("Results/plot_general_aggregation_framework.pdf"),
    width = 7, height = 6)
print(retention_plot)
dev.off()

#######################

