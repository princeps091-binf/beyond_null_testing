library(tidyverse)
library(ggplot2)
library(MASS)

theme_set(theme_minimal(base_size = 14))
# %%
encode_df <- read_delim('./data/encode_k562_promoters.csv',delim=',')

data_tbl <- encode_df |> dplyr::filter(active_ip_signal > 0 )
# %%
g <- data_tbl |> 
	ggplot(aes(x = as.factor(dnase_peak),y = active_ip_signal,fill=as.factor(dnase_peak))) +
	geom_boxplot()
g


# %%
# ==============================================================================
# LIVE-CODING DEMO: From t-tests to Linear Models (Real ENCODE Data)
# ==============================================================================


# ==============================================================================
# PART 1: MINIMAL PROOF OF CONCEPT (Using Real Subsets)
# ==============================================================================
# Extract two real vectors based on ATAC peak status
open_promoters   <- data_tbl %>% filter(dnase_peak == 1) %>% pull(active_ip_signal)
closed_promoters <- data_tbl %>% filter(dnase_peak == 0) %>% pull(active_ip_signal)

# 1. Classical Two-Sample t-test passing two vectors
t_res <- t.test(open_promoters, closed_promoters,var.equal = TRUE)

# 2. Equivalent Linear Model (OLS) using the full tidy table
lm_res <- lm(active_ip_signal ~ dnase_peak, data = data_tbl)

# Verify identity
t_res$statistic
summary(lm_res)$coefficients["dnase_peak", "t value"]

# %%

set.seed(42)

# 1. Define a range of effect sizes (mean shifts) to test
effect_sizes <- seq(0, 3, by = 0.2)
# Generate distribution data across selected deltas
dist_grid_data <- map_df(effect_sizes, function(delta) {
  tibble(
    true_effect = paste0("Delta = ", delta),
    delta_val = delta,
    control = list(rnorm(200, mean = 5, sd = 2)),
    treatment = list(rnorm(200, mean = 5 + delta, sd = 2))
  )
}) |>
  unnest(cols = c(control, treatment)) |>
  pivot_longer(
    cols = c(control, treatment),
    names_to = "group",
    values_to = "value"
  ) |>
  mutate(
    group = str_to_title(group),
    # Ensure factor order for clean plotting
    true_effect = factor(true_effect, levels = paste0("Delta = ", effect_sizes))
  )

# Plotting the density distributions using facet_grid
g <- ggplot(dist_grid_data, aes(x = value, fill = group, color = group)) +
  geom_density(alpha = 0.4, linewidth = 0.8) +
  facet_grid(. ~ true_effect) +
  scale_fill_manual(values = c("Control" = "#e76f51", "Treatment" = "#2a9d8f")) +
  scale_color_manual(values = c("Control" = "#e76f51", "Treatment" = "#2a9d8f")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution Separation Across Effect Sizes",
    subtitle = "As delta increases, overlap shrinks, driving both the Mann-Whitney U test and OLS rank regression",
    x = "Simulated Signal Value",
    y = "Density",
    fill = "Group",
    color = "Group"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12)
  )
g
# %%

# 2. Run simulations across all effect sizes
simulation_results <- map_df(effect_sizes, function(delta) {
  # Repeat simulation a few times per effect size for smooth tracking
  map_df(1:30, function(sim) {
    sim_data <- tibble(
      group = rep(c("Control", "Treatment"), each = 50),
      value = c(rnorm(50, mean = 5, sd = 2), rnorm(50, mean = 5 + delta, sd = 2))
    )
    
    # Run Wilcoxon rank-sum test
    wilcox_res <- wilcox.test(sim_data|>filter(group == "Treatment")|>pull(value),sim_data|>filter(group == "Control")|>pull(value), exact = FALSE)
    
    # Run OLS Linear Model on rank-transformed data
    sim_data <- sim_data |> mutate(ranked_val = rank(value))
    lm_res <- summary(lm(ranked_val ~ group, data = sim_data))
    
    tibble(
      true_effect = delta,
      sim_id = sim,
      wilcox_stat = wilcox_res$statistic,
      wilcox_p = wilcox_res$p.value,
      lm_t_stat = lm_res$coefficients["groupTreatment", "t value"],
      lm_p = lm_res$coefficients["groupTreatment", "Pr(>|t|)"]
    )
  })
})

# %%
# Reshape data into long format for faceting
sim_long <- simulation_results |>
  pivot_longer(
    cols = c(lm_t_stat, wilcox_stat),
    names_to = "test_type",
    values_to = "statistic_value"
  ) |>
  mutate(
    test_type = recode(
      test_type,
      "lm_t_stat" = "OLS on Ranks (t-statistic)",
      "wilcox_stat" = "Mann-Whitney U (W Statistic)"
    )
  )

# Plot using facet_wrap with free y-scales
ggplot(sim_long, aes(x = true_effect, y = statistic_value, color = test_type)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "loess", color = "#1d3557", se = FALSE, linetype = "dashed") +
  facet_wrap(~ test_type, scales = "free_y") +
  scale_color_manual(values = c("OLS on Ranks (t-statistic)" = "#2b5c8f", 
                                "Mann-Whitney U (W Statistic)" = "#2a9d8f")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Asymptotic Evolution of Test Statistics",
    subtitle = "Comparing OLS on Ranks against the Mann-Whitney U test across effect sizes",
    x = "True Effect Size (delta)",
    y = "Test Statistic Value"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank(),
    legend.position = "none" # Hide legend since facets are already labeled
  )
# %%
# 3. Visualize the tracking relationship
ggplot(simulation_results, aes(x = lm_t_stat, y = as.numeric(wilcox_stat))) +
  geom_point(alpha = 0.6, color = "#2b5c8f", size = 2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Asymptotic Equivalence: Wilcoxon vs. OLS on Ranks",
    subtitle = "Tracking test statistics across varying effect sizes (small to large differences)",
    x = "OLS Model t-statistic (on Ranks)",
    y = "Wilcoxon Rank-Sum Test Statistic"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank()
  )
# %%
# ==============================================================================
# PART 2: MULTIPLE REGRESSION WITH REAL CONTROL TRACKS & CONFOUNDERS
# ==============================================================================

# --- Naive model ---
fit_naive <- lm(active_ip_signal ~ dnase_peak, data = data_tbl)
summary(fit_naive)$coefficients

# %%

g <- data_tbl |> 
	ggplot(aes(x = active_ip_signal,y = gc_content,color=as.factor(dnase_peak))) +
	geom_point()+
  scale_x_log10() 
g

# %%

# --- Controlled model: Using real control track & GC content as covariates ---
fit_controlled <- lm(
  active_ip_signal ~ dnase_peak + gc_content, 
  data = data_tbl
)
summary(fit_controlled)$coefficients


# %%

cor(residuals(fit_naive), data_tbl$gc_content)
cor(residuals(fit_controlled), data_tbl$gc_content)

# %%

data_tbl|>mutate(ctrl_residuals = fit_controlled$residuals,ctrl_pred = fit_controlled$fitted.values,naive_residuals = fit_naive$residuals, naive_pred = fit_naive$fitted.values)|>ggplot(aes(x=gc_content,y=ctrl_residuals,color=dnase_peak))+geom_point()

# %%

fit_interaction <- lm(
  active_ip_signal ~ dnase_peak + gc_content + dnase_peak * gc_content, 
  data = data_tbl
)
summary(fit_interaction)$coefficients

# %%

data_tbl|>mutate(inter_residuals = fit_interaction$residuals,inter_pred = fit_interaction$fitted.values,naive_residuals = fit_naive$residuals, naive_pred = fit_naive$fitted.values)|>ggplot(aes(x=gc_content,y=inter_residuals,color=dnase_peak))+geom_point()
# %%
data_tbl|>mutate(residuals = fit_controlled$residuals,pred = fit_controlled$fitted.values)|>ggplot(aes(x=pred,y=residuals))+geom_point()

# %%

data_tbl <- data_tbl |> mutate(logy = log(active_ip_signal))


fit_log_controlled <- lm(
  logy ~ dnase_peak + gc_content + dnase_peak * gc_content, 
  data = data_tbl
)
summary(fit_log_controlled)$coefficients
# %%

data_tbl|>mutate(residuals = fit_log_controlled$residuals,pred = fit_log_controlled$fitted.values)|>ggplot(aes(x=pred,y=residuals))+geom_point()

# %%

data_tbl|>mutate(log_residuals = fit_log_controlled$residuals,log_pred = fit_log_controlled$fitted.values,naive_residuals = fit_naive$residuals, naive_pred = fit_naive$fitted.values)|>ggplot(aes(x=gc_content,y=log_residuals,color=dnase_peak))+geom_point()

# %%

data_tbl <- data_tbl |> mutate(gc_centered = gc_content - mean(gc_content, na.rm = TRUE))

fit_centered <- lm(logy ~ dnase_peak + gc_centered + dnase_peak * gc_centered, data = data_tbl)
summary(fit_centered)

# %%


data_tbl|>mutate(residuals = fit_centered$residuals,pred = fit_centered$fitted.values)|>ggplot(aes(x=gc_centered,y=residuals))+geom_point()

# ==============================================================================
# PART 3: MCELREATH WORKFLOW — PARAMETER UNCERTAINTY & PREDICTIONS
# ==============================================================================

mu_hat <- coef(fit_centered)
sigma_hat <- vcov(fit_centered)

# Draw 1,000 parameter realizations from real MLE variance-covariance matrix
param_draws <- mvrnorm(n = 1000, mu = mu_hat, Sigma = sigma_hat) %>% as_tibble()

# %%

# 2. Compute group-wise mean predictions per draw using the actual GC distribution
# 2. Compute group-wise mean predictions per draw using the actual GC distribution & interaction
sim_means_df <- map_df(1:nrow(param_draws), function(i) {
  b0         <- as.numeric(param_draws[i, "(Intercept)"])
  b_dnase    <- as.numeric(param_draws[i, "dnase_peak"])
  b_gc       <- as.numeric(param_draws[i, "gc_centered"])
  b_inter    <- as.numeric(param_draws[i, "dnase_peak:gc_centered"])
  
  # Predict across all actual promoter rows for this parameter draw
  # Formula: b0 + b_dnase*X + b_gc*GC + b_inter*(X * GC)
  pred_y <- b0 + (b_dnase * data_tbl$dnase_peak) + 
                (b_gc * data_tbl$gc_centered) + 
                (b_inter * data_tbl$dnase_peak * data_tbl$gc_centered)
  
  # Calculate empirical group means using the actual covariate distribution
  mean_closed <- mean(pred_y[data_tbl$dnase_peak == 0], na.rm = TRUE)
  mean_open   <- mean(pred_y[data_tbl$dnase_peak == 1], na.rm = TRUE)
  
  tibble(
    sample_id = i,
    "Closed" = mean_closed,
    "Open" = mean_open,
    diff = as.numeric(mean_open - mean_closed)
  )
}) %>%  
  pivot_longer(
    cols = c("Closed", "Open"),
    names_to = "group",
    values_to = "mean_prediction"
  )

# %%
sim_means_df$group <- factor(sim_means_df$group, levels = c("Closed", "Open"))
data_tbl$dnase_label <- factor(ifelse(data_tbl$dnase_peak == 1, "Open", "Closed"), levels = c("Closed", "Open"))

# 3. Build the visualization combining raw boxplots and connected derived means
ggplot() +
  # Layer 1: Raw data violins/boxplots in the background
  geom_violin(
    data = data_tbl, 
    aes(x = dnase_label, y = logy), 
    fill = "gray95", color = "gray50", alpha = 0.5, 
    width = 0.35
  ) +
  # Layer 2: Model-implied lines connecting derived group means, colored by effect size
  geom_line(
    data = sim_means_df, 
    aes(x = group, y = mean_prediction, group = sample_id, color = diff), 
    alpha = 0.4, linewidth = 0.8
  ) +
  # Layer 3: Endpoints of the lines representing the group means for each draw
  geom_point(
    data = sim_means_df, 
    aes(x = group, y = mean_prediction, color = diff), 
    alpha = 0.5, size = 1.8
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Empirical Distribution Meets Interaction-Aware Counterfactuals",
    subtitle = "Lines connect derived group-wise means accounting for GC-centering and interaction across 1000 parameter draws",
    x = "Chromatin State",
    y = "log(H3K4me3 Fold-Change)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
