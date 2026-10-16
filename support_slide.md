---
marp: true
theme: gaia
_class: lead
paginate: true
backgroundColor: #f5f5f5
color: #333
style: |
  section {
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
  }
  code {
    background: #eef2f5;
    padding: 0.1em 0.3em;
    border-radius: 4px;
    font-size: 0.85em;
  }
---

# Null hypothesis testing is just regression: Linear Modeling & Parameter Uncertainty in Genomics

**Speaker:** Vipin Kumar  

---

**Our Goal Today:** 
* 1. Prove that familiar tools ($t$-tests, rank-sum tests) are actually linear models.
* 2. Escape "hacky" normalizations by controlling confounders directly.
* 3. Replace point estimates with parameter simulations.

---

# Part 1: A toy model of genomics

- K562 Dnase and H3Kme4 ChIP data
- Is open chromatin affecting the magnitude of the H3K4me3 Fold change?

---

## $t$-Tests *Are* Linear Models

- A classical Student’s $t$-test is mathematically identical to a binary ordinary least squares (OLS) regression.

- **Takeaway:** You aren't learning new math; you are just removing the training wheels.

---

# Non-Parametric Bridge

- Because real genomic data introduces a variety of exotism (extreme values, ties, etc) let us simulate data


---

- **Simulation Grid Code:**
```r
dist_grid_data <- map_df(effect_sizes, function(delta) {
  tibble(true_effect = paste0("Delta = ", delta),
         control = list(rnorm(200, mean = 5, sd = 2)),
         treatment = list(rnorm(200, mean = 5 + delta, sd = 2)))
}) |> unnest(cols = c(control, treatment)) |> ...

```

---

## Asymptotic Equivalence: Wilcoxon vs. OLS on Ranks

- Comparing Mann-Whitney U statistics against OLS $t$-statistics on rank-transformed data.

---

- **Simulation :**
```r
simulation_results <- map_df(effect_sizes, function(delta) {
  map_df(1:30, function(sim) {
    wilcox_res <- wilcox.test(..., exact = FALSE)
    lm_res <- summary(lm(rank(value) ~ group, data = sim_data))
    ...
  })
})

```

---

- **Takeaway:** Rank regression yields the exact same statistical power as non-parametric tests while unlocking multi-variable confounder control.

---

# Part 2: Multiple Regression with Real Confounders

## The Naive Model vs. GC Content Bias

- **The perils of ignoring GC content:** Our naive model ignores sequence composition (gc_content), attributing technical capture bias entirely to biological activity.

---

- **Live-Coding Implementation:**

```r
# --- Naive model ---
fit_naive <- lm(active_ip_signal ~ dnase_peak, data = data_tbl)

# --- Controlled model: Adding GC content covariate ---
fit_controlled <- lm(
  active_ip_signal ~ dnase_peak + gc_content, 
  data = data_tbl
)

```

---

**Takeaway:** Multiple regression lets us control technical confounders in plain sight without mutating our target data matrix.


---

# Part 2: Residual Diagnostics

## The Funnel of Uncertainty & The Log Transformation

- Plotting residuals against fitted values reveals heteroscedasticity (megaphone funnel).

---


- **The Intervention:** Stabilizing variance via log transformation.


```r
data_tbl <- data_tbl |> mutate(logy = log(active_ip_signal))

fit_log_controlled <- lm(
  logy ~ dnase_peak + gc_content, 
  data = data_tbl
)

```

---

- **Diagnostic Check:** Residual scatter plots confirm the necessity of addressing variance scaling before downstream modeling.



---

# Part 3: The simulated worlds

## Parameter Uncertainty & Predictions

- Moving beyond fixed point estimates to sample from the joint parameter distribution (mvrnorm).

---


**Simulation Setup:**
```r
mu_hat <- coef(fit_log_controlled)
sigma_hat <- vcov(fit_log_controlled)

# Draw 1,000 parameter realizations from MLE variance-covariance matrix
param_draws <- mvrnorm(n = 1000, mu = mu_hat, Sigma = sigma_hat) %>% 
  as_tibble()

```


- **The Logic:** Every draw represents a plausible universe consistent with our empirical observations.



---

# Part 3: Counterfactual Visualization

## Empirical Violins Meet Model-Implied Counterfactuals

---

- **Code Implementation:**
```r
sim_means_df <- map_df(1:nrow(param_draws), function(i) {
  # Compute group predictions per draw
  ...
})

ggplot() +
  geom_violin(data = data_tbl, aes(x = dnase_label, y = logy), ...) +
  geom_line(data = sim_means_df, aes(x = group, y = mean_prediction, group = sample_id, color = diff), ...)

```

---

**The Payoff:** We visualize the exact distribution of plausible differences instead of relying on an arbitrary $p$-value threshold.

---

# Conclusion: Summary & Next Steps

## The Modern Genomic Analyst's Toolkit

**What We Achieved Today:**

* 1. **Demystified Statistics:** Proved $t$-tests and rank tests are special cases of linear models.


* 2. **Clean Confounder Control:** Replaced black-box preprocessing hacks with explicit regression covariates.


* 3. **Counterfactual Simulation:** Used multivariate normal sampling to quantify parameter uncertainty.

---

**Thank you!** 

Q&A.



