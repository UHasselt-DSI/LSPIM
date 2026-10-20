## Processing new simulations

## RI model
getwd()
load("Revision/Sims_RI.Rdata")


test_runs <- dplyr::bind_rows(test_store)
contrast_runs <- dplyr::bind_rows(contrast_store)

table(test_runs$N)
table(contrast_runs$N)

model_levels <- c("GEE", "PGEE", "PIM", "PIM_adj")

procedure_levels <- c(
  "Asymptotic_Global", "Asymptotic_Holm",
  "Finite_Global", "Finite_Holm"
)

# ------------------------------------------------------------
# 1. Variance comparison
# ------------------------------------------------------------

contrast_runs <- bind_rows(contrast_store) |>
  mutate(
    valid = is.finite(estimate) &
      is.finite(variance) &
      variance > 0
  )

# Detect accidental duplicate storage
stopifnot(
  !anyDuplicated(
    contrast_runs[c("N", "run", "model", "hypothesis", "contrast")]
  )
)

variance_comparison <- contrast_runs |>
  group_by(N, model, hypothesis, contrast) |>
  summarise(
    n_total = n(),
    n_valid = sum(valid),
    n_failed = sum(!valid),
    
    mean_estimate = if (any(valid)) {
      mean(estimate[valid])
    } else {
      NA_real_
    },
    
    empirical_variance = if (sum(valid) >= 2L) {
      var(estimate[valid])
    } else {
      NA_real_
    },
    
    average_estimated_variance = if (any(valid)) {
      mean(variance[valid])
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) |>
  mutate(
    failure_rate = n_failed / n_total,
    empirical_SD = sqrt(empirical_variance),
    estimated_SD = sqrt(average_estimated_variance),
    variance_ratio =
      empirical_variance / average_estimated_variance
  )

print(variance_comparison, n = Inf)

# ------------------------------------------------------------
# 2. Rejection rates
# ------------------------------------------------------------

test_runs <- bind_rows(test_store)

stopifnot(
  !anyDuplicated(
    test_runs[c("N", "run", "model", "hypothesis", "procedure")]
  )
)

rejection_summary <- test_runs |>
  group_by(N, model, hypothesis, procedure) |>
  summarise(
    n_total = n(),
    n_valid = sum(!is.na(reject)),
    n_failed = sum(is.na(reject)),
    n_rejected = sum(reject, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    failure_rate = n_failed / n_total,
    rejection_rate = n_rejected / na_if(n_valid, 0L),
    
    MCSE = sqrt(
      rejection_rate * (1 - rejection_rate) /
        na_if(n_valid, 0L)
    )
  )

# Wilson intervals: also work when there are zero rejections
z <- qnorm(0.975)

rejection_summary <- rejection_summary |>
  mutate(
    denominator = na_if(n_valid, 0L),
    interval_center =
      (rejection_rate + z^2 / (2 * denominator)) /
      (1 + z^2 / denominator),
    
    interval_halfwidth =
      z * sqrt(
        rejection_rate * (1 - rejection_rate) / denominator +
          z^2 / (4 * denominator^2)
      ) / (1 + z^2 / denominator),
    
    lower = pmax(0, interval_center - interval_halfwidth),
    upper = pmin(1, interval_center + interval_halfwidth)
  ) |>
  dplyr::select(-denominator, -interval_center, -interval_halfwidth)

print(rejection_summary, n = Inf)

# ------------------------------------------------------------
# 3. Variance-ratio plot
# ------------------------------------------------------------

variance_plot_data <- variance_comparison |>
  mutate(
    model = factor(model, levels = model_levels),
    hypothesis = factor(
      hypothesis,
      levels = c("time", "treatment"),
      labels = c("Equal temporal effects", "Treatment-effect constancy")
    ),
    contrast_label = paste("Contrast", contrast)
  )

p_variance <- ggplot(
  variance_plot_data,
  aes(
    x = N,
    y = variance_ratio,
    shape = model,
    linetype = model,
    group = model
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_line() +
  geom_point(size = 2.5) +
  facet_wrap(vars(hypothesis, contrast_label), ncol = 2) +
  scale_shape_manual(
    values = c(GEE = 16, PGEE = 17, PIM = 1, PIM_adj = 2)
  ) +
  scale_linetype_manual(
    values = c(
      GEE = "solid",
      PGEE = "dashed",
      PIM = "dotted",
      PIM_adj = "dotdash"
    )
  ) +
  scale_x_continuous(breaks = sample_sizes) +
  labs(
    x = "Number of subjects",
    y = "Empirical variance / mean estimated variance",
    shape = "Model",
    linetype = "Model"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_variance)

# ------------------------------------------------------------
# 4. Rejection-rate plot
# ------------------------------------------------------------

rejection_plot_data <- rejection_summary |>
  mutate(
    model = factor(model, levels = model_levels),
    procedure = factor(procedure, levels = procedure_levels)
  )

procedure_labels <- c(
  Asymptotic_Global = "Global: chi-square",
  Asymptotic_Holm   = "Holm: z",
  Finite_Global    = "Global: F",
  Finite_Holm      = "Holm: t"
)

pd <- position_dodge(width = 0.8)

p_rejection <- ggplot(
  rejection_plot_data,
  aes(
    x = N,
    y = rejection_rate,
    shape = procedure,
    linetype = procedure,
    group = procedure
  )
) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_line(position = pd) +
  geom_point(position = pd, size = 2) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = pd,
    width = 0.5
  ) +
  facet_grid(hypothesis ~ model) +
  scale_shape_manual(
    values = setNames(c(16, 17, 1, 2), procedure_levels),
    labels = procedure_labels
  ) +
  scale_linetype_manual(
    values = setNames(
      c("solid", "dashed", "dotted", "dotdash"),
      procedure_levels
    ),
    labels = procedure_labels
  ) +
  scale_x_continuous(breaks = sample_sizes) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    x = "Number of subjects",
    y = "Rejection rate",
    shape = "Procedure",
    linetype = "Procedure",
    caption = "Error bars: 95% Wilson Monte Carlo intervals."
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_rejection)

# ------------------------------------------------------------
# 5. Test-failure plot
# ------------------------------------------------------------

p_failure <- ggplot(
  rejection_plot_data,
  aes(
    x = N,
    y = failure_rate,
    shape = procedure,
    linetype = procedure,
    group = procedure
  )
) +
  geom_line() +
  geom_point(size = 2) +
  facet_grid(hypothesis ~ model) +
  scale_shape_manual(
    values = setNames(c(16, 17, 1, 2), procedure_levels),
    labels = procedure_labels
  ) +
  scale_linetype_manual(
    values = setNames(
      c("solid", "dashed", "dotted", "dotdash"),
      procedure_levels
    ),
    labels = procedure_labels
  ) +
  scale_x_continuous(breaks = sample_sizes) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    x = "Number of subjects",
    y = "Test failure rate",
    shape = "Procedure",
    linetype = "Procedure"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p_failure)




library(dplyr)
library(ggplot2)

M <- nr_of_iters
alpha <- 0.05

mc_lower <- qbinom(0.025, M, alpha) / M
mc_upper <- qbinom(0.975, M, alpha) / M

plot_data <- rejection_summary |>
  mutate(
    model = factor(
      model, levels = c("GEE", "PGEE", "PIM", "PIM_adj")
    ),
    procedure = factor(
      procedure,
      levels = c(
        "Asymptotic_Global", "Asymptotic_Holm",
        "Finite_Global", "Finite_Holm"
      ),
      labels = c(
        "Global: chi-square", "Holm: z",
        "Global: F", "Holm: t"
      )
    )
  )

p_rejection <- ggplot(
  plot_data,
  aes(
    x = N, y = rejection_rate,
    shape = model, linetype = model, group = model
  )
) +
  geom_hline(
    yintercept = alpha,
    linetype = "dashed", colour = "grey40"
  ) +
  geom_hline(
    yintercept = c(mc_lower, mc_upper),
    linetype = "dotted", colour = "grey50"
  ) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5) +
  facet_wrap(~ procedure, ncol = 2) +
  scale_shape_manual(
    values = c(GEE = 16, PGEE = 17, PIM = 1, PIM_adj = 2)
  ) +
  scale_linetype_manual(
    values = c(
      GEE = "solid", PGEE = "dashed",
      PIM = "dotted", PIM_adj = "dotdash"
    )
  ) +
  scale_x_continuous(
    breaks = sort(unique(plot_data$N)),
    guide = guide_axis(angle = 45)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, NA)
  ) +
  labs(
    x = "Number of subjects",
    y = "Rejection rate",
    shape = "Model",
    linetype = "Model",
    caption = paste0(
      "Dashed: 5%. Dotted: pointwise 95% Monte Carlo range; M = ", M, "."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p_treatment <- p_rejection %+%
  filter(plot_data, hypothesis == "treatment") +
  labs(title = "Treatment-effect constancy")

p_time <- p_rejection %+%
  filter(plot_data, hypothesis == "time") +
  labs(title = "Equal temporal effects")

print(p_treatment)
print(p_time)






library(dplyr)

parameter_runs <- bind_rows(
  lapply(model_store, function(x) {
    
    b <- x$coefficients
    V <- x$vcov
    
    parameter_names <- names(b)
    if (is.null(parameter_names)) {
      parameter_names <- paste0("beta", seq_along(b))
    }
    
    stopifnot(
      nrow(V) == length(b),
      ncol(V) == length(b)
    )
    
    data.frame(
      N = x$N,
      run = x$run,
      model = x$model,
      parameter = seq_along(b),
      parameter_name = parameter_names,
      estimate = as.numeric(b),
      model_variance = diag(V),
      row.names = NULL
    )
  })
) |>
  mutate(
    valid = is.finite(estimate) &
      is.finite(model_variance) &
      model_variance > 0
  )

parameter_variance_comparison <- parameter_runs |>
  group_by(N, model, parameter, parameter_name) |>
  summarise(
    n_total = n(),
    n_valid = sum(valid),
    n_failed = sum(!valid),
    
    mean_estimate = if (any(valid)) {
      mean(estimate[valid])
    } else {
      NA_real_
    },
    
    empirical_variance = if (sum(valid) >= 2L) {
      var(estimate[valid])
    } else {
      NA_real_
    },
    
    average_model_variance = if (any(valid)) {
      mean(model_variance[valid])
    } else {
      NA_real_
    },
    
    .groups = "drop"
  ) |>
  mutate(
    failure_rate = n_failed / n_total,
    empirical_SD = sqrt(empirical_variance),
    model_SD = sqrt(average_model_variance),
    variance_ratio = empirical_variance / average_model_variance
  ) |>
  arrange(N, model, parameter)

print(parameter_variance_comparison, n = Inf)




library(ggplot2)

ggplot(
  parameter_variance_comparison,
  aes(
    x = N,
    y = variance_ratio,
    shape = model,
    linetype = model,
    group = model
  )
) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line() +
  geom_point(size = 2.5) +
  facet_wrap(~ parameter, ncol = 2, labeller = label_both) +
  scale_shape_manual(
    values = c(GEE = 16, PGEE = 17, PIM = 1, PIM_adj = 2)
  ) +
  scale_linetype_manual(
    values = c(
      GEE = "solid", PGEE = "dashed",
      PIM = "dotted", PIM_adj = "dotdash"
    )
  ) +
  scale_x_continuous(
    breaks = sort(unique(parameter_variance_comparison$N)),
    guide = guide_axis(angle = 45)
  ) +
  labs(
    x = "Number of subjects",
    y = "Empirical variance / mean model variance",
    shape = "Model",
    linetype = "Model"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")


benchmark_runs <- dplyr::bind_rows(benchmark_store)

benchmark_summary <- benchmark_runs |>
  dplyr::group_by(N, hypothesis, procedure) |>
  dplyr::summarise(
    n_total = dplyr::n(),
    n_valid = sum(!is.na(reject)),
    n_failed = sum(is.na(reject)),
    n_rejected = sum(reject, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    failure_rate = n_failed / n_total,
    rejection_rate = n_rejected / dplyr::na_if(n_valid, 0L),
    MCSE = sqrt(
      rejection_rate * (1 - rejection_rate) /
        dplyr::na_if(n_valid, 0L)
    )
  )

print(benchmark_summary, n = Inf)



## Appendix benchmark plots

library(dplyr)
library(ggplot2)

M <- nr_of_iters
alpha <- 0.05

mc_lower <- qbinom(0.025, M, alpha) / M
mc_upper <- qbinom(0.975, M, alpha) / M

# PIM treatment-constancy results
pim_plot_data <- rejection_summary |>
  filter(hypothesis == "treatment") |>
  transmute(
    N,
    rejection_rate,
    method = as.character(model),
    panel = recode(
      as.character(procedure),
      Asymptotic_Global = "Global: chi-square",
      Asymptotic_Holm   = "Holm: z",
      Finite_Global    = "Global: F",
      Finite_Holm      = "Holm: t"
    )
  )

# nparLD interaction results
benchmark_plot_data <- benchmark_summary |>
  transmute(
    N,
    rejection_rate,
    method = recode(
      as.character(procedure),
      nparLD_WTS = "nparLD WTS",
      nparLD_ATS = "nparLD ATS"
    ),
    panel = "nparLD: RTE interaction"
  )

appendix_plot_data <- bind_rows(
  pim_plot_data,
  benchmark_plot_data
) |>
  mutate(
    panel = factor(
      panel,
      levels = c(
        "Global: chi-square", "Holm: z",
        "Global: F", "Holm: t",
        "nparLD: RTE interaction"
      )
    )
  )

p_appendix <- ggplot(
  appendix_plot_data,
  aes(
    x = N,
    y = rejection_rate,
    shape = method,
    linetype = method,
    group = method
  )
) +
  geom_hline(
    yintercept = alpha,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_hline(
    yintercept = c(mc_lower, mc_upper),
    linetype = "dotted",
    colour = "grey50"
  ) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5) +
  facet_wrap(~ panel, ncol = 2) +
  scale_shape_manual(
    values = c(
      GEE = 16, PGEE = 17, PIM = 1, PIM_adj = 2,
      "nparLD WTS" = 15, "nparLD ATS" = 0
    )
  ) +
  scale_linetype_manual(
    values = c(
      GEE = "solid",
      PGEE = "dashed",
      PIM = "dotted",
      PIM_adj = "dotdash",
      "nparLD WTS" = "solid",
      "nparLD ATS" = "dashed"
    )
  ) +
  scale_x_continuous(
    breaks = sort(unique(appendix_plot_data$N)),
    guide = guide_axis(angle = 45)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, NA)
  ) +
  labs(
    x = "Number of subjects",
    y = "Rejection rate",
    shape = "Method",
    linetype = "Method",
    caption = paste0(
      "PIM panels: treatment-effect constancy; nparLD panel: RTE interaction.\n",
      "Dashed: 5%. Dotted: pointwise 95% binomial reference range; M = ",
      M, "."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  guides(
    shape = guide_legend(nrow = 2),
    linetype = guide_legend(nrow = 2)
  )

print(p_appendix)




library(dplyr)
library(ggplot2)

# These must match the degrees of freedom used in your tests
df_trt  <- function(N) N - 1
df_time <- function(N) N - 1

density_data <- bind_rows(contrast_store) |>
  filter(
    is.finite(estimate),
    is.finite(variance),
    variance > 0
  ) |>
  mutate(
    statistic = estimate / sqrt(variance),
    nu = if_else(
      hypothesis == "treatment",
      df_trt(N),
      df_time(N)
    ),
    contrast_label = paste("Contrast", contrast)
  )

# Reference densities for each panel
reference_data <- density_data |>
  distinct(N, hypothesis, contrast_label, nu) |>
  group_by(N, hypothesis, contrast_label, nu) |>
  reframe(
    x = seq(-5, 5, length.out = 501),
    t_density = dt(x, df = first(nu)),
    normal_density = dnorm(x)
  )

plot_contrast_density <- function(hypothesis_name) {
  ggplot(
    filter(density_data, hypothesis == hypothesis_name),
    aes(x = statistic, colour = model)
  ) +
    geom_density(linewidth = 0.7) +
    geom_line(
      data = filter(reference_data, hypothesis == hypothesis_name),
      aes(x = x, y = t_density),
      inherit.aes = FALSE,
      colour = "black",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    geom_line(
      data = filter(reference_data, hypothesis == hypothesis_name),
      aes(x = x, y = normal_density),
      inherit.aes = FALSE,
      colour = "black",
      linetype = "dotted",
      linewidth = 0.8
    ) +
    facet_grid(contrast_label ~ N, labeller = label_both) +
    coord_cartesian(xlim = c(-5, 5)) +
    labs(
      x = "Estimated contrast / estimated SE",
      y = "Density",
      colour = "Model",
      caption = "Dashed: t reference; dotted: standard normal reference."
    ) +
    theme_bw()
}

plot_contrast_density("treatment")
plot_contrast_density("time")


# Add Lachin benchmark

library(dplyr)
library(ggplot2)

# Lachin: one result per dataset
lachin_summary <- bind_rows(lachin_store) |>
  group_by(N) |>
  summarise(
    n_valid = sum(!is.na(reject)),
    n_failed = sum(is.na(reject)),
    n_rejected = sum(reject, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    rejection_rate = n_rejected / na_if(n_valid, 0L)
  )

# Repeat the same Lachin benchmark in each treatment procedure panel
procedures <- rejection_summary |>
  distinct(procedure)

lachin_overlay <- cross_join(lachin_summary, procedures) |>
  mutate(
    hypothesis = "treatment",
    model = "Lachin"
  )

plot_data <- bind_rows(
  rejection_summary |>
    filter(hypothesis %in% c("treatment", "time")) |>
    dplyr::select(N, hypothesis, procedure, model, rejection_rate),
  
  lachin_overlay |>
    dplyr::select(N, hypothesis, procedure, model, rejection_rate)
) |>
  mutate(
    hypothesis = factor(
      hypothesis,
      levels = c("treatment", "time"),
      labels = c(
        "Treatment-effect constancy",
        "Equal temporal effects"
      )
    ),
    model = factor(
      model,
      levels = c("GEE", "PGEE", "PIM", "PIM_adj", "Lachin")
    )
  )

M <- nr_of_iters
mc_lower <- qbinom(0.025, M, 0.05) / M
mc_upper <- qbinom(0.975, M, 0.05) / M

ggplot(
  plot_data,
  aes(
    x = N, y = rejection_rate,
    group = model, shape = model, linetype = model
  )
) +
  geom_line(colour = "black", linewidth = 0.6) +
  geom_point(colour = "black", size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_hline(
    yintercept = c(mc_lower, mc_upper),
    linetype = "dotted"
  ) +
  scale_shape_manual(values = c(
    GEE = 16, PGEE = 17, PIM = 1, PIM_adj = 2, Lachin = 8
  )) +
  scale_linetype_manual(values = c(
    GEE = "solid",
    PGEE = "dashed",
    PIM = "dotted",
    PIM_adj = "dotdash",
    Lachin = "longdash"
  )) +
  scale_x_continuous(breaks = sort(unique(plot_data$N))) +
  facet_grid(hypothesis ~ procedure) +
  labs(
    x = "Number of subjects",
    y = "Rejection rate",
    shape = "Method",
    linetype = "Method",
    caption = paste(
      "Lachin's treatment-constancy benchmark is repeated unchanged",
      "across the treatment panels.",
      "Dotted horizontal lines: 95% Monte Carlo reference bounds at 5%."
    )
  ) +
  theme_bw() +
  theme(legend.position = "bottom")







library(dplyr)
library(tidyr)

test_runs <- bind_rows(test_store)

expected <- expand_grid(
  N = sample_sizes,
  run = seq_len(nr_of_iters),
  model = c("PIM", "PIM_adj", "GEE", "PGEE"),
  hypothesis = c("treatment", "time"),
  procedure = c(
    "Asymptotic_Global", "Asymptotic_Holm",
    "Finite_Global", "Finite_Holm"
  )
)

keys <- c("N", "run", "model", "hypothesis", "procedure")

missing_tests <- anti_join(expected, test_runs, by = keys)

duplicate_tests <- test_runs |>
  group_by(across(all_of(keys))) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n > 1L)

failure_audit <- test_runs |>
  group_by(N, model, hypothesis, procedure) |>
  summarise(
    n_stored = n(),
    n_failed = sum(is.na(reject)),
    n_rejected = sum(reject, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    n_missing = nr_of_iters - n_stored,
    n_valid = n_stored - n_failed,
    rejection_valid = n_rejected / na_if(n_valid, 0L),
    
    # Bounds over all intended runs:
    # unresolved runs treated as all non-rejections or all rejections
    rejection_min = n_rejected / nr_of_iters,
    rejection_max =
      (n_rejected + n_failed + n_missing) / nr_of_iters
  )

nrow(missing_tests)
duplicate_tests

failure_audit |>
  filter(N == 24) |>
  print(n = Inf)

test_runs |>
  filter(N == 24, is.na(reject)) |>
  count(model, hypothesis, procedure, error, sort = TRUE) |>
  print(n = Inf)








mod_use <- mod1
mod_use$covb <- V_PGEE

g <- multcomp::glht(mod_use, linfct = L_trt)

p_old <- as.numeric(
  summary(g, test = multcomp::adjusted("holm"))$test$pvalues
)

theta <- drop(L_trt %*% b_PGEE)
SE <- sqrt(diag(L_trt %*% V_PGEE %*% t(L_trt)))

p_new <- p.adjust(
  2 * pnorm(-abs(theta / SE)),
  method = "holm"
)

if (!isTRUE(all.equal(p_old, p_new, tolerance = 1e-8))) {
  print(data.frame(
    N = N,
    run = l,
    contrast = seq_along(p_old),
    p_old = p_old,
    p_new = p_new
  ))
  stop("glht and manual asymptotic Holm disagree.")
}


