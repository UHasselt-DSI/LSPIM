###############################################################################
# summarise_simulation_results.R
#
# Post-processing script for the LSPIM paper.
#
# Main figures generated for the manuscript:
#   fig_type1_main_<scenario>.png
#       Type I error for both interaction hypotheses.
#       Black-and-white: symbol = hypothesis; linetype = method.
#       Includes: Asymptotic LPIM, LSPIM-Omnibus, LSPIM-Holm.
#
#   fig_type1_benchmark_<scenario>.png
#       Supplementary benchmark figure including nparLD. Faceted by hypothesis.
#
#   fig_variance_<scenario>.png
#       Variance calibration for treatment-effect parameters.
#
# The nparLD procedures are kept as external benchmarks because they test
# hypotheses formulated in terms of marginal relative treatment effects, not the
# same PIM parameters as the LSPIM-Holm.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
})

## ---------------------------------------------------------------------------
## User settings
## ---------------------------------------------------------------------------
results_dir <- "Simulation_output"          # folder containing the .Rdata files
out_dir <- file.path(results_dir, "simulation_summary_test_SIM_paper")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

alpha <- 0.05
nr_of_iters <- 10000
mc_band <- 1.96 * sqrt(alpha * (1 - alpha) / nr_of_iters)

sample_sizes <- c(16, 18, 20, 24, 26, 30, 40, 50, 60)

simulation_files <- tibble::tribble(
  ~scenario,      ~display_name,                         ~pattern,                                                                                ~n_visit,
  "PO_5",        "Proportional odds model, 5 visits",    "small_sample_type1_pim_gee_%s_correct_PGEE_POmodel_testwithmeanhypothesis.Rdata",       5,
  "Gaussian_4",  "Gaussian random intercept model, 4 visits",  "small_sample_type1_pim_gee_%s_correct_PGEE_score_unrounded_testwithmeanhypothesis_4timepoints.Rdata", 4
)

method_map <- c(
  "PIM_wald.test"   = "Asymptotic LPIM",
  "GEE_wald.test"   = "LSPIM-Omnibus",
  "GEE_glht"        = "LSPIM-Holm",
  "GEE_glht_global" = "Global F-test",
  "nparLD_W"        = "nparLD (WTS)",
  "nparLD_A"        = "nparLD (ATS)",
  "PIM"             = "Original PIM",
  "PIM_adjusted"    = "Asymptotic LPIM",
  "GEE"             = "LSPIM-Omnibus",
  "GEE_score"       = "LSPIM-Holm",
  "nparLD"          = "nparLD"
)

method_levels_main <- c("Asymptotic LPIM", "LSPIM-Omnibus", "LSPIM-Holm")
method_levels_benchmark <- c(method_levels_main, "nparLD (WTS)", "nparLD (ATS)", "Global F-test")
hypothesis_levels <- c("Constant Treatment Effect", "Constant Time Effect", "Marginal RTE interaction")

## ---------------------------------------------------------------------------
## Helper functions
## ---------------------------------------------------------------------------

load_rdata_env <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  e
}

get0e <- function(e, name) {
  if (exists(name, envir = e, inherits = FALSE)) get(name, envir = e) else NULL
}

mean_na <- function(x) mean(x, na.rm = TRUE)
var_na <- function(x) stats::var(x, na.rm = TRUE)

recode_method <- function(x) {
  out <- unname(method_map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

make_par_names <- function(p) {
  out <- c(paste0("T", seq_len(min(2, p))), paste0("A", seq_len(max(0, p - 2))))
  out[seq_len(p)]
}

# Interaction-test extraction follows the simulation code:
#   column 3  = PIM group/treatment-effect constancy, between-within df
#   column 7  = PIM time-effect constancy, between-within df
#   column 11 = GEE group/treatment-effect constancy, between-within df
#   column 15 = GEE time-effect constancy, between-within df
#   column 17 = nparLD WTS interaction
#   column 18 = nparLD ATS interaction
#
# score_inter_group/time are Holm deviation-contrast rejection indicators.
# score_inter_group_global/time_global are optional global F-test p-values.
summarise_interaction <- function(e, scenario, N) {
  store <- get0e(e, "store_interaction_tests")
  if (is.null(store)) return(tibble::tibble())
  
  type1 <- apply(store, 2, function(x) mean(x < alpha, na.rm = TRUE))
  n_col <- length(type1)
  get_type1 <- function(j) if (n_col >= j) type1[j] else NA_real_
  
  out <- tibble::tibble(
    scenario = scenario,
    sampsize = N,
    method = c("PIM_wald.test", "PIM_wald.test",
               "GEE_wald.test", "GEE_wald.test",
               "nparLD_W", "nparLD_A"),
    hypothesis = c("Constant Treatment Effect", "Constant Time Effect",
                   "Constant Treatment Effect", "Constant Time Effect",
                   "Marginal RTE interaction", "Marginal RTE interaction"),
    endpoint = c("between-within Wald/F", "between-within Wald/F",
                 "between-within Wald/F", "between-within Wald/F",
                 "rank-based benchmark", "rank-based benchmark"),
    TypeI = c(get_type1(3), get_type1(7), get_type1(11), get_type1(15), get_type1(17), get_type1(18))
  )
  
  score_inter_group <- get0e(e, "score_inter_group")
  score_inter_time <- get0e(e, "score_inter_time")
  
  if (!is.null(score_inter_group)) {
    out <- dplyr::bind_rows(out, tibble::tibble(
      scenario = scenario, sampsize = N, method = "GEE_glht",
      hypothesis = "Constant Treatment Effect",
      endpoint = "Holm deviation contrasts",
      TypeI = mean(score_inter_group, na.rm = TRUE)
    ))
  }
  if (!is.null(score_inter_time)) {
    out <- dplyr::bind_rows(out, tibble::tibble(
      scenario = scenario, sampsize = N, method = "GEE_glht",
      hypothesis = "Constant Time Effect",
      endpoint = "Holm deviation contrasts",
      TypeI = mean(score_inter_time, na.rm = TRUE)
    ))
  }
  
  score_inter_group_global <- get0e(e, "score_inter_group_global")
  score_inter_time_global <- get0e(e, "score_inter_time_global")
  
  if (!is.null(score_inter_group_global)) {
    out <- dplyr::bind_rows(out, tibble::tibble(
      scenario = scenario, sampsize = N, method = "GEE_glht_global",
      hypothesis = "Constant Treatment Effect",
      endpoint = "global F-test from glht",
      TypeI = mean(score_inter_group_global <= alpha, na.rm = TRUE)
    ))
  }
  if (!is.null(score_inter_time_global)) {
    out <- dplyr::bind_rows(out, tibble::tibble(
      scenario = scenario, sampsize = N, method = "GEE_glht_global",
      hypothesis = "Constant Time Effect",
      endpoint = "global F-test from glht",
      TypeI = mean(score_inter_time_global <= alpha, na.rm = TRUE)
    ))
  }
  
  out |>
    dplyr::mutate(
      contrast = hypothesis, # kept for backward compatibility with older plotting code
      method_label = recode_method(method),
      method_label = factor(method_label, levels = method_levels_benchmark),
      hypothesis = factor(hypothesis, levels = hypothesis_levels),
      contrast = factor(contrast, levels = hypothesis_levels)
    ) |>
    dplyr::filter(!is.na(TypeI))
}

summarise_variance_bias <- function(e, scenario, N) {
  store_parms_mod2 <- get0e(e, "store_parms_mod2")
  store_parms_mod2_adj <- store_parms_mod2
  store_parms_gee <- get0e(e, "store_parms_gee")
  store_se_mod2 <- get0e(e, "store_se_mod2")
  store_se_mod2_adj <- get0e(e, "store_se_mod2_adj")
  store_se_gee <- get0e(e, "store_se_gee")
  
  empty <- list(bias = tibble::tibble(), variance = tibble::tibble())
  if (is.null(store_parms_mod2) || is.null(store_parms_gee)) return(empty)
  
  p_pim <- ncol(store_parms_mod2)
  p_gee <- ncol(store_parms_gee)
  par_pim <- make_par_names(p_pim)
  par_gee <- make_par_names(p_gee)
  
  bias <- dplyr::bind_rows(
    tibble::tibble(scenario = scenario, sampsize = N, method = "PIM", par = par_pim,
                   estimate_mean = apply(store_parms_mod2, 2, mean_na)),
    tibble::tibble(scenario = scenario, sampsize = N, method = "GEE", par = par_gee,
                   estimate_mean = apply(store_parms_gee, 2, mean_na))
  ) |>
    dplyr::mutate(method_label = recode_method(method))
  
  variance <- tibble::tibble()
  if (!is.null(store_se_mod2)) {
    variance <- dplyr::bind_rows(variance, tibble::tibble(
      scenario = scenario, sampsize = N, method = "PIM", par = par_pim,
      model_variance = apply(store_se_mod2, 2, mean_na),
      empirical_variance = apply(store_parms_mod2, 2, var_na)
    ))
  }
  if (!is.null(store_se_mod2_adj) && !is.null(store_parms_mod2_adj)) {
    variance <- dplyr::bind_rows(variance, tibble::tibble(
      scenario = scenario, sampsize = N, method = "PIM_adjusted", par = par_pim,
      model_variance = apply(store_se_mod2_adj, 2, mean_na),
      empirical_variance = apply(store_parms_mod2_adj, 2, var_na)
    ))
  }
  if (!is.null(store_se_gee)) {
    variance <- dplyr::bind_rows(variance, tibble::tibble(
      scenario = scenario, sampsize = N, method = "GEE", par = par_gee,
      model_variance = apply(store_se_gee, 2, mean_na),
      empirical_variance = apply(store_parms_gee, 2, var_na)
    ))
  }
  variance <- variance |>
    dplyr::mutate(
      variance_ratio = model_variance / empirical_variance,
      method_label = recode_method(method),
      method_label = factor(method_label, levels = c("Asymptotic LPIM", "Adjusted PIM covariance", "LSPIM-Omnibus"))
    )
  
  list(bias = bias, variance = variance)
}

summarise_file <- function(path, scenario, N) {
  e <- load_rdata_env(path)
  list(
    interaction = summarise_interaction(e, scenario, N),
    vb = summarise_variance_bias(e, scenario, N)
  )
}

## ---------------------------------------------------------------------------
## Read files and create summaries
## ---------------------------------------------------------------------------

all_interaction <- list()
all_bias <- list()
all_variance <- list()
missing_files <- character()

for (i in seq_len(nrow(simulation_files))) {
  scenario_i <- simulation_files$scenario[i]
  pattern_i <- simulation_files$pattern[i]
  message("Processing scenario: ", scenario_i)
  
  for (N in sample_sizes) {
    path <- file.path(results_dir, sprintf(pattern_i, N))
    if (!file.exists(path)) {
      missing_files <- c(missing_files, path)
      next
    }
    res <- summarise_file(path, scenario_i, N)
    all_interaction[[length(all_interaction) + 1]] <- res$interaction
    all_bias[[length(all_bias) + 1]] <- res$vb$bias
    all_variance[[length(all_variance) + 1]] <- res$vb$variance
  }
}

interaction_results <- dplyr::bind_rows(all_interaction)
bias_results <- dplyr::bind_rows(all_bias)
variance_results <- dplyr::bind_rows(all_variance)

if (length(missing_files) > 0) {
  writeLines(missing_files, con = file.path(out_dir, "missing_files.txt"))
  warning(length(missing_files), " expected .Rdata files were not found. See missing_files.txt.")
}

## ---------------------------------------------------------------------------
## Figures
## ---------------------------------------------------------------------------

scenario_title <- function(scenario_name) {
  x <- simulation_files$display_name[match(scenario_name, simulation_files$scenario)]
  ifelse(is.na(x), scenario_name, x)
}

base_type1_layers <- function() {
  list(
    ggplot2::geom_hline(yintercept = alpha),
    ggplot2::geom_hline(yintercept = c(alpha - mc_band, alpha + mc_band), linetype = "dashed",color="lightgrey"),
    ggplot2::scale_x_continuous(breaks = sample_sizes),
    ggplot2::labs(
      x = "Total sample size (balanced between treatment groups)",
      y = "Empirical Type I error",
      linetype = "Method",
      shape = "Hypothesis"
    ),
    ggplot2::theme_bw(base_size = 12),
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold")
    )
  )
}

plot_type1_main <- function(dat, scenario_name) {
  dat_sub <- dat |>
    dplyr::filter(scenario == scenario_name) |>
    dplyr::filter(method_label %in% method_levels_main) |>
    dplyr::filter(hypothesis %in% c("Constant Treatment Effect", "Constant Time Effect")) |>
    dplyr::mutate(
      method_label = factor(as.character(method_label), levels = method_levels_main),
      hypothesis = factor(as.character(hypothesis), levels = c("Constant Treatment Effect", "Constant Time Effect"))
    )
  
  if (nrow(dat_sub) == 0) return(NULL)
  
  p <- ggplot2::ggplot(dat_sub, ggplot2::aes(x = sampsize, y = TypeI,
                                             linetype = method_label,
                                             shape = hypothesis,
                                             group = interaction(method_label, hypothesis))) +
    ggplot2::geom_line(colour = "black", linewidth = 0.6) +
    ggplot2::geom_point(size = 2.4, colour = "black") +
    ggplot2::scale_shape_manual(values = c(
      "Constant Treatment Effect" = 16,
      "Constant Time Effect" = 17
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "Asymptotic LPIM" = "solid",
      "LSPIM-Omnibus" = "dashed",
      "LSPIM-Holm" = "dotdash"
    )) +
    base_type1_layers() +
    ggplot2::ggtitle(paste("Type I error:", scenario_title(scenario_name)))
  
  ggplot2::ggsave(file.path(out_dir, paste0("fig_type1_main_", scenario_name, ".png")), p,
                  width = 8.5, height = 5.5, dpi = 300)
  p
}

plot_type1_benchmark <- function(dat, scenario_name) {
  dat_sub <- dat |>
    dplyr::filter(scenario == scenario_name) |>
    dplyr::filter(method_label %in% c(method_levels_main, "nparLD (WTS)", "nparLD (ATS)")) |>
    dplyr::mutate(
      method_label = factor(as.character(method_label), levels = c(method_levels_main, "nparLD (WTS)", "nparLD (ATS)")),
      hypothesis = factor(as.character(hypothesis), levels = hypothesis_levels)
    )
  
  if (nrow(dat_sub) == 0) return(NULL)
  
  p <- ggplot2::ggplot(dat_sub, ggplot2::aes(x = sampsize, y = TypeI,
                                             linetype = method_label,
                                             shape = hypothesis,
                                             group = interaction(method_label, hypothesis))) +
    ggplot2::geom_line(colour = "black", linewidth = 0.6) +
    ggplot2::geom_point(size = 2.2, colour = "black") +
    ggplot2::scale_shape_manual(values = c(
      "Constant Treatment Effect" = 16,
      "Constant Time Effect" = 17,
      "Marginal RTE interaction" = 15
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "Asymptotic LPIM" = "solid",
      "LSPIM-Omnibus" = "dashed",
      "LSPIM-Holm" = "dotdash",
      "nparLD (WTS)" = "longdash",
      "nparLD (ATS)" = "twodash"
    )) +
    base_type1_layers() +
    ggplot2::ggtitle(paste("Benchmark Type I error:", scenario_title(scenario_name)))
  
  ggplot2::ggsave(file.path(out_dir, paste0("fig_type1_benchmark_", scenario_name, ".png")), p,
                  width = 10, height = 5.5, dpi = 300)
  p
}

plot_variance <- function(dat, scenario_name) {
  dat_sub <- dat |>
    dplyr::filter(scenario == scenario_name) |>
    dplyr::filter(par %in%c("T1","A1")) |>
    dplyr::filter(method_label %in% c("Asymptotic LPIM", "LSPIM-Omnibus"))%>%
    mutate(
      method_label = factor(
        method_label,
        levels = c("LSPIM-Omnibus","Asymptotic LPIM"),
        labels = c(
          "LSPIM",
          "Asymptotic LPIM"
        )
      ),
      par = factor(
        par,
        levels = c("T1","A1"),
        labels = c(
          "Time",
          "Treatment"
        )
      )
    )
  
  if (nrow(dat_sub) == 0) return(NULL)
  
  p <- ggplot2::ggplot(dat_sub, ggplot2::aes(x = sampsize, y = variance_ratio,
                                             colour = par,
                                             linetype = method_label,
                                             shape = par,
                                             group = interaction(method, par)))+
    geom_line(colour = "black", linewidth = 0.8) +
    geom_point(colour = "black", size = 2.5) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    scale_shape_manual(
      values = c(
        "Treatment" = 16,  # filled circle
        "Time" = 17   # filled triangle
      )
    ) +
    labs(      x = "Sample size",
               y = "Model variance / empirical variance",
               linetype = "Method",
               shape = "Parameter"
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    )+
    ggplot2::ggtitle(paste(scenario_title(scenario_name)))
  
  ggplot2::ggsave(file.path(out_dir, paste0("fig_variance_", scenario_name, ".png")), p,
                  width = 8.5, height = 5.5, dpi = 300)
  p
}

if (nrow(interaction_results) > 0) {
  for (sc in unique(interaction_results$scenario)) {
    plot_type1_main(interaction_results, sc)
    plot_type1_benchmark(interaction_results, sc)
  }
}
if (nrow(variance_results) > 0) {
  for (sc in unique(variance_results$scenario)) plot_variance(variance_results, sc)
}

message("Done. Results written to: ", normalizePath(out_dir))






