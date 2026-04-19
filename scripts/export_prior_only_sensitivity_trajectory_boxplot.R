#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
setwd(repo_root)

source("data/Albacore_Data.R")

models <- c("ct_fmax", "ct_softmax", "nct_fmax", "nct_softmax")
pt_vars <- paste0("Pt[", seq_len(ntimes + 1), "]")
param_vars <- c("K", "r", "q", "sigma2", "tau2")

model_palette <- c(
  "ct_fmax" = "#E66100",
  "ct_softmax" = "#009E73",
  "nct_fmax" = "#0072B2",
  "nct_softmax" = "#CC79A7"
)

latest_stem <- function(model_name) {
  csvs <- list.files(file.path("fitcsv", model_name), pattern = "\\.csv$", full.names = TRUE)
  if (length(csvs) == 0) {
    stop(
      "No CSV files found for ", model_name,
      ". Run `make sensitivity-softmax` first to generate the required fitcsv output."
    )
  }
  stems <- sub("-[0-9]+-[0-9]+\\.csv$", "", basename(csvs))
  info <- file.info(csvs)
  stems[order(info$mtime, decreasing = TRUE)][1]
}

csv_group <- function(model_name) {
  stem <- latest_stem(model_name)
  files <- list.files(
    file.path("fitcsv", model_name),
    pattern = paste0("^", stem, "-[0-9]+-[0-9]+\\.csv$"),
    full.names = TRUE
  )
  if (length(files) == 0) stop("Could not resolve latest CSV group for ", model_name)
  sort(files)
}

box_stats <- function(x) {
  qs <- quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE, na.rm = TRUE)
  iqr <- qs[3] - qs[1]
  lo_lim <- qs[1] - 1.5 * iqr
  hi_lim <- qs[3] + 1.5 * iqr
  inside <- x[x >= lo_lim & x <= hi_lim]
  c(
    ymin = min(inside, na.rm = TRUE),
    lower = qs[1],
    middle = qs[2],
    upper = qs[3],
    ymax = max(inside, na.rm = TRUE)
  )
}

message("Reading latest prior-only fit CSVs for the four-model sensitivity check...")

stats_df <- bind_rows(lapply(models, function(model_name) {
  files <- csv_group(model_name)
  fit <- read_cmdstan_csv(files)
  draws_df <- posterior::as_draws_df(fit$post_warmup_draws)
  draws_df$.draw_id <- seq_len(nrow(draws_df))

  long_df <- draws_df %>%
    select(all_of(pt_vars), .draw_id) %>%
    pivot_longer(
      cols = all_of(pt_vars),
      names_to = "param",
      values_to = "value"
    ) %>%
    mutate(
      model = model_name,
      t_idx = as.integer(gsub("^Pt\\[|\\]$", "", param))
    )

  stat_mat <- long_df %>%
    group_by(model, t_idx) %>%
    summarise(stats = list(box_stats(value)), .groups = "drop") %>%
    mutate(stats = lapply(stats, as.list)) %>%
    tidyr::unnest_wider(stats)

  stat_mat
}))

param_df <- bind_rows(lapply(models, function(model_name) {
  files <- csv_group(model_name)
  fit <- read_cmdstan_csv(files)
  draws_df <- posterior::as_draws_df(fit$post_warmup_draws)

  draws_df %>%
    select(all_of(param_vars)) %>%
    mutate(model = model_name) %>%
    pivot_longer(
      cols = all_of(param_vars),
      names_to = "param",
      values_to = "value"
    ) %>%
    mutate(
      value = ifelse(param == "q", log(value), value),
      param = dplyr::case_when(
        param == "K" ~ "K",
        param == "r" ~ "r",
        param == "q" ~ "log(q)",
        param == "sigma2" ~ "sigma^2",
        param == "tau2" ~ "tau^2",
        TRUE ~ param
      )
    )
}))

param_df <- param_df %>%
  mutate(
    model = factor(model, levels = models),
    param = factor(param, levels = c("K", "r", "log(q)", "sigma^2", "tau^2"))
  )

state_labels <- c(as.character(min(years) - 1), as.character(years))
stats_df <- stats_df %>%
  mutate(
    model = factor(model, levels = models),
    state = factor(state_labels[t_idx], levels = state_labels)
  )

p <- ggplot(
  stats_df,
  aes(
    x = state,
    ymin = ymin,
    lower = lower,
    middle = middle,
    upper = upper,
    ymax = ymax,
    fill = model,
    colour = model
  )
) +
  geom_boxplot(
    stat = "identity",
    position = position_dodge2(width = 0.82, preserve = "single"),
    width = 0.72,
    linewidth = 0.38,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = model_palette, name = "Model", guide = "none") +
  scale_colour_manual(values = model_palette, name = "Model", guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Year",
    y = "Relative biomass"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8.5),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12)
  )

p_param <- ggplot(
  param_df,
  aes(x = value, colour = model)
) +
  geom_density(linewidth = 1.15) +
  facet_wrap(~ param, scales = "free", nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = model_palette, name = "Model") +
  labs(
    x = NULL,
    y = "Density"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12)
  )

p_combined <- p_param / p +
  plot_layout(heights = c(1.1, 2.4), guides = "collect") &
  theme(legend.position = "bottom")

out_dir <- "figures_generated"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "prior_only_fmax_softmax_trajectory_boxplot.png")
out_file_combined <- file.path(out_dir, "prior_only_fmax_softmax_combined_overlay_boxplot.png")

ggsave(
  filename = out_file,
  plot = p,
  width = 15,
  height = 6.8,
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)

message("Saved: ", out_file)

ggsave(
  filename = out_file_combined,
  plot = p_combined,
  width = 17,
  height = 10.4,
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)

message("Saved: ", out_file_combined)
