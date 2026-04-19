#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_path <- sub("^--file=", "", script_arg)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

args <- commandArgs(trailingOnly = TRUE)

list_plot_rds <- function(dir_path, pattern) {
  if (!dir.exists(dir_path)) return(character())
  files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(character())
  info <- file.info(files)
  files[order(info$mtime, decreasing = TRUE)]
}

latest_four_model_rds <- function() {
  candidates <- list_plot_rds(
    file.path("results", "prior_only"),
    "^prior_only_plot_data_.*\\.rds$"
  )
  if (length(candidates) == 0) return(NA_character_)
  hits <- candidates[
    grepl("ct_fmax", basename(candidates)) &
      grepl("ct_softmax", basename(candidates)) &
      grepl("nct_fmax", basename(candidates)) &
      grepl("nct_softmax", basename(candidates))
  ]
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

rds_path <- if (length(args) >= 1) args[[1]] else latest_four_model_rds()
if (is.na(rds_path) || !file.exists(rds_path)) {
  stop("Could not find a prior-only RDS bundle for the four-model fmax/softmax sensitivity check.")
}

bundle <- readRDS(rds_path)
post_long <- bundle$post_long
sim_long <- bundle$sim_long

if (is.null(post_long) || is.null(sim_long)) {
  stop("The supplied RDS bundle does not contain post_long and sim_long.")
}

output_dir <- "figures_generated"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot_obj, filename, width, height, dpi = 300) {
  ggsave(
    filename = file.path(output_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    limitsize = FALSE,
    bg = "white"
  )
}

model_levels <- c("ct_fmax", "ct_softmax", "nct_fmax", "nct_softmax")
model_palette <- c(
  "ct_fmax" = "#D55E00",
  "ct_softmax" = "#A3A500",
  "nct_fmax" = "#0072B2",
  "nct_softmax" = "#CC79A7"
)

param_label_parse <- function(x) {
  labels <- c(
    "K" = "K",
    "r" = "r",
    "q" = "q",
    "sigma2" = "sigma^2",
    "tau2" = "tau^2",
    "Pt[21]" = "P[1987]",
    "Pt[22]" = "P[1988]",
    "Pt[23]" = "P[1989]"
  )
  unname(labels[x])
}

dens_long_fn <- function(df, out_col, x_limits, n_grid = 512) {
  parts <- split(df, list(df$model, df$param), drop = TRUE)
  out <- lapply(parts, function(dd) {
    if (nrow(dd) == 0) return(NULL)
    d <- density(dd$value, from = x_limits[1], to = x_limits[2], n = n_grid, na.rm = TRUE)
    tmp <- data.frame(
      model = as.character(dd$model[1]),
      param = as.character(dd$param[1]),
      idx = seq_along(d$x),
      x = d$x,
      dens = d$y
    )
    names(tmp)[names(tmp) == "dens"] <- out_col
    tmp
  })
  bind_rows(out)
}

dens_long_fn_param_limits <- function(df, out_col, param_limits, n_grid = 512) {
  parts <- split(df, list(df$model, df$param), drop = TRUE)
  out <- lapply(parts, function(dd) {
    if (nrow(dd) == 0) return(NULL)
    p <- as.character(dd$param[1])
    lim <- param_limits[[p]]
    if (is.null(lim) || any(!is.finite(lim))) return(NULL)
    d <- density(dd$value, from = lim[1], to = lim[2], n = n_grid, na.rm = TRUE)
    tmp <- data.frame(
      model = as.character(dd$model[1]),
      param = as.character(dd$param[1]),
      idx = seq_along(d$x),
      x = d$x,
      dens = d$y
    )
    names(tmp)[names(tmp) == "dens"] <- out_col
    tmp
  })
  bind_rows(out)
}

build_prior_diff_plot <- function(post_df, sim_df, params, filename,
                                  x_limits = c(0, 1), free_x_by_param = FALSE,
                                  width = 14, height = 8) {
  post_sub <- post_df %>% filter(param %in% params)
  sim_sub <- sim_df %>% filter(param %in% params)

  if (free_x_by_param && "q" %in% params) {
    post_sub <- post_sub %>%
      mutate(value = ifelse(as.character(param) == "q", log(pmax(value, .Machine$double.eps)), value))
    sim_sub <- sim_sub %>%
      mutate(value = ifelse(as.character(param) == "q", log(pmax(value, .Machine$double.eps)), value))
  }

  if (nrow(post_sub) == 0 || nrow(sim_sub) == 0) {
    stop("The supplied prior-only RDS does not contain the requested parameters.")
  }

  observed_models <- sort(unique(as.character(post_sub$model)))
  ct_models <- sort(unique(as.character(post_sub$model[grepl("^ct", as.character(post_sub$model))])))
  nct_models <- sort(unique(as.character(post_sub$model[grepl("^nct", as.character(post_sub$model))])))
  model_levels <- c(ct_models, nct_models)
  model_levels <- unique(c(model_levels, observed_models))
  ncol_prior <- max(1, max(length(ct_models), length(nct_models)))

  post_sub <- post_sub %>%
    mutate(
      model = factor(as.character(model), levels = model_levels),
      param = factor(as.character(param), levels = params)
    )
  sim_sub <- sim_sub %>%
    mutate(
      model = factor(as.character(model), levels = levels(post_sub$model)),
      param = factor(as.character(param), levels = levels(post_sub$param))
    )

  param_levels <- levels(post_sub$param)
  param_label_values <- param_label_parse(param_levels)
  param_label_values[param_levels == "q"] <- "log(q)"

  if (free_x_by_param) {
    param_limits_tbl <- bind_rows(
      post_sub %>% transmute(param = as.character(param), value = value),
      sim_sub %>% transmute(param = as.character(param), value = value)
    ) %>%
      group_by(param) %>%
      summarise(
        lo = min(value, na.rm = TRUE),
        hi = max(value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        span = hi - lo,
        pad = ifelse(span > 0, 0.03 * span, pmax(abs(lo), 1) * 0.03),
        lo = lo - pad,
        hi = hi + pad
      )
    param_limits <- setNames(
      lapply(seq_len(nrow(param_limits_tbl)), function(i) c(param_limits_tbl$lo[i], param_limits_tbl$hi[i])),
      param_limits_tbl$param
    )
    post_den <- dens_long_fn_param_limits(post_sub, "dens_post", param_limits = param_limits)
  } else {
    post_den <- dens_long_fn(post_sub, "dens_post", x_limits = x_limits)
  }

  sim_truth_base <- sim_sub %>% mutate(model = "__truth__")
  truth_den_base <- if (free_x_by_param) {
    dens_long_fn_param_limits(sim_truth_base, "dens_truth", param_limits = param_limits)
  } else {
    dens_long_fn(sim_truth_base, "dens_truth", x_limits = x_limits)
  }
  truth_den_base <- truth_den_base %>%
    transmute(param = as.character(param), idx = idx, x_truth = x, dens_truth = dens_truth)

  truth_models <- as.character(levels(post_sub$model))
  truth_den <- truth_den_base[rep(seq_len(nrow(truth_den_base)), each = length(truth_models)), , drop = FALSE]
  truth_den$model <- rep(truth_models, times = nrow(truth_den_base))
  truth_den <- as_tibble(truth_den) %>%
    select(model, param, idx, x_truth, dens_truth)

  den_df <- left_join(
    post_den,
    truth_den %>% select(model, param, idx, x_truth, dens_truth),
    by = c("model", "param", "idx")
  ) %>%
    mutate(
      dens_post = ifelse(is.na(dens_post), 0, dens_post),
      dens_truth = ifelse(is.na(dens_truth), 0, dens_truth),
      x_truth = ifelse(is.na(x_truth), x, x_truth),
      param = factor(param, levels = levels(post_sub$param)),
      model = factor(model, levels = levels(post_sub$model)),
      y0 = as.numeric(param)
    ) %>%
    group_by(model, param) %>%
    mutate(
      dens_post_s = dens_post / max(dens_post, na.rm = TRUE),
      dens_truth_s = dens_truth / max(dens_truth, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      y_post = y0 + 0.85 * dens_post_s,
      y_truth = y0 + 0.85 * dens_truth_s,
      ymin = pmin(y_post, y_truth),
      ymax = pmax(y_post, y_truth),
      diff_dir = ifelse(dens_post_s >= dens_truth_s, "post_gt_truth", "truth_gt_post"),
      diff_dir = factor(
        diff_dir,
        levels = c("post_gt_truth", "truth_gt_post"),
        labels = c("Realised prior > Truth", "Truth > Realised prior")
      )
    ) %>%
    group_by(model, param) %>%
    mutate(
      diff_run = cumsum(diff_dir != dplyr::lag(diff_dir, default = first(diff_dir)))
    ) %>%
    ungroup()

  diff_palette <- c(
    "Realised prior > Truth" = "#0B57A3",
    "Truth > Realised prior" = "#D04A02"
  )
  post_line_colour <- "#3A3A3A"
  truth_line_colour <- "#B80F0A"

  if (free_x_by_param) {
    den_df <- den_df %>%
      mutate(
        param_label = factor(
          ifelse(as.character(param) == "q", "log(q)", param_label_parse(as.character(param))),
          levels = param_label_values
        )
      )

    p <- ggplot() +
      geom_ribbon(
        data = den_df,
        aes(x = x, ymin = 0, ymax = dens_post_s, group = interaction(model, param)),
        fill = "grey70",
        alpha = 0.14,
        colour = NA
      ) +
      geom_ribbon(
        data = den_df,
        aes(
          x = x,
          ymin = pmin(dens_post_s, dens_truth_s),
          ymax = pmax(dens_post_s, dens_truth_s),
          fill = diff_dir,
          group = interaction(model, param, diff_run)
        ),
        alpha = 0.82,
        colour = NA
      ) +
      geom_line(
        data = den_df,
        aes(x = x, y = dens_post_s, group = interaction(model, param)),
        colour = post_line_colour,
        linewidth = 0.35,
        alpha = 0.75
      ) +
      geom_line(
        data = den_df,
        aes(x = x_truth, y = dens_truth_s, group = interaction(model, param)),
        colour = truth_line_colour,
        linewidth = 0.28,
        alpha = 0.75
      ) +
      facet_grid(
        model ~ param_label,
        scales = "free_x",
        labeller = labeller(param_label = label_parsed)
      ) +
      scale_fill_manual(values = diff_palette, name = "Difference") +
      labs(x = "Value", y = "Scaled density") +
      theme_bw() +
      theme(
        strip.text.x = element_text(angle = 0, size = 13, face = "bold"),
        strip.text.y = element_text(angle = 0, size = 13, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, margin = margin(t = 3)),
        axis.text.y = element_text(size = 12),
        panel.spacing.x = grid::unit(1.1, "lines"),
        legend.position = "bottom",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12)
      )
  } else {
    p <- ggplot() +
      geom_ribbon(
        data = den_df,
        aes(x = x, ymin = y0, ymax = y_post, group = interaction(model, param)),
        fill = "grey70",
        alpha = 0.14,
        colour = NA
      ) +
      geom_ribbon(
        data = den_df,
        aes(x = x, ymin = ymin, ymax = ymax, fill = diff_dir, group = interaction(model, param, diff_run)),
        alpha = 0.82,
        colour = NA
      ) +
      geom_line(
        data = den_df,
        aes(x = x, y = y_post, group = interaction(model, param)),
        colour = post_line_colour,
        linewidth = 0.35,
        alpha = 0.75
      ) +
      geom_line(
        data = den_df,
        aes(x = x_truth, y = y_truth, group = interaction(model, param)),
        colour = truth_line_colour,
        linewidth = 0.28,
        alpha = 0.75
      ) +
      facet_wrap(~ model, scales = "free_x", nrow = 2, ncol = ncol_prior) +
      scale_fill_manual(values = diff_palette, name = "Difference") +
      scale_y_continuous(
        breaks = seq_along(levels(post_sub$param)),
        labels = parse(text = param_label_values)
      ) +
      coord_cartesian(xlim = x_limits) +
      labs(x = "Value", y = NULL) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        strip.text = element_text(size = 13, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, margin = margin(t = 3)),
        axis.text.y = element_text(size = 12),
        panel.spacing.x = grid::unit(1.1, "lines"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12)
      )
  }

  save_plot(p, filename, width = width, height = height)
}

build_prior_diff_plot(
  post_df = post_long,
  sim_df = sim_long,
  params = paste0("Pt[", 21:23, "]"),
  filename = "prior_only_fmax_softmax_ridge_difference.png",
  x_limits = c(0, 1),
  free_x_by_param = FALSE,
  width = 14,
  height = 9.5
)

build_prior_diff_plot(
  post_df = post_long,
  sim_df = sim_long,
  params = c("K", "r", "q", "sigma2", "tau2"),
  filename = "prior_only_fmax_softmax_parameter_difference.png",
  x_limits = NULL,
  free_x_by_param = TRUE,
  width = 16,
  height = 10.5
)

build_prior_diff_plot(
  post_df = post_long,
  sim_df = sim_long,
  params = c("K", "r", "q", "sigma2", "tau2", "Pt[21]", "Pt[22]", "Pt[23]"),
  filename = "prior_only_fmax_softmax_combined_difference.png",
  x_limits = NULL,
  free_x_by_param = TRUE,
  width = 16,
  height = 12
)

post_overlay <- post_long %>%
  filter(param %in% paste0("Pt[", 21:23, "]")) %>%
  mutate(
    model = factor(as.character(model), levels = model_levels),
    param = factor(
      as.character(param),
      levels = c("Pt[23]", "Pt[22]", "Pt[21]"),
      labels = c("P[1989]", "P[1988]", "P[1987]")
    )
  )

overlay_den <- dens_long_fn(
  post_overlay,
  out_col = "density",
  x_limits = c(0, 1),
  n_grid = 512
) %>%
  mutate(
    model = factor(as.character(model), levels = model_levels),
    param = factor(
      as.character(param),
      levels = c("Pt[23]", "Pt[22]", "Pt[21]"),
      labels = c("P[1989]", "P[1988]", "P[1987]")
    )
  )

p_overlay <- ggplot(overlay_den, aes(x = x, y = density, colour = model, fill = model)) +
  geom_area(alpha = 0.12, position = "identity", linewidth = 0) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ param, ncol = 1, scales = "free_y", labeller = labeller(param = label_parsed)) +
  scale_colour_manual(values = model_palette, name = "Model") +
  scale_fill_manual(values = model_palette, name = "Model") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Value", y = "Density") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 13, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12)
  )

save_plot(
  p_overlay,
  "prior_only_fmax_softmax_overlay.png",
  width = 11,
  height = 9
)

message("Saved manuscript-style sensitivity figures from: ", rds_path)
