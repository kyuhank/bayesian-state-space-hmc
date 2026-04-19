
clean_filename <- function(x) {
  x <- gsub("/", "_", x)
  x <- sub("\\.R$", "", x)
  return(x)
}

format_parameter_label_plotmath <- function(x) {
  x <- as.character(x)
  label_map <- c(
    "K" = "K",
    "r" = "r",
    "q" = "q",
    "sigma2" = "sigma^2",
    "tau2" = "tau^2",
    "MSY" = "MSY",
    "Bmsy" = "B[MSY]",
    "Hmsy" = "H[MSY]",
    "loglik" = "log*L",
    "Pt[20]" = "P[1986]",
    "Pt[21]" = "P[1987]",
    "Pt[22]" = "P[1988]",
    "Pt[23]" = "P[1989]",
    "Bt[20]" = "B[1986]",
    "Bt[21]" = "B[1987]",
    "Bt[22]" = "B[1988]",
    "Bt[23]" = "B[1989]",
    "BtoverBmsy[23]" = "B[1989] / B[MSY]",
    "HtoverHmsy[23]" = "H[1989] / H[MSY]"
  )
  mapped <- ifelse(x %in% names(label_map), unname(label_map[x]), x)

  if (exists("years", inherits = TRUE)) {
    years_vec <- get("years", inherits = TRUE)
    pt_idx <- stringr::str_match(x, "^Pt\\[(\\d+)\\]$")
    bt_idx <- stringr::str_match(x, "^Bt\\[(\\d+)\\]$")
    bb_idx <- stringr::str_match(x, "^BtoverBmsy\\[(\\d+)\\]$")
    hh_idx <- stringr::str_match(x, "^HtoverHmsy\\[(\\d+)\\]$")

    if (!is.na(pt_idx[1, 2])) {
      idx <- as.integer(pt_idx[1, 2])
      if (idx >= 1 && idx <= length(years_vec)) mapped <- paste0("P[", years_vec[idx], "]")
    } else if (!is.na(bt_idx[1, 2])) {
      idx <- as.integer(bt_idx[1, 2])
      if (idx >= 1 && idx <= length(years_vec)) mapped <- paste0("B[", years_vec[idx], "]")
    } else if (!is.na(bb_idx[1, 2])) {
      idx <- as.integer(bb_idx[1, 2])
      if (idx >= 1 && idx <= length(years_vec)) mapped <- paste0("B[", years_vec[idx], "] / B[MSY]")
    } else if (!is.na(hh_idx[1, 2])) {
      idx <- as.integer(hh_idx[1, 2])
      if (idx >= 1 && idx <= length(years_vec)) mapped <- paste0("H[", years_vec[idx], "] / H[MSY]")
    }
  }

  mapped
}

parameter_display_order <- function() {
  c(
    "K", "r", "q", "sigma2", "tau2", "MSY",
    "Pt[23]", "BtoverBmsy[23]", "HtoverHmsy[23]",
    "Pt[21]", "Pt[22]",
    "Pt[20]", "Bt[20]", "Bt[21]", "Bt[22]", "Bt[23]",
    "Bmsy", "Hmsy",
    "loglik"
  )
}

format_sbc_run_label <- function(x) {
  # Expected stem formats:
  # sbc_<model>.stan_<adaptDelta>_<tau2PriorName>_<nSims>
  # <legacySBCtype>_<model>.stan_<adaptDelta>_<tau2PriorName>_<nSims>
  m <- stringr::str_match(x, "^(?:(?:sbc)|(?:\\w*SBC(?:dep)?))_(.+?)\\.stan_[^_]+_([^_]+)_.+$")
  if (is.na(m[1, 1])) {
    m <- stringr::str_match(x, "^(.+?)\\.stan_[^_]+_([^_]+)_.+$")
  }
  if (is.na(m[1, 1])) {
    return(x)
  }
  model_name <- m[1, 2]
  cv_value <- m[1, 3]
  if (identical(cv_value, "no")) {
    paste0(model_name, "_no")
  } else {
    paste0(model_name, "_cv", cv_value)
  }
}




prepare_ecdf_diff_overlay_data <- function(sbc_results_list, variables = NULL,
                                           prob = 0.95, colors = NULL,
                                           drop_loglik = TRUE, ...) {
  # Set up color palette with better distinction than viridis
  if (is.null(colors)) {
    n_results <- length(sbc_results_list)
    
    # Use RColorBrewer qualitative palettes for better color distinction
    if (n_results <= 3) {
      colors <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1:n_results]
    } else if (n_results <= 8) {
      colors <- RColorBrewer::brewer.pal(n = n_results, name = "Set1")
    } else if (n_results <= 12) {
      colors <- RColorBrewer::brewer.pal(n = 12, name = "Paired")[1:n_results]
    } else {
      # For many results, combine multiple palettes
      colors <- c(
        RColorBrewer::brewer.pal(8, "Set1"),
        RColorBrewer::brewer.pal(8, "Dark2"),
        RColorBrewer::brewer.pal(8, "Set2")
      )[1:n_results]
    }
  } else if (length(colors) < length(sbc_results_list)) {
    # Extend with RColorBrewer Set1 if not enough colors provided
    n_needed <- length(sbc_results_list)
    if (n_needed <= 8) {
      colors <- RColorBrewer::brewer.pal(n = max(3, n_needed), name = "Set1")
    } else {
      colors <- scales::hue_pal()(length(sbc_results_list))
    }
  }
  
  # Extract data from each SBC_results object
  all_ecdf_data <- list()
  limits_df_trans <- NULL
  
  for (i in seq_along(sbc_results_list)) {
    # Get result name from list names
    result_name <- names(sbc_results_list)[i]
    if (is.null(result_name) || result_name == "") {
      result_name <- paste0("Result_", i)
    }
    
    # Extract ECDF data using SBC internal function
    ecdf_data <- SBC:::data_for_ecdf_plots(
      sbc_results_list[[i]], 
      variables = variables, 
      prob = prob, 
      ...
    )
    
    # Calculate ECDF differences from theoretical uniform distribution
    ecdf_df <- dplyr::mutate(ecdf_data$ecdf_df, 
                             z_diff = ecdf - z, 
                             result_name = result_name,
                             result_id = i)
    
    # Use theoretical limits from first result only (should be identical across results)
  if (i == 1) {
      limits_df_trans <- dplyr::mutate(ecdf_data$limits_df, 
                                       ymax = upper - uniform_val, 
                                       ymin = lower - uniform_val)
    }
    
    all_ecdf_data[[i]] <- ecdf_df
  }
  
  # Combine all ECDF data
  combined_ecdf_df <- dplyr::bind_rows(all_ecdf_data)
  if (drop_loglik && "group" %in% names(combined_ecdf_df)) {
    combined_ecdf_df <- combined_ecdf_df %>%
      dplyr::filter(as.character(group) != "loglik")
  }
  if (!is.null(limits_df_trans) && !("group" %in% names(limits_df_trans)) && ("variable" %in% names(limits_df_trans))) {
    limits_df_trans <- limits_df_trans %>%
      dplyr::mutate(group = variable)
  }
  if (drop_loglik && !is.null(limits_df_trans) && "group" %in% names(limits_df_trans)) {
    limits_df_trans <- limits_df_trans %>%
      dplyr::filter(as.character(group) != "loglik")
  }
  raw_result_names <- names(sbc_results_list)
  pretty_result_names <- vapply(raw_result_names, format_sbc_run_label, character(1))
  pretty_result_names <- make.unique(pretty_result_names, sep = "_")
  legend_labels <- sub("^.*_(cv[^_]+|no)$", "\\1", pretty_result_names)
  legend_labels <- ifelse(grepl("^cv", legend_labels), sub("^cv", "", legend_labels), legend_labels)
  legend_order <- c("0.1", "0.3", "0.5", "no")
  legend_levels <- c(intersect(legend_order, unique(legend_labels)),
                     setdiff(unique(legend_labels), legend_order))
  label_map <- stats::setNames(legend_labels, raw_result_names)
  
  # Convert result_name to factor for legend ordering
  combined_ecdf_df$result_name <- factor(
    label_map[combined_ecdf_df$result_name],
    levels = legend_levels
  )
  if ("group" %in% names(combined_ecdf_df)) {
    raw_groups <- as.character(combined_ecdf_df$group)
    group_levels <- c(intersect(parameter_display_order(), unique(raw_groups)),
                      setdiff(unique(raw_groups), parameter_display_order()))
    label_levels <- vapply(group_levels, format_parameter_label_plotmath, character(1))
    combined_ecdf_df$group <- factor(
      vapply(raw_groups, format_parameter_label_plotmath, character(1)),
      levels = unique(label_levels)
    )
  }
  if (!is.null(limits_df_trans) && "group" %in% names(limits_df_trans)) {
    raw_groups <- as.character(limits_df_trans$group)
    group_levels <- c(intersect(parameter_display_order(), unique(raw_groups)),
                      setdiff(unique(raw_groups), parameter_display_order()))
    label_levels <- vapply(group_levels, format_parameter_label_plotmath, character(1))
    limits_df_trans$group <- factor(
      vapply(raw_groups, format_parameter_label_plotmath, character(1)),
      levels = unique(label_levels)
    )
  } else if (!is.null(limits_df_trans) &&
             "group" %in% names(combined_ecdf_df)) {
    # SBC's theoretical ECDF limits can be returned once for all variables.
    # Replicate the same envelope across facet groups so each panel gets the
    # same reference oval/band.
    group_levels <- levels(combined_ecdf_df$group)
    if (is.null(group_levels)) {
      group_levels <- unique(as.character(combined_ecdf_df$group))
    }
    limits_df_trans <- dplyr::bind_rows(lapply(group_levels, function(g) {
      limits_df_trans %>%
        dplyr::mutate(group = factor(g, levels = group_levels))
    }))
  }
  
  # Create color mapping
  if (length(colors) < length(legend_levels)) {
    colors <- scales::hue_pal()(length(legend_levels))
  }
  color_mapping <- setNames(colors[seq_along(legend_levels)], legend_levels)

  list(
    combined_ecdf_df = combined_ecdf_df,
    limits_df_trans = limits_df_trans,
    color_mapping = color_mapping,
    legend_levels = legend_levels
  )
}

plot_ecdf_diff_overlay <- function(sbc_results_list, variables = NULL, 
                                   line_alpha = 0.85, line_size = 0.85, 
                                   alpha = 0.33, prob = 0.95, 
                                   colors = NULL, 
                                   show_oval_band = FALSE,
                                   oval_fill = "lightblue",
                                   oval_alpha = 0.33,
                                   show_cdf_border = TRUE,
                                   show_points = TRUE,
                                   point_alpha = 0.95,
                                   point_size = 1.9,
                                   point_stroke = 0.55,
                                   point_n = 22,
                                   cdf_border_color = "blue",
                                   cdf_border_alpha = 0.7,
                                   cdf_border_size = 0.6,
                                   facet_ncol = 5,
                                   base_size = 13,
                                   legend_ncol = NULL,
                                   ...) {
  ecdf_prepped <- prepare_ecdf_diff_overlay_data(
    sbc_results_list = sbc_results_list,
    variables = variables,
    prob = prob,
    colors = colors,
    ...
  )
  combined_ecdf_df <- ecdf_prepped$combined_ecdf_df
  limits_df_trans <- ecdf_prepped$limits_df_trans
  color_mapping <- ecdf_prepped$color_mapping
  legend_levels <- ecdf_prepped$legend_levels
  ecdf_points_df <- combined_ecdf_df %>%
    dplyr::group_by(group, result_name) %>%
    dplyr::slice(unique(round(seq(1, dplyr::n(), length.out = min(point_n, dplyr::n()))))) %>%
    dplyr::ungroup()

  p <- ggplot()

  if (!is.null(limits_df_trans)) {
    p <- p +
      geom_ribbon(
        data = limits_df_trans,
        aes(x = x, ymax = ymax, ymin = ymin, group = group),
        fill = "lightblue",
        alpha = alpha
      )

    if (show_cdf_border) {
      p <- p +
        geom_line(
          data = limits_df_trans,
          aes(x = x, y = ymax, group = group),
          color = cdf_border_color,
          alpha = cdf_border_alpha,
          size = cdf_border_size,
          linetype = "solid"
        ) +
        geom_line(
          data = limits_df_trans,
          aes(x = x, y = ymin, group = group),
          color = cdf_border_color,
          alpha = cdf_border_alpha,
          size = cdf_border_size,
          linetype = "solid"
        ) +
        geom_hline(
          yintercept = 0,
          color = cdf_border_color,
          alpha = cdf_border_alpha * 1.2,
          size = cdf_border_size * 1.2,
          linetype = "solid"
        )
    }
  }

  p <- p +
    geom_step(
      data = combined_ecdf_df,
      aes(
        x = z, y = z_diff,
        color = result_name,
        group = interaction(variable, result_name)
      ),
      alpha = line_alpha, size = line_size
    )

  if (show_points) {
    p <- p +
      geom_point(
        data = ecdf_points_df,
        aes(
          x = z, y = z_diff,
          color = result_name,
          group = interaction(group, result_name)
        ),
        alpha = point_alpha,
        size = point_size,
        shape = 1,
        stroke = point_stroke
      )
  }

  p <- p +
    scale_color_manual(name = "CV", values = color_mapping, breaks = legend_levels) +
    scale_x_continuous(
      breaks = c(0, 0.5, 1.0),
      labels = c("0.0", "0.5", "1.0"),
      limits = c(0, 1)
    ) +
    facet_wrap(~group, scales = "free_y", ncol = facet_ncol, labeller = ggplot2::label_parsed) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    theme_minimal(base_size = base_size) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = base_size + 1),
      legend.text = element_text(size = base_size),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size + 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.title = element_text(size = base_size + 2),
      axis.text = element_text(size = base_size),
      axis.text.x = element_text(margin = margin(t = 2)),
      panel.spacing = grid::unit(0.9, "lines"),
      legend.key.width = grid::unit(1.4, "lines"),
      legend.key.height = grid::unit(0.9, "lines")
    ) +
    guides(color = guide_legend(
      override.aes = list(size = 1.8, alpha = 1.0),
      ncol = if (is.null(legend_ncol)) min(4, length(sbc_results_list)) else legend_ncol
    )) +
    xlab("PIT") +
    ylab("ECDF Difference")

  return(p)
}
