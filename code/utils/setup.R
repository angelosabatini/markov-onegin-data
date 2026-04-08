## ============================================================
## File: setup.R
##
## Purpose:
##   Centralized setup for the analysis environment.
##   Loads required packages, defines global plotting themes,
##   color palettes, and shared helper options.
##
## Scope:
##   - Sourced at the beginning of driver scripts
##   - Provides shared visual and stylistic conventions
##
## Notes:
##   - No analysis logic is defined here
##   - Modifying this file affects all figures
## ============================================================

## Load libraries
## Core tidyverse (includes dplyr, ggplot2, readr, etc.)
library(tidyverse)

## String and regex handling
library(stringi)        # Full Unicode/string control
library(rebus)          # Human-readable regex

## Text processing and NLP
library(tidytext)
library(tokenizers)
library(udpipe)

## Statistical tools
library(ppcor)
library(Kendall)
library(gtools)

## Visualization helpers
library(ggpubr)
library(patchwork)

## Reporting and formatting
library(knitr)
library(kableExtra)
library(broom)

## ---- Setup: Theme, Palette, Plotting Utilities ----

## Define a custom ggplot2 theme

my_theme <- function(output = c("pdf", "png"), text_scale = 1) {
  output <- match.arg(output)

  if (output == "pdf") {
    axis_text_size <- 9 * text_scale     # prima 8
    axis_title_size <- 10 * text_scale
  } else {
    axis_text_size <- 11 * text_scale    # prima 10
    axis_title_size <- 12 * text_scale
  }

  theme_minimal() +
    theme(
      axis.ticks = element_line(),
      axis.ticks.length = unit(0.15, "cm"),
      axis.ticks.y = element_line(),
      axis.ticks.x = element_line(),
      axis.line = element_line(),
      axis.text = element_text(size = axis_text_size),
      axis.title = element_text(size = axis_title_size),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none")
}

## Define a custom color palette
my_palette <- c(
  blue_dark  = "#1f78b4",
  blue_light = "#a6cee3",
  red_dark   = "#e31a1c",
  red_light  = "#fca39c"
)

## Useful to plot cf_complex vs block

plot_ci_by_block <- function(param = NULL,
                             bootstrap_rus, bootstrap_ita,
                             params_rus, params_ita,
                             y_label = NULL,
                             return_data = FALSE,
                             add_smooth = FALSE,
                             theme_output = "pdf",
                             text_scale = 1) {

  blue_dark  <- "#1f78b4"
  blue_light <- "#a6cee3"
  red_dark   <- "#e31a1c"
  red_light  <- "#fca39c"

  # Helper function for CI computation and merge
  compute_ci <- function(bootstrap, params, param, source_label) {
    bootstrap %>%
      group_by(window_id) %>%
      summarise(
        lower = quantile(.data[[param]], 0.025),
        upper = quantile(.data[[param]], 0.975),
        .groups = "drop"
      ) %>%
      left_join(params, by = "window_id") %>%
      dplyr::select(source, window_id, !!sym(param), lower, upper)
  }

  ci_rus <- compute_ci(bootstrap_rus, params_rus, param, "Russian")
  ci_ita <- compute_ci(bootstrap_ita, params_ita, param, "Italian")
  ci_combined <- bind_rows(ci_rus, ci_ita)

  # Build the plot
  p <- ggplot(ci_combined, aes(x = window_id, y = .data[[param]], color = source, fill = source)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA) +
    geom_point(size = 1) +
    {if (add_smooth) geom_smooth(method = "lm", se = FALSE, linewidth = 1)} +
    scale_color_manual(values = c(Russian = red_dark, Italian = blue_dark)) +
    scale_fill_manual(values = c(Russian = red_light, Italian = blue_light)) +
    scale_x_continuous(breaks = seq(1, max(ci_combined$window_id), 1)) +
    labs(x = "Block", y = y_label, color = "Source", fill = "Source") +
    my_theme(theme_output, text_scale)

  if (return_data) {
    return(list(plot = p, data = ci_combined))
  } else {
    return(p)
  }
}
