## ============================================================
## Script: make_all_figures.R
##
## Purpose:
##   Regenerate all paper figures (Fig2–Fig7) from processed datasets.
##   This script does NOT recompute the full pipeline; it only reads
##   objects in data/processed and writes PNGs to results/figures.
##
## Inputs (processed):
##   - data/processed/foreign_counts_rus_by_part.RdS
##   - data/processed/bootstrap_rus.RdS
##   - data/processed/bootstrap_ita.RdS
##   - data/processed/params_rus.RdS
##   - data/processed/params_ita.RdS
##   (plus any other processed objects used by Fig4–Fig7)
##
## Dependencies:
##   - code/utils/setup.R
##   - code/utils/simulate_vc_chain.R        (only if Fig5 simulates)
##   - code/utils/analyze_vc_sequence.R      (only if Fig5 simulates)
##   - tidyverse, patchwork
##
## Outputs:
##   - results/figures/Fig2.png ... results/figures/Fig7.png
##
## How to run:
##   - From repository root:
##       source("code/analysis/make_all_figures.R")
##
## Notes:
##   - If a figure is already present and identical, overwriting is harmless.
## ============================================================

## --- Setup (libraries, theme, palette, helpers) ------------------------
source("code/utils/setup.R")

## --- Output directories -------------------------------------
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

## ============================================================
## FIGURE 2
## ============================================================

message("Generating Figure 2...")

df <- readRDS("data/processed/Evgenij_Onegin_rus_only.RdS") %>%
  filter(type == "stanza", !is.na(stanza))

foreign_counts <- df %>%
  dplyr::select(part, stanza, line, text_rus) %>%
  mutate(verse_lines = str_split(text_rus, "\n")) %>%
  unnest(verse_lines) %>%
  unnest_tokens(word, verse_lines, token = "words", to_lower = FALSE) %>%
  filter(str_detect(word, "[\\p{Latin}]"), nchar(word) > 3) %>%
  count(part, name = "n_foreign")

p2 <- ggplot(foreign_counts, aes(x = factor(part), y = n_foreign)) +
  geom_col(fill = my_palette["blue_dark"]) +
  labs(
    x = "Part",
    y = "Count of Latin-script words (length > 3)"
  ) +
  my_theme(output = "png", text_scale = 1.2)

ggsave(
  filename = "results/figures/Fig2.png",
  plot = p2,
  width = 6, height = 4.5, units = "in",
  dpi = 300, bg = "white"
)

## ============================================================
## FIGURE 3
## ============================================================

message("Generating Figure 3...")

bootstrap_rus <- readRDS("data/processed/bootstrap_rus.RdS")
bootstrap_ita <- readRDS("data/processed/bootstrap_ita.RdS")
params_rus    <- readRDS("data/processed/params_rus.RdS")
params_ita    <- readRDS("data/processed/params_ita.RdS")

bootstrap_rus <- bootstrap_rus %>% mutate(q_00 = 1 - p_00, MD = 1 - cf_complex)
bootstrap_ita <- bootstrap_ita %>% mutate(q_00 = 1 - p_00, MD = 1 - cf_complex)
params_rus    <- params_rus    %>% mutate(q_00 = 1 - p_00, MD = 1 - cf_complex)
params_ita    <- params_ita    %>% mutate(q_00 = 1 - p_00, MD = 1 - cf_complex)

theme_output <- "png"
text_scale   <- 0.9

p_p  <- plot_ci_by_block("p_stationary", bootstrap_rus, bootstrap_ita,
                         params_rus, params_ita,
                         y_label = expression(p),
                         theme_output = theme_output,
                         text_scale = text_scale)

p_0  <- plot_ci_by_block("p_0", bootstrap_rus, bootstrap_ita,
                         params_rus, params_ita,
                         y_label = expression(p[0] == P(V~"|"~C)),
                         theme_output = theme_output,
                         text_scale = text_scale)

p_1  <- plot_ci_by_block("p_1", bootstrap_rus, bootstrap_ita,
                         params_rus, params_ita,
                         y_label = expression(p[1] == P(V~"|"~V)),
                         theme_output = theme_output,
                         text_scale = text_scale)

q_00 <- plot_ci_by_block("q_00", bootstrap_rus, bootstrap_ita,
                         params_rus, params_ita,
                         y_label = expression(q["00"] == P(C~"|"~CC)),
                         theme_output = theme_output,
                         text_scale = text_scale)

p_11 <- plot_ci_by_block("p_11", bootstrap_rus, bootstrap_ita,
                         params_rus, params_ita,
                         y_label = expression(p["11"] == P(V~"|"~VV)),
                         theme_output = theme_output,
                         text_scale = text_scale)

p_md <- plot_ci_by_block("MD", bootstrap_rus, bootstrap_ita,
                         params_rus, params_ita,
                         y_label = expression(MD),
                         theme_output = theme_output,
                         text_scale = text_scale)

p3 <- ((p_0 | p_1) / (p_11 | q_00)) / (p_p | p_md)

ggsave(
  filename = "results/figures/Fig3.png",
  plot = p3,
  width = 6, height = 6, units = "in",
  dpi = 300, bg = "white"
)

## ============================================================
## FIGURE 4
## ============================================================

message("Generating Figure 4...")

# Russian (window_id == 1)
bootstrap_rus_w1 <- readRDS("data/processed/bootstrap_rus.RdS") %>%
  filter(window_id == 1)

params_rus_w1 <- readRDS("data/processed/params_rus.RdS") %>%
  filter(window_id == 1)

acf_summary_rus <- bootstrap_rus_w1 %>%
  pivot_longer(cols = starts_with("acf"), names_to = "lag", values_to = "acf_value") %>%
  group_by(lag) %>%
  summarise(
    lower = quantile(acf_value, 0.025, na.rm = TRUE),
    upper = quantile(acf_value, 0.975, na.rm = TRUE),
    mean  = mean(acf_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(lag = as.integer(stringr::str_remove(lag, "acf")))

acf_empirical_rus <- params_rus_w1 %>%
  slice(1) %>%
  dplyr::select(starts_with("acf")) %>%
  pivot_longer(cols = everything(), names_to = "lag", values_to = "empirical_value") %>%
  mutate(lag = as.integer(stringr::str_remove(lag, "acf")))

acf_combined_rus <- left_join(acf_summary_rus, acf_empirical_rus, by = "lag")

n <- 10000
conf_int <- 1.96 / sqrt(n)

p4_rus <- ggplot(acf_combined_rus, aes(x = lag)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "gray40") +
  geom_point(aes(y = empirical_value), color = my_palette["red_dark"], size = 1) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_hline(yintercept = c(-conf_int, conf_int), linetype = "dashed",
             color = my_palette["red_light"]) +
  labs(x = "Lag", y = "ACF") +
  annotate("text", x = Inf, y = Inf, label = "(A)", hjust = 1.2, vjust = 1.2, size = 3.5) +
  scale_y_continuous(
    breaks = seq(-0.6, 0.2, by = 0.2),
    limits = c(-0.6, 0.2),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  my_theme(output = "png", text_scale = 0.8)

# Italian (window_id == 1)
bootstrap_ita_w1 <- readRDS("data/processed/bootstrap_ita.RdS") %>%
  filter(window_id == 1)

params_ita_w1 <- readRDS("data/processed/params_ita.RdS") %>%
  filter(window_id == 1)

acf_summary_ita <- bootstrap_ita_w1 %>%
  pivot_longer(cols = starts_with("acf"), names_to = "lag", values_to = "acf_value") %>%
  group_by(lag) %>%
  summarise(
    lower = quantile(acf_value, 0.025, na.rm = TRUE),
    upper = quantile(acf_value, 0.975, na.rm = TRUE),
    mean  = mean(acf_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(lag = as.integer(stringr::str_remove(lag, "acf")))

acf_empirical_ita <- params_ita_w1 %>%
  slice(1) %>%
  dplyr::select(starts_with("acf")) %>%
  pivot_longer(cols = everything(), names_to = "lag", values_to = "empirical_value") %>%
  mutate(lag = as.integer(stringr::str_remove(lag, "acf")))

acf_combined_ita <- left_join(acf_summary_ita, acf_empirical_ita, by = "lag")

p4_ita <- ggplot(acf_combined_ita, aes(x = lag)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "gray40") +
  geom_point(aes(y = empirical_value), color = my_palette["blue_dark"], size = 1) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_hline(yintercept = c(-conf_int, conf_int), linetype = "dashed",
             color = my_palette["blue_light"]) +
  labs(x = "Lag", y = "ACF") +
  annotate("text", x = Inf, y = Inf, label = "(B)", hjust = 1.2, vjust = 1.2, size = 3.5) +
  scale_y_continuous(
    breaks = seq(-0.6, 0.2, by = 0.2),
    limits = c(-0.6, 0.2),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  my_theme(output = "png", text_scale = 0.8)

p4 <- p4_rus | p4_ita

ggsave(
  filename = "results/figures/Fig4.png",
  plot = p4,
  width = 6, height = 3, units = "in",
  dpi = 300, bg = "white"
)

## ============================================================
## FIGURE 5
## ============================================================

message("Generating Figure 5")

source("code/utils/simulate_vc_chain.R")
source("code/utils/analyze_vc_sequence.R")

results_empirical <- readRDS("data/processed/params_rus.RdS") %>%
  dplyr::select(p_stationary, p_00, p_01, p_10, p_11, cf_complex) %>%
  mutate(MD = 1 - cf_complex)

p_00 <- results_empirical$p_00[1]
p_01 <- results_empirical$p_01[1]
p_10 <- results_empirical$p_10[1]
p_11 <- results_empirical$p_11[1]
cf_complex <- results_empirical$cf_complex[1]
MD_rus <- results_empirical$MD[1]

# Build transition matrix
build_transition_matrix <- function(p_00, p_01, p_10, p_11) {
  matrix(c(
    p_11, 1 - p_11, 0,    0,
    0,    0,    p_10, 1 - p_10,
    p_01, 1 - p_01, 0,    0,
    0,    0,    p_00, 1 - p_00
  ), nrow = 4, byrow = TRUE)
}

# Build transition matrix
trans_matrix <- build_transition_matrix(p_00, p_01, p_10, p_11)

# Compute stationary distribution of a Markov chain
stationary_distribution <- function(P) {
  d <- nrow(P)
  A <- t(P) - diag(d)
  A[d, ] <- 1  # enforce sum(pi) = 1
  b <- c(rep(0, d - 1), 1)
  solve(A, b)
}
pi_stat <- stationary_distribution(trans_matrix)
names(pi_stat) <- c("VV", "VC", "CV", "CC")

# Simulate 500 sequences of 10000 characters
n_simulations <- 500
set.seed(123)
simulated_results_list <- map(1:n_simulations, ~{
  sim_seq <- simulate_vc_chain(n = 10000, trans_matrix = trans_matrix, init_probs = pi_stat)
  analyze_vc_sequence(sim_seq, window_sizes = c(10000), source = paste0("Simulated_", .x))
})

# Combine results
results_simulated <- bind_rows(simulated_results_list)

results_simulated_rus <- results_simulated %>%
  mutate(MD = 1 - cf_complex)

results_empirical <- readRDS("data/processed/params_ita.RdS") %>%
  dplyr::select(p_stationary, p_00, p_01, p_10, p_11, cf_complex) %>%
  mutate(MD = 1 - cf_complex)

p_00 <- results_empirical$p_00[1]
p_01 <- results_empirical$p_01[1]
p_10 <- results_empirical$p_10[1]
p_11 <- results_empirical$p_11[1]
cf_complex <- results_empirical$cf_complex[1]
MD_ita <- results_empirical$MD[1]

# Build transition matrix
build_transition_matrix <- function(p_00, p_01, p_10, p_11) {
  matrix(c(
    p_11, 1 - p_11, 0,    0,
    0,    0,    p_10, 1 - p_10,
    p_01, 1 - p_01, 0,    0,
    0,    0,    p_00, 1 - p_00
  ), nrow = 4, byrow = TRUE)
}

# Build transition matrix
trans_matrix <- build_transition_matrix(p_00, p_01, p_10, p_11)

# Compute stationary distribution of a Markov chain
stationary_distribution <- function(P) {
  d <- nrow(P)
  A <- t(P) - diag(d)
  A[d, ] <- 1  # enforce sum(pi) = 1
  b <- c(rep(0, d - 1), 1)
  solve(A, b)
}
pi_stat <- stationary_distribution(trans_matrix)
names(pi_stat) <- c("VV", "VC", "CV", "CC")

# Simulate 500 sequences of 10000 characters
n_simulations <- 500
set.seed(123)
simulated_results_list <- map(1:n_simulations, ~{
  sim_seq <- simulate_vc_chain(n = 10000, trans_matrix = trans_matrix, init_probs = pi_stat)
  analyze_vc_sequence(sim_seq, window_sizes = c(10000), source = paste0("Simulated_", .x))
})

# Combine results
results_simulated <- bind_rows(simulated_results_list)

results_simulated_ita <- results_simulated %>%
  mutate(MD = 1 - cf_complex)

results_simulated_ita <- results_simulated_ita %>% mutate(language = "Italian")
results_simulated_rus <- results_simulated_rus %>% mutate(language = "Russian")
results_simulated_tot <- bind_rows(results_simulated_rus, results_simulated_ita)

results_simulated_tot <- results_simulated_tot %>%
  mutate(language = as.factor(as.character(language))) %>%  # forza character → factor pulito
  mutate(language = factor(language, levels = c("Russian", "Italian")))  # impone ordine

print(unique(results_simulated_tot$language))  # controllo finale

blue_dark  <- "#1f78b4"
blue_light <- "#a6cee3"
red_dark   <- "#e31a1c"
red_light  <- "#fca39c"

p5 <- results_simulated_tot %>%
  ggplot(aes(x = MD, fill = language)) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = c("Russian" = red_light, "Italian" = blue_light)) +
  geom_vline(xintercept = MD_rus, color = red_light, linetype = "dashed", size = 1) +
  geom_vline(xintercept = MD_ita, color = blue_light, linetype = "dashed", size = 1) +
  scale_x_continuous(limits = c(0.7, 0.9)) +
  labs(x = "MD", y = "Density", fill = "Source") +
  my_theme(output = "png", text_scale = 1.2)

ggsave(
  filename = "results/figures/Fig5.png",
  plot = p5,
  width = 6, height = 4.5, units = "in",
  dpi = 300, bg = "white"
)

## ============================================================
## FIGURE 6
## ============================================================

message("Generating Figure 6...")

bootstrap_rus <- readRDS("data/processed/bootstrap_rus.RdS")
bootstrap_ita <- readRDS("data/processed/bootstrap_ita.RdS")

combined_bootstrap_all <- bind_rows(bootstrap_rus, bootstrap_ita) %>%
  mutate(MD = 1 - cf_complex) %>%
  dplyr::select(resample_id, block = window_id, MD, source)

bootstrap_params <- combined_bootstrap_all %>%
  group_by(resample_id) %>%
  summarise(model = list(lm(MD ~ block * source, data = pick(everything()))), .groups = "drop") %>%
  mutate(coefs = purrr::map(model, ~ broom::tidy(.x, conf.int = TRUE))) %>%
  tidyr::unnest(coefs)

bootstrap_params_summary <- bootstrap_params %>%
  group_by(term) %>%
  summarise(
    mean_estimate = mean(estimate),
    ci_lower = quantile(estimate, 0.025),
    ci_upper = quantile(estimate, 0.975),
    .groups = "drop"
  )

vlines_df <- bootstrap_params_summary %>%
  rename(mean = mean_estimate, lower = ci_lower, upper = ci_upper)

p6 <- ggplot(bootstrap_params, aes(x = estimate, fill = term)) +
  geom_density(alpha = 0.2, color = "black") +
  geom_vline(data = vlines_df, aes(xintercept = mean),
             color = my_palette["blue_dark"], linetype = "dashed", linewidth = 1) +
  geom_vline(data = vlines_df, aes(xintercept = lower),
             color = my_palette["red_dark"], linetype = "dotted") +
  geom_vline(data = vlines_df, aes(xintercept = upper),
             color = my_palette["red_dark"], linetype = "dotted") +
  facet_wrap(~ term, scales = "free") +
  labs(x = "Estimate", y = "Density") +
  my_theme(output = "png", text_scale = 0.8)

ggsave(
  filename = "results/figures/Fig6.png",
  plot = p6,
  width = 6, height = 4, units = "in",
  dpi = 300, bg = "white"
)

## ============================================================
## FIGURE 7
## ============================================================

message("Generating Figure 7...")

bootstrap_rus <- readRDS("data/processed/bootstrap_rus.RdS") %>%
  mutate(q_00 = 1 - p_00, MD = 1 - cf_complex)

bootstrap_ita <- readRDS("data/processed/bootstrap_ita.RdS") %>%
  mutate(q_00 = 1 - p_00, MD = 1 - cf_complex)

params_rus <- readRDS("data/processed/params_rus.RdS") %>%
  mutate(q_00 = 1 - p_00, MD = 1 - cf_complex)

params_ita <- readRDS("data/processed/params_ita.RdS") %>%
  mutate(q_00 = 1 - p_00, MD = 1 - cf_complex)

p7 <- plot_ci_by_block(
  "MD",
  bootstrap_rus, bootstrap_ita,
  params_rus, params_ita,
  y_label = expression(MD),
  add_smooth = TRUE,
  theme_output = "png",
  text_scale = 1.2
)

ggsave(
  filename = "results/figures/Fig7.png",
  plot = p7,
  width = 6, height = 4, units = "in",
  dpi = 300, bg = "white"
)

## ============================================================
## DONE
## ============================================================

message("All figures generated successfully.")
