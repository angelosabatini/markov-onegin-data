## ============================================================
## Script: driver_phonological_probing.R
##
## Purpose:
##   Reproduce the phonological probe analysis for the Russian trigram "вст":
##     (1) Load trigram matches with context (Russian only)
##     (2) UDPipe annotation of contexts (lemma-level)
##     (3) Topic labeling (encounter / emotion / other)
##     (4) Block-level aggregation and trend tests vs block index
##     (5) Join with Memory Depth (MD) from params_rus
##
## Inputs (processed):
##   - data/processed/df_trigrams_rus.RdS
##   - data/processed/params_rus.RdS
## Optional input:
##   - data/processed/candidate_probes_rus.RdS
##
## Dependencies:
##   - udpipe (requires local model file)
##   - code/utils/setup.R (optional: theme/helpers)
##   - code/utils/ (if reusing helper functions)
##
## External resources:
##   - model/russian-gsd-ud-2.5-191206.udpipe
##
## Outputs (processed):
##   - data/processed/df_corr_vst_rus.RdS
##   (and/or printed summary tables for the paper)
##
## How to run:
##   - From repository root:
##       source("code/analysis/driver_phonological_probing.R")
##
## Notes:
##   - Russian only (the phonological probe analysis is not run on Italian).
##   - The topic mapping is rule-based and mirrors the paper.
## ============================================================


source("code/utils/setup.R")

## --- inputs --------------------------------------------------
target_patterns <- NULL   # e.g. "вст"; if NULL -> candidate list only
min_n <- 2                # min occurrences per block to keep a (trigram, pattern_str, block)
top_n <- 30               # how many patterns to screen per trigram type (CCC/CCV/VVV/VVC)

## --- load data ----------------------------------------------
df_trigrams <- readRDS("data/processed/df_trigrams_rus.RdS") %>%
  mutate(context = str_to_lower(context))

## --- helper: trigram_visualizer --------
trigram_visualizer <- function(df_trigram,
                               min_n = 5,
                               top_n = 20,
                               top_trigrams = NULL,
                               return_data = TRUE) {

  context_freq <- df_trigram %>%
    count(trigram, pattern_str, block_id) %>%
    filter(n >= min_n)

  totals <- context_freq %>%
    group_by(trigram, pattern_str) %>%
    summarise(total = sum(n), .groups = "drop")

  if (!is.null(top_trigrams)) {
    stopifnot(is.character(top_trigrams))
    totals <- totals %>% filter(pattern_str %in% top_trigrams)
  } else {
    totals <- totals %>%
      arrange(desc(total)) %>%
      group_by(trigram) %>%
      slice_head(n = top_n)
  }

  plot_data <- context_freq %>%
    semi_join(totals, by = c("trigram", "pattern_str")) %>%
    left_join(totals, by = c("trigram", "pattern_str")) %>%
    mutate(pattern_str = fct_reorder(pattern_str, total))

  if (return_data) return(list(data = plot_data)) else return(plot_data)
}

## --- helper: trend test (Spearman + p) -----------------------
trend_spearman <- function(block_id, counts) {
  if (length(unique(counts)) <= 1) {
    return(tibble(rho = NA_real_, p_value = NA_real_))
  }
  ct <- suppressWarnings(cor.test(block_id, counts, method = "spearman", exact = FALSE))
  tibble(rho = unname(ct$estimate), p_value = ct$p.value)
}

## --- 1) Candidate screening data (counts per block) ----------
viz <- trigram_visualizer(
  df_trigram = df_trigrams,
  min_n = min_n,
  top_n = top_n,
  top_trigrams = if (is.null(target_patterns)) NULL else target_patterns,
  return_data = TRUE
)

screen_data <- viz$data %>%
  dplyr::select(trigram, pattern_str, block_id, n)

## --- 2) Candidate list (Tab4 backbone) -----------------------
## build per (trigram, pattern_str) a full vector of counts across blocks
all_blocks <- sort(unique(df_trigrams$block_id))

cand <- screen_data %>%
  group_by(trigram, pattern_str) %>%
  summarise(
    total = sum(n),
    counts = list({
      tmp <- tibble(block_id = all_blocks) %>%
        left_join(pick(everything()) %>% dplyr::select(block_id, n), by = "block_id") %>%
        mutate(n = replace_na(n, 0)) %>%
        arrange(block_id)
      tmp$n
    }),
    .groups = "drop"
  ) %>%
  mutate(
    test = purrr::map(counts, ~ trend_spearman(all_blocks, .x)),
    rho = purrr::map_dbl(test, "rho"),
    p_value = purrr::map_dbl(test, "p_value"),
    direction = case_when(
      is.na(rho) ~ NA_character_,
      rho > 0 ~ "increasing",
      rho < 0 ~ "decreasing",
      TRUE ~ "flat"
    ),
    pattern = case_when(
      trigram %in% c("CCC", "VVV") ~ "persistent",
      trigram %in% c("CCV", "VVC") ~ "alternating",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(p_value)

# shortlist: e.g. p < 0.01 and total >= 5
short <- cand %>% filter(!is.na(p_value), p_value < 0.01, total >= 5)

## ============================================================
## 3) Probe-specific UDPIPE block (Tab5 style) for target pattern
## ============================================================

df_focus <- df_trigrams %>%
  filter(pattern_str %in% "вст") %>%
  mutate(match_id = row_number())

## Load MD (block-wise)
MD_df <- readRDS("data/processed/params_rus.RdS") %>%
  mutate(MD = 1 - cf_complex) %>%
  transmute(block_id = window_id, MD)

## UDPIPE annotate contexts
ud_model <- udpipe_load_model("model/russian-gsd-ud-2.5-191206.udpipe")

anno <- udpipe_annotate(ud_model, x = df_focus$context, doc_id = df_focus$match_id)
lemma_df <- as_tibble(anno) %>%
  mutate(match_id = as.integer(doc_id)) %>%
  left_join(df_focus %>% dplyr::select(match_id, block_id), by = "match_id") %>%
  mutate(context_id = paste(block_id, match_id, sep = "_"))

## keep monoword contexts only
lemma_df_mono <- lemma_df %>%
  group_by(context_id) %>%
  #filter(n() == 1) %>%
  ungroup()

## thematic coding (your exact rules)
lemma_df_mono <- lemma_df_mono %>%
  mutate(
    lemma_lower = str_to_lower(lemma),
    topic = case_when(
      str_detect(lemma_lower, "встр") ~ "encounter",
      str_detect(lemma_lower, "вступ") ~ "encounter",
      str_detect(lemma_lower, "здравст") ~ "encounter",
      str_detect(lemma_lower, "бесчувст") ~ "emotion",
      str_detect(lemma_lower, "вста") ~ "emotion",
      str_detect(lemma_lower, "девств") ~ "emotion",
      str_detect(lemma_lower, "предчувств") ~ "emotion",
      str_detect(lemma_lower, "чувст") ~ "emotion",
      TRUE ~ "other"
    )
  )

## block-wise counts including other
df_topic_block <- lemma_df_mono %>%
  group_by(block_id) %>%
  summarise(
    emotion = sum(topic == "emotion"),
    encounter = sum(topic == "encounter"),
    occurrence = n(), # includes OTHER
    .groups = "drop"
  ) %>%
  left_join(MD_df, by = "block_id") %>%
  arrange(block_id)

## --- Spearman trend tests vs block_id ------------------------

x <- df_topic_block$block_id

trend_table <- bind_rows(
  trend_spearman(x, df_topic_block$emotion) %>%
    mutate(category = "emotion", n = sum(df_topic_block$emotion)),

  trend_spearman(x, df_topic_block$encounter) %>%
    mutate(category = "encounter", n = sum(df_topic_block$encounter)),

  trend_spearman(x, df_topic_block$emotion + df_topic_block$encounter) %>%
    mutate(category = "emotion + encounter", n = sum(df_topic_block$emotion + df_topic_block$encounter)),
) %>%
  dplyr::select(category, n, rho, p_value)

trend_table

saveRDS(df_topic_block, "data/processed/df_corr_vst_rus.RdS")

message("Probe-specific UDPIPE summary saved.")
