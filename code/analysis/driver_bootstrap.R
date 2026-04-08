## ============================================================
## Script: driver_bootstrap.R
##
## Purpose:
##   Illustrate the core quantitative pipeline on the V/C sequences:
##     (1) V/C encoding (Russian text) -> vc_sequence_rus
##     (2) Markov summary statistics by 10k blocks -> params_rus
##         including the final residual block (~7k) as an extra window
##     (3) Moving-block bootstrap by window (Russian + Italian) -> bootstrap_*
##
## Inputs:
##   - data/processed/Evgenij_Onegin_rus_only.RdS
##   - data/processed/Evgenij_Onegin_vc_sequences.RdS   (optional shortcut)
##   - code/utils/vc_alphabet.R
##   - code/utils/vc_encode.R
##   - code/utils/analyze_vc_sequence.R
##   - code/utils/bootstrap_by_block_global_pool.R
##
## Outputs (processed):
##   - data/processed/vc_sequence_rus.RdS        (optional)
##   - data/raw/vc_sequence_rus.txt              (optional human-readable)
##   - data/processed/params_rus.RdS
##   - data/processed/bootstrap_rus.RdS
##   - data/processed/bootstrap_ita.RdS
##
## How to run:
##   - From repository root:
##       source("code/analysis/driver_bootstrap.R")
##
## Notes:
##   - This driver is intentionally explicit and self-contained.
##   - Long-running step: bootstrap (n_resamples = 1000).
##   - Russian sequence includes a final residual window (< 10k chars),
##     stored as an extra window_id in params_rus.
## ============================================================

source("code/utils/setup.R")
source("code/utils/vc_alphabet.R")
source("code/utils/vc_encode.R")
source("code/utils/analyze_vc_sequence.R")
source("code/utils/bootstrap_by_block_global_pool.R")

# ============================================================
# RUS: build VC from rus-only dataset (no epigrafe)
# ============================================================

df_rus <- readRDS("data/processed/Evgenij_Onegin_rus_only.RdS") %>%
  filter(type != "epigrafe")

vc_vec_rus <- df_rus$text_rus %>%
  paste(collapse = "") %>%
  vc_encode(
    vowels = vowels_all,
    consonants = consonants_all,
    exclude_chars = exclude_all,
    preserve_case = FALSE
  )

# ============================================================
# PARAMS RUS (10k blocks + final residual block)
# ============================================================

params_rus <- analyze_vc_sequence(
  vc_sequence = vc_vec_rus[1:100000],
  window_sizes = 10000,
  source = "Russian"
)
params_tail <- analyze_vc_sequence(
    vc_sequence = vc_vec_rus[100001:107000],
    window_sizes = 7000,
    source = "Russian"
  ) %>%
    mutate(window_id = 11)

params_rus <- bind_rows(params_rus, params_tail)

saveRDS(params_rus, "data/processed/params_rus.RdS")

# ============================================================
# ITA: take VC directly from vc_sequences dataset (no text)
# ============================================================

df_vc <- readRDS("data/processed/Evgenij_Onegin_vc_sequences.RdS")

vc_vec_ita <- unlist(df_vc$vc_seq_ita)

params_ita <- analyze_vc_sequence(vc_vec_ita[1:120000], window_sizes = 10000, source = "Italian")
saveRDS(params_ita, "data/processed/params_ita.RdS")

# ============================================================
# BOOTSTRAP (global donor pool; rus includes tail automatically)
# ============================================================

set.seed(123)
bootstrap_rus <- bootstrap_by_block_global_pool(
  vc_sequence  = vc_vec_rus[1:107000],
  window_size  = 10000,
  block_size   = 250,
  n_resamples  = 1000,
  source_label = "Russian"
)

set.seed(123)
bootstrap_ita <- bootstrap_by_block_global_pool(
  vc_sequence  = vc_vec_ita[1:120000],
  window_size  = 10000,
  block_size   = 250,
  n_resamples  = 1000,
  source_label = "Italian"
)

saveRDS(bootstrap_rus, "data/processed/bootstrap_rus.RdS")
saveRDS(bootstrap_ita, "data/processed/bootstrap_ita.RdS")
