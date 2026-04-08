## ============================================================
## Script: driver_trigrams_rus.R
##
## Purpose:
##   Extract vowel–consonant (V/C) trigrams from the Russian text of
##   *Evgenij Onegin*, preserving positional alignment with the original
##   character stream and providing short textual context for each match.
##
##   This driver performs:
##     - text cleaning and V/C encoding
##     - extraction of CCC, CCV, VVV, and VVC trigrams
##     - assignment of 10k-character block IDs
##     - verification of each trigram against the cleaned V/C stream
##
## Inputs:
##   - data/processed/Evgenij_Onegin_rus_only.RdS
##
## Dependencies:
##   - code/utils/setup.R
##   - code/utils/vc_alphabet.R
##   - code/utils/vc_encode.R
##   - code/utils/extract_match_context.R
##
## Outputs (processed):
##   - data/processed/df_trigrams_rus.RdS
##
## How to run:
##   - From repository root:
##       source("code/analysis/driver_trigrams_rus.R")
##
## Notes:
##   - Russian only. Italian trigrams are not extracted in this driver.
##   - The output is reused by:
##       * phonological probing (вст)
##       * candidate probe screening
##       * trigram trend visualization
##   - No statistical modeling is performed here.
## ============================================================

source("code/utils/setup.R")
source("code/utils/vc_alphabet.R")
source("code/utils/vc_encode.R")
source("code/utils/extract_match_context.R")

df <- readRDS("data/processed/Evgenij_Onegin_rus_only.RdS") %>%
  filter(type == "stanza", !is.na(stanza))

# Clean punctuation and format text
text <- str_c(df$text_rus, collapse = " ") %>%
  str_replace_all("\n", " ") %>%
  str_replace_all("[()\"…,.!?;:-]", " ") %>%
  str_squish()

# Lowercase version for encoding, original casing for context
char_seq  <- strsplit(text, "", fixed = TRUE)[[1]]
char_low  <- strsplit(str_to_lower(text), "", fixed = TRUE)[[1]]

# VC sequence aligned with char_seq: " " for non-letters, V/C for letters
vc_seq <- rep(" ", length(char_low))
vc_seq[char_low %in% vowels_all]     <- "V"
vc_seq[char_low %in% consonants_all] <- "C"

# Extract trigram matches, label each group and assign trigram labels
df_CCC <- extract_match_context(char_seq, vc_seq, pattern = "CCC") %>% mutate(trigram = "CCC")
df_CCV <- extract_match_context(char_seq, vc_seq, pattern = "CCV") %>% mutate(trigram = "CCV")
df_VVV <- extract_match_context(char_seq, vc_seq, pattern = "VVV") %>% mutate(trigram = "VVV")
df_VVC <- extract_match_context(char_seq, vc_seq, pattern = "VVC") %>% mutate(trigram = "VVC")

# Combine all trigram matches and assign 10k-character block IDs
all_matches <- bind_rows(df_CCC, df_CCV, df_VVV, df_VVC) %>%
  mutate(block_id = floor((vc_index_mid - 1) / 10000) + 1)

# Verify trigram sequence from cleaned V/C stream (excluding spaces)
vc_seq_nospace <- vc_seq[vc_seq != " "]
all_matches <- all_matches %>%
  mutate(trigram_check = map_chr(vc_index_mid, ~ paste(vc_seq_nospace[(.x - 1):(.x + 1)], collapse = "")))

saveRDS(all_matches, "data/processed/df_trigrams_rus.RdS")
