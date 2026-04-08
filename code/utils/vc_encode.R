## ============================================================
## File: vc_encode.R
##
## Purpose:
##   Encode a character string into a vowel/consonant (V/C)
##   symbolic sequence using predefined alphabets.
##
## Functionality:
##   - Preserves or normalizes case
##   - Filters excluded characters
##   - Returns a clean V/C sequence
##
## Notes:
##   - Stateless and deterministic
##   - No text segmentation is performed here
## ============================================================

vc_encode <- function(text,
                      vowels,
                      consonants,
                      exclude_chars = NULL,
                      preserve_case = TRUE) {

  stopifnot(is.character(text))

  # Normalize text
  txt <- if (preserve_case) text else str_to_lower(text)

  # Tokenize characters
  tokens <- unlist(tokenize_character_shingles(txt, n = 1, n_min = 1))

  # Filter valid letters
  valid <- c(vowels, consonants)

  if (!is.null(exclude_chars)) {
    tokens <- tokens[!tokens %in% exclude_chars]
  }
  tokens <- tokens[tokens %in% valid]

  # Encode V/C
  vc <- case_when(
    tokens %in% vowels ~ "V",
    tokens %in% consonants ~ "C",
    TRUE ~ NA_character_
  )

  vc[!is.na(vc)]
}
