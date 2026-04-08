## ============================================================
## File: vc_alphabet.R
##
## Purpose:
##   Define vowel and consonant sets for V/C encoding across
##   Cyrillic (Russian) and Latin alphabets.
##
## Contents:
##   - Vowel sets
##   - Consonant sets
##   - Characters to exclude from encoding
##
## Notes:
##   - Language-agnostic by design
##   - Used by vc_encode() and all downstream analyses
## ============================================================

# Extended alphabet
vowels_cyr <- c("а", "э", "ы", "у", "о", "я", "е", "ё", "ю", "и", "й")
consonants_cyr <- c("б", "в", "г", "д", "ж", "з", "к", "л", "м", "н",
                    "п", "р", "с", "т", "ф", "х", "ц", "ч", "ш", "щ")
exclude_cyr <- c("ь", "ъ")

vowels_lat <- c("a", "à", "á", "â", "ä", "e", "è", "é", "ê", "ë", "i", "ì", "í", "î", "ï",
                "o", "ò", "ó", "ô", "ö", "u", "ù", "ú", "û", "ü")
consonants_lat <- c("b", "c", "d", "f", "g", "h", "j", "k", "l", "m", "n",
                    "p", "q", "r", "s", "t", "v", "w", "x", "y", "z",
                    "č", "š", "ž", "ỳ")
exclude_lat <- "'"

vowels_all     <- unique(c(vowels_cyr, vowels_lat))
consonants_all <- unique(c(consonants_cyr, consonants_lat))
exclude_all    <- unique(c(exclude_cyr, exclude_lat))
