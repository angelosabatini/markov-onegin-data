## ============================================================
## File: extract_match_context.R
##
## Purpose:
##   Extract short or extended textual contexts around matched
##   V/C trigram patterns (e.g., CCC, VVV).
##
## Functionality:
##   - Pattern matching on compressed V/C stream
##   - Mapping back to original character positions
##   - Context extraction bounded by whitespace
##
## Notes:
##   - Central to phonological probing
##   - Assumes alignment between character and V/C sequences
## ============================================================

extract_match_context <- function(char_seq, vc_seq, pattern = "CCC") {
  char_pos <- seq_along(char_seq)
  text_seq <- char_seq

  is_valid <- vc_seq != " "
  vc_ns    <- vc_seq[is_valid]
  char_ns  <- char_seq[is_valid]
  pos_ns   <- char_pos[is_valid]

  vc_string <- paste(vc_ns, collapse = "")
  n <- nchar(vc_string)

  starts <- which(
    sapply(1:(n - 2), function(i) substr(vc_string, i, i + 2) == pattern)
  )

  if (length(starts) == 0) return(tibble())

  results <- vector("list", length(starts))

  for (i in seq_along(starts)) {
    start_ns <- starts[i]
    mid_ns   <- start_ns + 1
    end_ns   <- start_ns + 2

    start_pos <- pos_ns[start_ns]
    mid_pos   <- pos_ns[mid_ns]
    end_pos   <- pos_ns[end_ns]

    match_start <- which(char_pos == start_pos)[1]
    match_mid   <- which(char_pos == mid_pos)[1]
    match_end   <- which(char_pos == end_pos)[1]

    left_space <- if (match_start == 1) 0 else
      suppressWarnings(max(which(vc_seq[1:(match_start - 1)] == " ")))
    if (is.infinite(left_space)) left_space <- 0

    right_space <- if (match_end == length(vc_seq)) length(vc_seq) + 1 else
      suppressWarnings(min(which(vc_seq[(match_end + 1):length(vc_seq)] == " ")) + match_end)
    if (is.infinite(right_space)) right_space <- length(vc_seq) + 1

    context_range <- (left_space + 1):right_space
    context_str <- stringr::str_squish(paste0(text_seq[context_range], collapse = ""))

    pattern_str <- str_to_lower(paste0(text_seq[c(start_pos, mid_pos, end_pos)], collapse = ""))

    results[[i]] <- tibble(
      match_id         = i,
      match_center_pos = mid_pos,
      vc_index_mid     = mid_ns,
      context          = context_str,
      pattern_str      = pattern_str
    )
  }

  bind_rows(results)
}
