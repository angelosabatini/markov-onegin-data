## ============================================================
## File: bootstrap_globally.R
##
## Purpose:
##   Experimental or legacy bootstrap routines operating on
##   full V/C sequences without window segmentation.
##
## Notes:
##   - Not used in final analysis pipeline
##   - Retained for completeness and reproducibility
## ============================================================

bootstrap_globally <- function(vc_sequence, window_size = 10000,
                               block_size = 250, n_resamples = 1000,
                               source_label = "Surrogate") {

  n_total <- floor(length(vc_sequence) / window_size) * window_size

  # Full sequence to be used
  vc_trimmed <- vc_sequence[1:n_total]
  n_blocks <- n_total / window_size

  bootstrap_summary <- purrr::map_dfr(1:n_resamples, function(i) {

    n_subblocks <- ceiling(n_total / block_size)
    starts <- sample.int(length(vc_trimmed) - block_size + 1, size = n_subblocks, replace = TRUE)
    resampled_sequence <- unlist(purrr::map(starts, ~ vc_trimmed[.x:(.x + block_size - 1)]))
    resampled_sequence <- resampled_sequence[1:n_total]  # truncate if too long

    # Now split into window-sized blocks
    windows <- split(resampled_sequence, rep(1:n_blocks, each = window_size))

    purrr::map_dfr(seq_along(windows), function(w_id) {
      result <- analyze_vc_sequence(windows[[w_id]], window_sizes = window_size, source = source_label)
      result %>% mutate(resample_id = i, window_id = w_id)
    })
  })

  return(bootstrap_summary)
}
