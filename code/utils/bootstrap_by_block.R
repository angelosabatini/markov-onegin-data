## ============================================================
## File: bootstrap_by_block.R
##
## Purpose:
##   Perform moving block bootstrap on V/C sequences within
##   fixed-size windows.
##
## Functionality:
##   - Window-wise resampling
##   - Re-analysis via analyze_vc_sequence()
##
## Notes:
##   - Preserves local dependence structure
##   - Original (non-global) bootstrap implementation
## ============================================================

bootstrap_by_block <- function(vc_sequence, window_size = 10000,
                               block_size = 250, n_resamples = 1000,
                               source_label = "Empirical") {

  # Segment full sequence into non-overlapping blocks
  num_blocks <- floor(length(vc_sequence) / window_size)
  block_list <- split(vc_sequence[1:(num_blocks * window_size)],
                      rep(1:num_blocks, each = window_size))

  # Resampling with moving blocks and applying analysis
  bootstrap_summary <- purrr::map_dfr(1:n_resamples, function(i) {
    purrr::map_dfr(seq_along(block_list), function(w_id) {

      block <- block_list[[w_id]]
      n_subblocks <- ceiling(window_size / block_size)

      # Random starting points for subblocks
      starts <- sample.int(length(block) - block_size + 1, size = n_subblocks, replace = TRUE)
      resampled <- unlist(purrr::map(starts, ~ block[.x:(.x + block_size - 1)]))
      resampled <- resampled[1:window_size]  # truncate in case of overrun

      # Analyze resampled sequence
      result <- analyze_vc_sequence(resampled, window_sizes = window_size, source = source_label)

      result %>% mutate(resample_id = i, window_id = w_id)
    })
  })

  return(bootstrap_summary)
}
