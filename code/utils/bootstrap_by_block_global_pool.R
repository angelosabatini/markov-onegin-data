## ============================================================
## File: bootstrap_by_block_global_pool.R
##
## Purpose:
##   Perform moving block bootstrap using a global pool of
##   sub-blocks across the entire V/C sequence.
##
## Functionality:
##   - Resampling from shared block pool
##   - Improved stability for edge blocks
##
## Notes:
##   - Used for final bootstrap results in the paper
##   - Supersedes bootstrap_by_block.R for main analyses
## ============================================================

bootstrap_by_block_global_pool <- function(vc_sequence,
                                           window_size = 10000,
                                           block_size = 250,
                                           n_resamples = 1000,
                                           source_label = "Empirical") {

  # one seed outside, not here (driver sets it)
  n <- length(vc_sequence)
  num_blocks <- floor(n / window_size)          # full 10k blocks
  tail_size  <- n - num_blocks * window_size    # residual (e.g., 7000)

  # window sizes: full blocks + optional tail block
  win_sizes <- c(rep(window_size, num_blocks), if (tail_size > 0) tail_size else NULL)

  purrr::map_dfr(1:n_resamples, function(i) {
    purrr::map_dfr(seq_along(win_sizes), function(w_id) {

      wlen <- win_sizes[[w_id]]
      n_subblocks <- ceiling(wlen / block_size)

      # starts anywhere in the FULL sequence (global pool)
      starts <- sample.int(n, size = n_subblocks, replace = TRUE)

      # circular block extraction from global pool
      resampled <- unlist(purrr::map(starts, function(s) {
        idx <- ((s:(s + block_size - 1)) - 1) %% n + 1
        vc_sequence[idx]
      }), use.names = FALSE)

      resampled <- resampled[1:wlen]

      analyze_vc_sequence(resampled, window_sizes = wlen, source = source_label) %>%
        mutate(resample_id = i, window_id = w_id)
    })
  })
}
