## ============================================================
## File: analyze_vc_sequence.R
##
## Purpose:
##   Quantitative analysis of V/C sequences over fixed-length
##   windows, including autocorrelation, transition probabilities,
##   and dispersion correction factors.
##
## Functionality:
##   - 2-state and 4-state Markov statistics
##   - ACF and Ljung–Box diagnostics
##   - Computation of cf_simple and cf_complex
##
## Notes:
##   - Core analytical engine of the paper
##   - Does not perform simulation or bootstrap
## ============================================================

analyze_vc_sequence <- function(vc_sequence, window_sizes = c(5000, 10000, 20000), source = "Empirical") {

  safe_divide <- function(a, b) ifelse(b == 0, NA_real_, a / b)

  results <- list()

  for (window_size in window_sizes) {

    num_windows <- floor(length(vc_sequence) / window_size)
    if (num_windows < 1) next

    windows <- split(vc_sequence[1:(num_windows * window_size)],
                     rep(1:num_windows, each = window_size))

    analysis <- purrr::map_dfr(seq_along(windows), ~{

      window_seq <- windows[[.x]]
      binom_seq <- if_else(window_seq == "V", 1, 0)
      n <- length(binom_seq)

      acf_result <- acf(binom_seq, lag.max = 10, plot = FALSE)
      acf_values <- head(as.numeric(acf_result$acf[-1]), 5)
      if (length(acf_values) < 5) {
        acf_values <- c(acf_values, rep(NA, 5 - length(acf_values)))
      }

      ljung_box <- Box.test(binom_seq, lag = 10, type = "Ljung-Box")
      lb_pvalue <- ljung_box$p.value

      p_stationary <- mean(binom_seq == 1)
      q_stationary <- 1 - p_stationary

      bigrams <- tokenize_character_shingles(paste(window_seq, collapse = ""),
                                             lowercase = FALSE, n = 2, n_min = 2)[[1]]
      trigrams <- tokenize_character_shingles(paste(window_seq, collapse = ""),
                                              lowercase = FALSE, n = 3, n_min = 3)[[1]]

      n_VV <- sum(bigrams == "VV")
      n_VC <- sum(bigrams == "VC")
      n_CV <- sum(bigrams == "CV")
      n_CC <- sum(bigrams == "CC")

      p_1 <- safe_divide(n_VV, n_VV + n_VC)
      p_0 <- safe_divide(n_CV, n_CV + n_CC)

      n_VVV <- sum(trigrams == "VVV")
      n_VVC <- sum(trigrams == "VVC")
      n_VCV <- sum(trigrams == "VCV")
      n_VCC <- sum(trigrams == "VCC")
      n_CVV <- sum(trigrams == "CVV")
      n_CVC <- sum(trigrams == "CVC")
      n_CCV <- sum(trigrams == "CCV")
      n_CCC <- sum(trigrams == "CCC")

      p_11 <- safe_divide(n_VVV, n_VVV + n_VVC)
      p_10 <- safe_divide(n_VCV, n_VCV + n_VCC)
      p_01 <- safe_divide(n_CVV, n_CVV + n_CVC)
      p_00 <- safe_divide(n_CCV, n_CCC + n_CCV)

      d <- safe_divide(p_1 - p_0, 1)
      cf_simple <- safe_divide(1 + d, 1 - d)

      eta <- safe_divide(p_11 - p_1, 1 - p_1)
      nu  <- safe_divide(p_0 - p_00, p_0)

      cf_complex <- 0.5 * cf_simple * ((1 + eta) / (1 - eta) + (1 + nu) / (1 - nu)) +
        (q_stationary - p_stationary) * (nu - eta) / ((1 - eta) * (1 - nu))

      tibble(
        source = source,
        window_size = window_size,
        window_id = .x,
        p_stationary = p_stationary,
        q_stationary = q_stationary,
        acf1 = acf_values[1],
        acf2 = acf_values[2],
        acf3 = acf_values[3],
        acf4 = acf_values[4],
        acf5 = acf_values[5],
        lb_pvalue = lb_pvalue,
        p_0 = p_0,
        p_1 = p_1,
        p_00 = p_00,
        p_01 = p_01,
        p_10 = p_10,
        p_11 = p_11,
        cf_simple = cf_simple,
        cf_complex = cf_complex
      )
    })

    results[[as.character(window_size)]] <- analysis
  }

  bind_rows(results)
}
