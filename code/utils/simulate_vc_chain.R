## ============================================================
## File: simulate_vc_chain.R
##
## Purpose:
##   Simulate V/C sequences from a specified Markov transition
##   matrix and initial state distribution.
##
## Functionality:
##   - Supports 4-state Markov chains
##   - Used for null and comparison models
##
## Notes:
##   - No estimation is performed here
##   - Used in Figure 5 simulations
## ============================================================

simulate_vc_chain <- function(n, trans_matrix, init_probs = NULL, labels = c("VV", "VC", "CV", "CC")) {

  n_states <- length(labels)
  if (is.null(init_probs)) {
    init_probs <- rep(1 / n_states, n_states)
  }

  sim_states <- numeric(n)
  sim_states[1] <- sample(1:n_states, size = 1, prob = init_probs)

  for (i in 2:n) {
    sim_states[i] <- sample(1:n_states, size = 1, prob = trans_matrix[sim_states[i - 1], ])
  }

  label_sequence <- labels[sim_states]
  first_char <- substr(label_sequence[1], 1, 1)
  second_chars <- substr(label_sequence, 2, 2)

  vc_sequence <- c(first_char, second_chars)
  return(vc_sequence)
}
