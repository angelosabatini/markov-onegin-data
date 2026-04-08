## Repository overview

This repository accompanies the paper
"Markov reads Puškin, again: A statistical journey into the poetic world of Evgenij Onegin".

It exposes the full analytical pipeline used in the paper, from vowel–consonant (V/C) encoding 
to Markovian modeling, bootstrap validation, and phonological probing via contextual trigram analysis.

The repository is **not** intended as a general-purpose package, but as a **transparent, inspectable 
research compendium**.

## Folder structure

- data/processed  
  Processed datasets used in the paper (V/C sequences, Markov parameters, bootstrap outputs, trigram data).

- code/utils  
  Reusable functions for V/C encoding, Markov analysis, bootstrap, trigram extraction.

- code/analysis  
  Driver scripts illustrating each analytical step.

- results/figures  
  Final figures used in the paper.
  
## Conceptual pipeline

The analysis follows a layered structure:

1. **Text preparation**
   - Russian text of *Evgenij Onegin* is segmented at stanza level
   - Punctuation and non-letter characters are handled conservatively

2. **V/C encoding**
   - Characters are mapped to V (vowel) or C (consonant)
   - Encoding is language-agnostic and alphabet-driven

3. **Markovian modeling**
   - 2-state and 4-state Markov models are estimated
   - Transition probabilities, autocorrelation, and dispersion factors are computed
   - Analysis is performed over fixed-length blocks (10,000 characters)

4. **Bootstrap validation**
   - Moving block bootstrap with a global block pool
   - Used to assess stability of Markovian parameters and Memory Depth (MD)

5. **Trigram analysis**
   - Extraction of V/C trigrams (CCC, CCV, VVV, VVC)
   - Contextual reconstruction at character level
   - Block-wise frequency analysis

6. **Phonological probing**
   - Linguistic annotation via UDPipe (Russian)
   - Lemma-based thematic categorization
   - Trend and correlation analysis

The trigram-based phonological analysis (including the probe “вст”) is illustrated in:
> source("code/analysis/driver_phonological_probing.R")

This script uses trigram data computed by `driver_trigrams_rus.R` and performs linguistic 
annotation via UDPipe (Russian model provided in model/).

Only the **key analytical steps** are exposed; intermediate exploratory code is intentionally excluded.

## Data policy

- All datasets required to reproduce the **figures and tables in the paper** are provided in `data/processed/`.
- The original text of *Evgenij Onegin* is not redistributed in this repository.
- All analyses are fully reproducible from the provided derived datasets.

This design avoids unnecessary recomputation while preserving transparency and ensuring compliance 
with copyright restrictions.

## How to reproduce results

### Regenerate figures (recommended)

Figures can be regenerated from processed data by running (from the repository root):
> source("code/analysis/make_all_figures.R")

This regenerates all figures in results/figures/ using the processed datasets.

### Regenerate processed data (advanced)

High-level driver scripts in `code/analysis/` reproduce specific stages of the pipeline:

- `driver_trigrams_rus.R` — trigram extraction
- `driver_phonological_probing.R` — phonological probe analysis
- `driver_bootstrap.R` — moving block bootstrap estimation

These scripts are independent and can be run selectively.

Contextual fragments included in trigram datasets are intentionally limited to short spans 
and do not allow reconstruction of the original text.

## Dependencies

Main R packages used:

- `tidyverse`
- `tokenizers`
- `udpipe`
- `Kendall`
- `ppcor`
- `ggplot2`

Exact versions are recorded in session_info.txt.
