# Introduction
This project implements a Gibbs sampler to identify short, shared DNA motifs across a set of sequences — a common problem in molecular biology where a regulatory signal (e.g., a transcription factor binding site or ribosome binding site) is known to exist somewhere in each sequence, but its exact position and composition aren't known in advance.

Because real genomic datasets are far too large to exhaustively check every possible motif position across every sequence, we use Gibbs Sampling - a Markov Chain Monte Carlo (MCMC) approach.

The algorithm starts from a random guess at the motif's location in each sequence, then iteratively refines those guesses: on each iteration, one sequence is set aside, a position weight matrix (PWM) is built from the current guesses in every other sequence, and the set-aside sequence's guess is updated by scoring all possible windows against that PWM and sampling a new position probabilistically (rather than always taking the best-scoring window). Repeating this thousands of times allows the guesses to converge on the sequences' shared motif, without ever exhaustively searching the full solution space.

We test our implementation on two datasets: 
(1) GCF_000009045.1_ASM904v1_genomic.fna and GCF_000009045.1_ASM904v1_genomic.gff
Promoter regions upstream of Bacillus subtilis coding sequences, pre-filtered for a fragment of the Shine-Dalgarno motif
(2) nrf1_gibbs.fa
NRF1 ChIP-seq peak sequences, to recover the motif associated with NRF1 transcription factor binding.

# Usage
Open `project03.ipynb` in Jupyter and run all cells in order. The notebook is split into:
- Core deliverable: **Implement Gibbs Sampler** — the `GibbsMotifFinder()` function
- Provided Programs, not to be modified:
**Driver Program** - runs `GibbsMotifFinder` on the *B. subtilis* promoter dataset and plots the resulting sequence logo. 
**NRF1 Driver Program** — runs `GibbsMotifFinder` on the NRF1 ChIP-seq peaks. Requires completing the data-ingest cell above it.
Input data files are NOT tracked in this repo, they are listed in Project Structure below and placed in local folder before running the notebook.

# Pseudocode for project03.ipynb
GibbsMotifFinder(seqs, k, seed, ic_window, max_iterations)
    0. SET UP
       seed the random generator
       stop with an error if ic_window or max_iterations < 1
       uppercase all seqs; drop any shorter than k
       N ← number of seqs
Note: The notebook's "Important considerations" explicitly list random.randint()/numpy.random.randint() and random.choices()/ numpy.random.choice() as equally valid options. We chose the `random` module for both steps, hence we disregard the provided: rng = np.random.default_rng(seed) since our implementation uses random.randint()/random.choices() throughout, making the rng object unused.
    1. INITIALIZE
       for each seq:
           motif ← random k-letter piece of seq
       Motifs ← list of all these motifs
    2. ITERATE (up to max_iterations times)
       a. i ← random sequence number (0 to N-1)
       b. PWM ← built from all Motifs except Motif_i
       c. for every k-letter window in seq_i (skip windows with N):
              score the window and its reverse complement against PWM
       d. A ← 2^score for each candidate
          P ← A / sum(A)
          Motif_i ← one candidate, picked at random weighted by P
       e. IC ← information content of all Motifs; record it
          if IC unchanged over the last ic_window rounds:
              stop early
    3. RETURN the PFM of all Motifs (4 × k)
## Dependencies
- Python 3.14.2
- numpy
- [bamnostic]
- [seqlogo] (required for the plotting cells, not for the core algorithm)
  
## Project Structure
\`\`\`
project03/
├── project03.ipynb        # main notebook — implements GibbsMotifFinder
├── data_readers.py          # FASTA/GFF file readers (do not modify)
├── seq_ops.py                # reverse complement + promoter extraction (do not modify)
├── motif_ops.py               # PFM/PWM building, scoring, information content (do not modify)
└── README.md
\`\`\`



```


```

# Successes

# Struggles
Trang's note: algo works I think, will need to do sth about Ghostscript for visualization, will mention that here. 
# Personal Reflections
## Group Leader

## Other members
Dianah: Project 3 is my first time being a collaborator instead of a project leader, so I am learning that side of GitHub as I go, trying to figure out forking, 
how to open a pull request into someone else’s branch instead of my own and what my responsibilities look like when I am not managing the whole repo. 
I am still getting my footing with it, and I think it’s helping me understand GitHub better.

# Generative AI Appendix
