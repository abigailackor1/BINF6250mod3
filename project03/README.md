# Introduction
This project implements a Gibbs sampler to identify short, shared DNA motifs across a set of sequences — a common problem in molecular biology where a regulatory signal (e.g., a transcription factor binding site or ribosome binding site) is known to exist somewhere in each sequence, but its exact position and composition aren't known in advance.

Because real genomic datasets are far too large to exhaustively check every possible motif position across every sequence, we use Gibbs Sampling - a Markov Chain Monte Carlo (MCMC) approach.

The algorithm starts from a random guess at the motif's location in each sequence, then iteratively refines those guesses: on each iteration, one sequence is set aside, a position weight matrix (PWM) is built from the current guesses in every other sequence, and the set-aside sequence's guess is updated by scoring all possible windows against that PWM and sampling a new position probabilistically (rather than always taking the best-scoring window). Repeating this thousands of times allows the guesses to converge on the sequences' shared motif, without ever exhaustively searching the full solution space.

We test our implementation on two datasets: 
1. GCF_000009045.1_ASM904v1_genomic.fna and GCF_000009045.1_ASM904v1_genomic.gff\
Promoter regions upstream of Bacillus subtilis coding sequences, pre-filtered for a fragment of the Shine-Dalgarno motif
2. nrf1_gibbs.fa\
NRF1 ChIP-seq peak sequences, to recover the motif associated with NRF1 transcription factor binding.

# Usage
Open `project03.ipynb` in Jupyter and run all cells in order. The notebook is split into:
- Core deliverable: **Implement Gibbs Sampler** — the `GibbsMotifFinder()` function
- Provided Programs, not to be modified:
**Driver Program** - runs `GibbsMotifFinder` on the *B. subtilis* promoter dataset and plots the resulting sequence logo.\
**NRF1 Driver Program** — runs `GibbsMotifFinder` on the NRF1 ChIP-seq peaks. Requires completing the data-ingest cell above it.
Input data files are NOT tracked in this repo, they are listed in Project Structure below and placed in local folder before running the notebook.

# Pseudocode for project03.ipynb
### GibbsMotifFinder
```
GibbsMotifFinder(seqs, k, seed, ic_window, max_iterations)
    0. SET UP
       seed the random generator
       stop with an error if ic_window or max_iterations < 1
       uppercase all seqs; drop any shorter than k
       N ← number of seqs
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
          print IC every 1000 iterations (sanity check)
          if IC unchanged over the last ic_window rounds:
              stop early
    3. RETURN the PFM of all Motifs (4 × k)
```
**Note #1: `random` module instead of rng object:**
The notebook's "Important considerations" list `random.randint()/numpy.random.randint()` and `random.choices()/ numpy.random.choice()` as equally valid options. We disregard the provided: `rng = np.random.default_rng(seed)` since our implementation uses `random.randint()/random.choices()` throughout, making the rng object unused.

Note #2: parameters instead of hardcoded values:
- `max_iterations` externalizes the pseudocode's "1 to 10000" cap, allowing
  shorter test runs (e.g., for debugging on very large datasets, like NRF1) without
  modifying the function itself.
- `ic_window` externalizes the convergence check ("or Motifs stops
  changing"), which the pseudocode leaves undefined numerically. We used a window of 100 iterations; as a parameter rather than hardcoding it, since the right convergence window may differ depending on dataset size and signal strength.
## Dependencies
- Python 3.14.2
- numpy
- [bamnostic]
- [seqlogo] 

## Project Structure
```
project03/
├── project03.ipynb        # main notebook — implements GibbsMotifFinder
├── data_readers.py          # FASTA/GFF file readers (do not modify)
├── seq_ops.py                # reverse complement + promoter extraction (do not modify)
├── motif_ops.py               # PFM/PWM building, scoring, information content (do not modify)
└── README.md
Note: Data files are not in the repo but stored locally

```
## Helper Modules

These modules are provided by the course and to be remain unmodified — `GibbsMotifFinder` is built against their function signatures.

**`data_readers.py`**
- `get_fasta(file)` — generator that lazily yields `(name, sequence)` tuples from a FASTA file (handles `.gz` compression).
- `get_gff(gff_file)` — generator that lazily yields `GffEntry` objects (seqid, type, start, end, strand) from a GFF file (handles `.gz compression).

**`seq_ops.py`**
- `reverse_complement(seq)` — returns the reverse complement of a DNA sequence string.
- `get_seq(seq, start, end, strand, size)` — extracts the promoter region upstream of a gene, correctly handling both `+` and `-` strand genes.

**`motif_ops.py`**
- `build_pfm(sequences, length)` — builds a Position Frequency Matrix (4 × length) from a list of equal-length motif strings.
- `build_pwm(pfm)` — converts a PFM into a Position Weight Matrix using log-odds scoring against a uniform (0.25) background.
- `score_kmer(seq, pwm)` — scores a single k-mer against a PWM.
- `pfm_ic(pfm)` — computes the information content (in bits) of a PFM; used to monitor convergence of the Gibbs sampler.
  
# Successes
- Successfully implemented the full Gibbs sampling loop (leave-one-out PWM construction, both-strand scoring, and probabilistic weighted selection) matching the pseudocode's formula.
- Ran `GibbsMotifFinder` on the real *B. subtilis* promoter dataset (800+ sequences, each 50bp) and had it complete successfully in approximately 1 minute, returning a valid 4×k PFM.
- Implemented a convergence check using a sliding window of information content (IC) values, exiting early once IC stabilizes across the window rather than always running the full 10,000-iteration ceiling — made the window size (`ic_window`) and iteration cap (`max_iterations`) configurable parameters rather than hardcoded values, for reusability.
- Diagnosed and resolved a data file naming mismatch between the driver program's hardcoded `.gz` paths and the actual (uncompressed) provided files.
- Diagnosed a Ghostscript dependency failure in the `seqlogo` plotting step and confirmed (per the assignment's own note) that this is a known, non-blocking issue separate from the core algorithm's correctness.

## Results Interpretation
- Promoter set:
`GibbsMotifFinder` converged to a total information content of 12.12 bits out of a theoretical maximum of 20 bits for a 10-position motif — roughly 60% of maximum, indicating that conservation is concentrated in a subset of positions rather than spread evenly across the full window. Per-position analysis confirms this directly: positions 4–9 show 98–99% consensus on a single base, spelling AGGAGG, while positions 1–3 and 10 remain weak and near-random (41–54%). This matches expectation, since the driver program pre-filtered promoters to only include sequences already containing this exact string (a fragment of the Shine-Dalgarno motif) — meaning this result primarily validates that the algorithm correctly recovers a known, real signal from random initialization, distinguishing it clearly from uninformative flanking sequence, rather than demonstrating discovery of a previously unknown motif.
- NRF1 ChIP-seq peak sequences: still running, to be updated 
# Struggles
TRang (im just putting my name here bc this was my struggle, please add yours and considalte them with mine, remove the name. Format: issue => what we found out => what was our debug action: 
- `.gz` path mismatch in driver cell\
→ `get_fasta()`/`get_gff()` only gzip-open when `.gz` is in the filename string; provided files were uncompressed, causing `FileNotFoundError`.\
→ Updated the two hardcoded path strings to match actual filenames.
- Ghostscript missing (`OSError`)\
→ `seqlogo` depends on `weblogo`, which requires the external Ghostscript program on PATH — not something `pip`/`conda` installs.\
→ Installed via `brew install ghostscript`, restarted kernel.
- Duplicate nested `project03/project03/` from a teammate's PR\
→ Caused by extracting the project zip inside an already-existing folder of the same name.\
→ Moved real edits to the correct path, deleted the duplicate.
- NRF1 runtime (90,061 sequences, ~1 sec/iteration)\
→ `build_pfm()` reprocesses nearly the full sequence list every iteration, so runtime scales directly with dataset size— a full 10,000-iteration run would take roughly 2.8 hours. Measured real per-iteration cost with a short test run(`max_iterations=20`) to get grounded timing instead of guessing; let the full run continue in the background given the Oct 7 deadline.
- First-time Terminal clone/push\
→ `git clone` targets the whole repo, not a branch; switching branches is always a separate step.\
→ Practiced the full clone → branch → pull → edit → commit → push cycle.\
# Personal Reflections
## Group Leader
Trang: Our approach was making sure every teammate understand and write the Gibbs sampler independently and cross-validating results against each other before merging. This allows us to struggle and learn, and I do think the efforts and time spent was worth it because each teammate own implementation works correctly without needing to follow the others, and not only that, we can explain what the output means and confidently say that we understand the project.
Project 3 is my first time being a group leader and owning repo. Navigating and untangling confusing as a repo owner on here was time-consuming to me, but then again, practice makes perfect and I really appreciate the chances to do this more often. 
## Other members
Dianah: Project 3 is my first time being a collaborator instead of a project leader, so I am learning that side of GitHub as I go, trying to figure out forking, 
how to open a pull request into someone else’s branch instead of my own and what my responsibilities look like when I am not managing the whole repo. 
I am still getting my footing with it, and I think it’s helping me understand GitHub better.

# Generative AI Appendix
**Tool used:** Claude (Anthropic), Claude Sonnet 5
**Prompt:** 
**Used for:**
