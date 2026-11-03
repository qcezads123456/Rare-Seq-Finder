# Rare Sequence Finder

Rare Sequence Finder is a C++17 command-line tool that uses an **alignment-free k-mer approach** to locate candidate regions in DNA sequences containing k-mers that are absent from, or uncommon among the reference collection.

**Only reference and target FASTA files are needed for detection.** Simulation is optional and is used for validation.


## Features

- Encoded DNA k-mer analysis with configurable `k` (`1–31`).
- Reference k-mer tracking across up to **64 reference sequences**.
- Adjustable reference-frequency and rare k-mer coverage thresholds.
- Multi-record reference and target FASTA support.
- Candidate coordinates and sequences in TXT/TSV reports, with JSON metadata.
- Optional simulation-based validation for SNPs, insertions, deletions, and mixed mutations.


## Requirements

**Core detector (Linux or WSL):** GNU Make, CMake 3.20+, and a C++17-compatible compiler.

**Optional validation:** Python 3.12 with `venv` and `pip`. Mutation-Simulator and the analysis dependencies are installed through the Makefile. Python is not needed for normal detection.


## Quick Start

## Validation Result

## Project Structure

```text
bioinfo_project/
├── Makefile                  # Build, run, and validation commands
├── CMakeLists.txt
├── src/                      # C++ implementation
├── include/                  # C++ headers
├── data/
│   ├── demo/                 # Built-in demo
│   ├── ref/                  # Reference input for default workflow
│   ├── selected/             # Selected simulation input
│   ├── MS_output/            # Simulated FASTA and VCF
│   └──validation_results/    # Validation figures and data
├── config/                   # Simulation configuration
├── scripts/                  # Simulation input preparation
├── analysis/                 # Validation summary
├── build/                    # Build output
└── result/                   # Detector output
```

