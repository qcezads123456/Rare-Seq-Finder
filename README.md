# Rare Sequence Finder

Rare Sequence Finder is a C++17 command-line tool that uses an **alignment-free k-mer approach** to locate candidate regions in DNA sequences containing k-mers that are absent from, or uncommon across, a user-provided reference collection.

The program compares each target FASTA record against a reference k-mer index and reports positions with concentrated **rare k-mer coverage**.


## Features

- Alignment-free scanning with encoded DNA k-mers (`k = 1–31`).
- Reference-wide k-mer presence tracking using a 64-bit bitmask (up to 64 reference sequences).
- Configurable k-mer size, maximum reference frequency, and minimum rare k-mer coverage.
- Support for multi-record reference and target FASTA files.
- Per-target reporting of candidate positions and their corresponding DNA sequences.
- Optional synthetic test-data generation with [Mutation-Simulator](https://github.com/mkpython3/Mutation-Simulator), using SNP, insertion, deletion, and combined-mutation profiles.
- Makefile commands for environment setup, test-data generation, compilation, and analysis.

## Requirements

Use a Linux or WSL environment with:

- GNU Make
- CMake **3.20+**
- A C++ compiler supporting **C++17**
- **Python 3.12** with `venv` and `pip` **only if generating synthetic test data** (the Makefile uses `python3.12` by default)

The C++ analyzer can be built and used without installing Mutation-Simulator.

## Quick Start


## Project Structure
```text
Rare-Seq-Finder/
├── README.md
├── CMakeLists.txt
├── Makefile                  
├── config/
│   └── mutation.rmt             # Mutation-Simulator configuration
├── data/
│   ├── demo/
│   │   ├── demo_ref.fasta
│   │   ├── demo_mutation.fasta
│   │   ├── demo_result.txt
│   │   └── demo_mutation_ground_truth.vcf
│   ├── ref/                    # User-provided reference data
│   ├── selected/               # Sequence-selection output
│   └── MS_output/              # Mutation-simulation output directory
├── include/                    # C++ headers
├── src/                        # C++ implementation
├── scripts/
│   └── select_seq.py
├── build/                      # Compiled executable
└── result/                     # Analysis report
```
## Performance Benchmark