# Rare Sequence Finder

Rare Sequence Finder is a C++ tool for identifying candidate DNA sequence positions enriched in **rare k-mers** relative to a user-provided reference collection with alignment-free approach.

## Features

- 
- 
- 
- 
- 
- 
- 

## Quick Start

### Requirements
Linux / wsl with:

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
│   │   ├── demo_result.fasta
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