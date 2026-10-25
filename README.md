# Rare Sequence Finder

description

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
Rare-Seq-Finder
├── CMakeLists.txt
├── README.md
├── build/              # Stores compiled executable files
├── config/
│   └── mutation.rmt
├── data/
│   ├── demo/
│   │   ├── default_mutation_ground_truth.vcf
│   │   ├── default_mutation.fasta
│   │   └── default_ref.fasta
│   ├── MS_output/
│   ├── ref/
│   └── selected/          
├── include/            # Header files
│   ├── main.h
│   ├── cli.h
│   ├── decoder.h
│   ├── Fasta_reader.h
│   ├── Kmer_encoder.h
│   └── result_processor.h
├── src/                # Source files
│   ├── main.cpp
│   ├── cli.cpp
│   ├── decoder.cpp
│   ├── Fasta_reader.cpp
│   ├── Kmer_encoder.cpp
│   └── result_processor.cpp
├── scripts/
|   └── select_seq.py  # Select one sequence randomly from reference FASTA file
└── result/           
```
## Performance Benchmark