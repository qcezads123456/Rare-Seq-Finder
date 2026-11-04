#ifndef KMER_ENCODER_H
#define KMER_ENCODER_H

#include "Fasta_reader.h"
#include "decoder.h"
#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <cstdint>
#include <cmath>
#include "cli.h"
std::unordered_map <uint64_t, std::vector <uint64_t>> build_kmers_index(const std::vector <DNA_info_t> &genome_isolate,const int kmer_size);
std::unordered_map <uint64_t, std::vector <uint64_t>> find_rare_kmers  
(const DNA_info_t &genome_mutant, int kmer_size,const std::unordered_map<uint64_t,std::vector<uint64_t>> &kmer_count,int isolate_amount,const cli_parameter &parameter_set);
#endif
