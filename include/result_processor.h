#ifndef RESULT_PROCESSOR_H
#define RESULT_PROCESSOR_H
#include <iostream>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include "Kmer_encoder.h"
#include "cli.h"
#include "Fasta_reader.h"
void result_process
(const std::vector <DNA_info_t> &genome_isolate 
,const std::vector <DNA_info_t> &genome_mutant
,const cli_parameter &parameter_set);
#endif