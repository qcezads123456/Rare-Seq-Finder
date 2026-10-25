#ifndef CLI_H
#define CLI_H
#include <string>
#include <iostream>
#include <filesystem>
#include <unordered_set>
#define compare() extend.find(entry.extension()) == extend.end()
#define max_ref_freq_condition() (arg=="-mrf" || arg== "--max-ref-freq")
#define coverage_threshold_condition() (arg=="-ct" || arg== "--coverage")
#define kmer_size_conditoin() (arg=="-k" || arg== "--kmer-size")
#define ref_destination() (arg=="-r" || arg== "--reference")
#define tested_destination() (arg=="-m" || arg== "--mutant")
#define help_CLI() (arg=="-h" || arg== "--help")
typedef struct cli_parameter{
    std::filesystem::path ref_des;
    std::filesystem::path mutant_des;
    int kmer_size=21;
    float max_ref_freq=0;
    int coverage_threshold=kmer_size; 
}cli_parameter;
bool help_parser(int argc, char* argv[]);
bool cli_parser(int argc, char* argv[],cli_parameter &parameter_set);
#endif