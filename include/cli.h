#ifndef CLI_H
#define CLI_H
#include <string>
#include <iostream>
#include <filesystem>
#include <unordered_set>
typedef struct cli_parameter{
    std::filesystem::path ref_des;
    std::filesystem::path mutant_des;
    int kmer_size=21;
    float max_ref_freq=0;
    int coverage_threshold=kmer_size; 
}cli_parameter;
bool help_parser(int argc, char* argv[]);
bool cli_parser(int argc, char* argv[],cli_parameter &parameter_set);
void print_help_line(const std::string prog);
#endif