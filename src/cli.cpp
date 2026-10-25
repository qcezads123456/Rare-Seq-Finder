#include "cli.h"
using namespace std;
namespace fs= std::filesystem;
unordered_set<string> extend={".fa", ".fasta", ".fna", ".fas"};
bool help_parser(int argc, char* argv[]){
    for(int i=0;i<argc;i++){
        string arg=argv[i];
        if(help_CLI()) return 1;
    }
    return 0;
}
bool cli_parser(int argc,char* argv[],cli_parameter &parameter_set){
    bool ref_set=0;
    bool mutant_set=0;
    bool ct_set=0;
    for(int i=1;i<argc;i++){
        string arg=argv[i];
        if (i + 1 >= argc) return false;
        if(max_ref_freq_condition()){
            float num=stof(argv[++i]);
            if(num<0 || num>=1.0) return 0;
            parameter_set.max_ref_freq=num;
        }
        else if(kmer_size_conditoin()){
            parameter_set.kmer_size=stoi(argv[++i]);
            if(parameter_set.kmer_size>=32 || parameter_set.kmer_size<=0) return 0;
        }
        else if(ref_destination()){
            fs::path entry=argv[++i];
            if(compare()) return 0;
            parameter_set.ref_des=entry;
            ref_set=1;
        }
        else if(tested_destination()){
            fs::path entry=argv[++i];
            if(compare()) return 0;
            parameter_set.mutant_des=entry;
            mutant_set=1;
        }
        else if(coverage_threshold_condition()){
            parameter_set.coverage_threshold=stoi(argv[++i]);
            ct_set=1;
        }
        else return 0;
    }
    if(!ct_set) parameter_set.coverage_threshold=parameter_set.kmer_size;
    if(parameter_set.coverage_threshold<=0 || parameter_set.coverage_threshold >parameter_set.kmer_size) return 0;
    return ref_set & mutant_set;
}