#include "main.h"

using namespace std;
namespace fs=std::filesystem;
int main(int argc,char* argv[]){
    vector <DNA_info_t> genome_isolate,genome_mutant;
    cli_parameter parameter_set;
    const char* env=getenv("RUN_VIA_MAKE");
    const string prog= env ? env : fs::path(argv[0]).filename().string();
    if(help_parser(argc,argv)){
        print_help_line(prog);
        return 0;
    }
    if (!cli_parser(argc, argv, parameter_set)) {
        cerr<<"ERROR: Invalid arguments.\n\n";
        print_help_line(prog);

        return EXIT_FAILURE;
    }
    genome_isolate=Fasta_reader(parameter_set.ref_des);
    if(genome_isolate.empty()) return EXIT_FAILURE;
    genome_mutant=Fasta_reader(parameter_set.mutant_des);
    if(genome_mutant.empty()) return EXIT_FAILURE;
    try{
        result_process(genome_isolate, genome_mutant,parameter_set);
    }
    catch (const exception & e){
        cerr<<"result process failure";
        return EXIT_FAILURE;
    }
    cout<<"analysis success\n";
    return 0;
}