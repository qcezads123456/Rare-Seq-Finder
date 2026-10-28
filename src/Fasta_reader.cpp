#include "Fasta_reader.h"
using namespace std;
namespace fs = std::filesystem;

vector <DNA_info_t> Fasta_reader(const fs::path entry){
    ifstream file(entry);
    if (!file.is_open()) {
        cerr << "Cannot open file!" << endl;
        return {};
    } 
    vector <DNA_info_t> genome_isolate;
    string line;
    while (getline(file, line)) {
        if (!line.empty() && line[0] == '>') {
            genome_isolate.emplace_back();
            genome_isolate.back().header = line.substr(1);
        }
        else if (!line.empty()) {
            genome_isolate.back().seq += line;
        }
    }
    return genome_isolate;
}