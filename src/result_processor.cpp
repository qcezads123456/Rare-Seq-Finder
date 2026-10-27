#include "result_processor.h"
using namespace std;

static bool result_file_create(fstream &file){
    file.open("./result/result.txt",ios::out);
    if(!file.is_open()){
        return 0;
    }
    else return 1;
}
static bool result_tsv_file_create(fstream &file){
    file.open("./result/result.tsv",ios::out);
    if(!file.is_open()){
        return 0;
    }
    else return 1;
}
void result_process
(const vector <DNA_info_t> ref_seq 
,const vector <DNA_info_t> mutant_seq
,const cli_parameter parameter_set){
    unordered_map <uint64_t, uint64_t> ref_kmer_count;
    unordered_map <uint64_t, vector<uint64_t>> rare_kmer_position;
    fstream file,file_tsv;
    if(!result_file_create(file)){
        throw runtime_error("file create failure\n");
    }
    if(!result_tsv_file_create(file_tsv)){
        throw runtime_error("file create failure\n");
    }
    ref_kmer_count=build_kmers_index(ref_seq,parameter_set.kmer_size);
    file<<"Analysis Report\n\n"
    <<"Reference sequences: "<<ref_seq.size()<<endl
    <<"k-mer size: "<<parameter_set.kmer_size
    <<"\ncoverage threshold: "<<parameter_set.coverage_threshold
    <<"\nmax reference frequency: "<<parameter_set.max_ref_freq*100<<"%\n\n"
    <<"------------------------------------------------------------------------\n\n";
    file_tsv<<"candidate_id\tseq_id\tstart\tend\tsequence\tpeak_coverage\n";
    for(int i=0;i<mutant_seq.size();i++){
        rare_kmer_position=find_rare_kmers(mutant_seq[i],parameter_set.kmer_size,ref_kmer_count,ref_seq.size(),parameter_set);
        const string& seq = mutant_seq[i].seq;
        int size = seq.length();
        vector<int> frequent(size, 0);
        file<<mutant_seq[i].header<<"\n\n";
        const string ID=mutant_seq[i].header.substr(1,mutant_seq[i].header.find(' ')-1);
        for (const auto& [kmer, position] : rare_kmer_position) {
            for(auto & pos : position){
                int start = pos - parameter_set.kmer_size;
                int end = pos;
                if (start < 0 || end > size)
                continue;

                for (int i = start; i < end; i++) {
                    frequent[i]++;
                }
            }
        }
        int count=1;
        for(int start=0; start<size;){
            int end=start;
            while(end+1 < size && frequent[end+1]==frequent[end]){
                end++;
            }
            int coverage=frequent[end];
            if(frequent[start]>=parameter_set.coverage_threshold && 
                (start==0|| frequent[start-1]< frequent[start])&&(end==size-1 || frequent[end+1]< frequent[end])){
                    file<<start+1<<" - "<< end+1<<"  ";
                    file<<seq.substr(start,end-start+1)<<endl;
                    file_tsv<<i+1<<"-"
                            <<count++<<"\t"
                            <<ID<<"\t"
                            <<start+1<<"\t"
                            <<end+1<<"\t"
                            <<seq.substr(start,end-start+1)<<"\t"
                            <<coverage<<"\n";
                            
                }
            start=end+1;
        }
        file<<"\n\n";
    }
}