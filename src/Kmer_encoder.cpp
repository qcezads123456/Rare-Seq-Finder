#include "Kmer_encoder.h"
using namespace std;

static int parser(char nucleotide){
    switch(nucleotide){
        case 'A':
        case 'a': return 0b00;
        case 'C':
        case 'c': return 0b01;
        case 'G':
        case 'g': return 0b10;
        case 'T':
        case 't': return 0b11;
        default : return 0b111;
    }
}

/********************************************************************
                            1-based
*********************************************************************/


unordered_map <uint64_t, uint64_t> build_kmers_index   
(const vector <DNA_info_t> &genome_isolate,const int kmer_size){
    unordered_map <uint64_t, uint64_t> kmer_index;
    uint64_t kmer_mask=1ULL;
    kmer_mask=(kmer_mask<<(kmer_size*2))-1;
    for(size_t i=0;i<genome_isolate.size();i++){
        stringstream ss(genome_isolate[i].seq);
        char ch;
        int count=0;
        uint64_t encode_DNA=0;
        uint64_t isolate_mask=1ULL<<i;
        while(!ss.eof()){
            if(count!= kmer_size){
                
                while(count!=kmer_size && ss>>ch ){
                    if(parser(ch)!=0b111){
                        count++;
                        if(count>0){
                            encode_DNA<<=2;
                        }
                        encode_DNA|=parser(ch);
                    }
                    else{                   //skip nucleotide which can't encode
                        count=0;
                        encode_DNA=0;
                        break;
                    }
                }
                if(count==kmer_size){
                    kmer_index[encode_DNA] |= isolate_mask;
                }
            }
            else{
                while(ss>>ch){
                    if(parser(ch)==0b111){
                        count=0;
                        encode_DNA=0;
                        break;
                    }
                    encode_DNA <<= 2;
                    encode_DNA |= parser(ch);
                    encode_DNA &=kmer_mask;
                    kmer_index[encode_DNA] |= isolate_mask;
                }
            }
        }
    }
    return kmer_index;
}

unordered_map <uint64_t, vector <uint64_t>> find_rare_kmers   // [kmer_encode, position(actual value)]
(const DNA_info_t &genome_mutant, int kmer_size 
,const unordered_map<uint64_t,uint64_t> &kmer_index 
,const int isolate_amount,const cli_parameter &parameter_set){  //probably need struct to encapsulate all parameters
    unordered_map <uint64_t, vector<uint64_t>> low_frequency_kmer;
    uint64_t kmer_mask=1ULL;
    kmer_mask=(kmer_mask<<(kmer_size*2))-1;
    int count=0;
    uint64_t encode_DNA=0;
    stringstream ss(genome_mutant.seq);
    char ch;
    uint64_t position=0;
    int threshold = static_cast<int>(floor(isolate_amount * parameter_set.max_ref_freq));
    while(!ss.eof()){
        if(count!= kmer_size){    
            while(count!=kmer_size && ss>>ch){
                if(parser(ch)!=0b111){
                    count++;
                    position++;
                    if(count>0){
                        encode_DNA<<=2;
                    }
                    encode_DNA|=parser(ch);
                }
                else{                   //skip nucleotide which can't encode
                    count=0;
                    encode_DNA=0;
                    position++;
                    break;
                }
            }
            if(count==kmer_size){
                if(kmer_index.find(encode_DNA)!=kmer_index.end()){
                    if(count_isolates(isolate_amount,kmer_index.at(encode_DNA))<=threshold){
                        low_frequency_kmer[encode_DNA].push_back(position);
                    }
                }
                else{
                    low_frequency_kmer[encode_DNA].push_back(position);
                }
            }
        }
        else{
            while(ss>>ch){
                if(parser(ch)==0b111){
                    count=0;
                    encode_DNA=0;
                    position++;
                    break;
                }
                encode_DNA <<= 2;
                encode_DNA |= parser(ch);
                encode_DNA &=kmer_mask;
                position++;
                if(kmer_index.find(encode_DNA)!=kmer_index.end()){
                    if(count_isolates(isolate_amount,kmer_index.at(encode_DNA))<=threshold){
                        low_frequency_kmer[encode_DNA].push_back(position);
                    }
                }
                else{
                    low_frequency_kmer[encode_DNA].push_back(position);
                }
            }
        }
    }
    return low_frequency_kmer;
}