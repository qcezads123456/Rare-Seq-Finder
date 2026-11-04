#ifndef DECODER_H
#define DECODER_H

#include <string>
#include <cstdint>
#include <algorithm>
#include <vector>

uint64_t count_isolates(int isolate_amount,const std::vector<uint64_t> &key);
std::string kmer_decode(int kmer_size,uint64_t kmer);

#endif