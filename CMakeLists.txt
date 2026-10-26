cmake_minimum_required(VERSION 3.20)
project(Rare-Seq-Finder LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(prog rare-seq-finder)
add_executable(${prog}
    src/cli.cpp
    src/decoder.cpp
    src/Fasta_reader.cpp
    src/Kmer_encoder.cpp
    src/main.cpp
    src/result_processor.cpp
)
target_include_directories(${prog}
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include    
)
set_target_properties(${prog}
    PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)