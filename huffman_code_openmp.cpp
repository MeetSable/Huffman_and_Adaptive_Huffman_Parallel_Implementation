#include <cstdint>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <cmath>

#include "huffman_code_core.h"

#define HUFFMAN_SUCCESS 0
#define HUFFMAN_FAILURE 1

uint8_t huffman_encode_parallel(std::string& in_filename, std::string& out_filename, int my_rank, int p)
{
    unsigned long long in_file_size = std::filesystem::file_size(in_filename);
    if(in_file_size == 0) {
        std::cout << "Input file doesn't exist" << std::endl;
        return HUFFMAN_FAILURE;
    }

    std::ifstream in_stream(in_filename, std::ios::binary | std::ios::in);
    if(!in_stream.is_open()){
        std::cout << "Error opening encoded file" << std::endl;
        return HUFFMAN_FAILURE;
    }
    unsigned long chunk_length, chunk_offset;
    if(in_file_size%p==0){
        chunk_length = (double)in_file_size/(double)p;
        chunk_offset = chunk_length * my_rank;
    }
    else {
        chunk_length = std::floor((double)in_file_size/(double)(p-1));
        chunk_offset = chunk_length * my_rank;
        if(my_rank == p-1) 
        {
            chunk_length = in_file_size % chunk_length == 0 ? chunk_length : in_file_size % chunk_length; // if the file size isn't a multiple of p
        }
    }
    // std::cout << my_rank << " " << chunk_length << " " << chunk_offset << std::endl;

    std::unordered_map<uint8_t, uint32_t> symbol_frequencies;
    pq frequency_queue = generate_frequency_queue(in_stream, chunk_offset, chunk_length, symbol_frequencies);

    Node* huffman_tree_root;
    build_huffman_tree_from_q(frequency_queue, &huffman_tree_root);

    int arr[MAX_TREE_HT];
    std::unordered_map<uint8_t, huffman_code> huffman_codes;
    build_huffman_codes(huffman_tree_root, arr, 0, huffman_codes);

    DeleteTree(huffman_tree_root);

    unsigned long total_length_bits = total_number_of_bits(huffman_codes, symbol_frequencies);
    std::ofstream out_stream(out_filename, std::ios::binary | std::ios::out | std::ios::trunc);
    write_huffman_codes(out_stream, huffman_codes, total_length_bits);
    encode_append(out_stream, in_stream, huffman_codes, chunk_offset, chunk_length);

    in_stream.close();
    out_stream.close();
    // std::cout << "Encoding end: " << my_rank << std::endl;
    return HUFFMAN_SUCCESS;
}

uint8_t huffman_decode_parallel(std::vector<std::string> &in_filenames, std::string& out_filename)
{
    // std::cout << "decode start: 0" << std::endl;
    std::unordered_map<uint8_t, huffman_code> huffman_codes;
    std::ifstream in_stream(in_filenames[0], std::ios::binary | std::ios::in);
    unsigned long total_length_bits;
    read_huffman_codes(in_stream, huffman_codes, total_length_bits);

    // build huffman tree
    Node* huffman_tree_root;
    build_huffman_tree_from_codes(huffman_codes, &huffman_tree_root);

    std::ofstream out_stream(out_filename, std::ios::binary | std::ios::out | std::ios::trunc);
    decode_file(in_stream, out_stream, huffman_tree_root, total_length_bits);

    huffman_codes.clear();
    DeleteTree(huffman_tree_root);
    out_stream.close();
    in_stream.close();
    // std::cout << "decode end: 0" << std::endl;

    #pragma omp parallel for ordered
    for(int i = 1 ; i < in_filenames.size() ; i++)
    {
        // std::cout << "decode start: " << i << std::endl;
        std::unordered_map<uint8_t, huffman_code> huffman_codes;
        std::ifstream in_stream(in_filenames[i], std::ios::binary | std::ios::in);
        read_huffman_codes(in_stream, huffman_codes, total_length_bits);

        // build huffman tree
        Node* huffman_tree_root;
        build_huffman_tree_from_codes(huffman_codes, &huffman_tree_root);

        #pragma omp ordered
        {
            std::ofstream out_stream(out_filename, std::ios::binary | std::ios::out | std::ios::app);
            decode_file(in_stream, out_stream, huffman_tree_root, total_length_bits);
            huffman_codes.clear();
            DeleteTree(huffman_tree_root);
            out_stream.close();
            in_stream.close();
            // std::cout << "decode end: " << i << std::endl;
        }
    }

    return HUFFMAN_SUCCESS;
}

int main(int argc, char** argv)
{
    std::string in_filename(argv[1]);
    std::string decoded_filename(argv[2]);
    std::vector<std::string> encoded_filenames;

    unsigned long in_file_size = std::filesystem::file_size(in_filename);

    std::chrono::high_resolution_clock::time_point t_start, t_end;
    t_start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel shared(encoded_filenames)
    {
        int nthreads = omp_get_num_threads();
        int tid = omp_get_thread_num();

        #pragma omp master
        {
            encoded_filenames.resize(nthreads);
        }
        #pragma omp barrier

        encoded_filenames[tid] = in_filename + "." + std::to_string(tid) + ".huff";

        huffman_encode_parallel(in_filename, encoded_filenames[tid], tid, nthreads);
    }
    t_end = std::chrono::high_resolution_clock::now();
    long encoding_time = std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count();
    // std::cout << "Encoding time: " << (double)encoding_time/(double)1000000 << std::endl;

    unsigned long total_encoded_file_size = 0;
    for(auto& encoded_filename : encoded_filenames)
        total_encoded_file_size += std::filesystem::file_size(encoded_filename);

    t_start = std::chrono::high_resolution_clock::now();
    huffman_decode_parallel(encoded_filenames, decoded_filename);
    t_end = std::chrono::high_resolution_clock::now();

    long decoding_time = std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count();
    // std::cout << "Decoding time: " << (double)decoding_time/(double)1000000 << std::endl;
    // std::cout << in_file_size << " " << total_encoded_file_size << " Compression ratio: " << (double)in_file_size/(double)total_encoded_file_size << std::endl;
    std::cout << in_file_size << ", " << (double)encoding_time/(double)1000000 << ", " << (double)decoding_time/(double)1000000 << ", " << (double)in_file_size/(double)total_encoded_file_size << std::endl;




    return 0;
}