#include <iostream>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <chrono>

#include "huffman_code_core.h"

#define HUFFMAN_SUCCESS 0
#define HUFFMAN_FAILURE 1

uint8_t huffman_encode_serial(std::string &in_filename, std::string &out_filename)
{
    unsigned long long in_file_size = std::filesystem::file_size(in_filename);
    if(in_file_size == 0){
        std::cout << "Input file doesn't exist" << std::endl;
        return HUFFMAN_FAILURE;
    }

    std::ifstream in_stream(in_filename, std::ios::binary | std::ios::in);
    if(!in_stream.is_open()){
        std::cout << "Error opening in file" << std::endl;
        return HUFFMAN_FAILURE;
    }

    // compute frequency map
    std::unordered_map<uint8_t, uint32_t> symbol_frequencies;
    pq frequency_queue = generate_frequency_queue(in_stream, 0, in_file_size, symbol_frequencies);

    // build huffman tree
    Node* huffman_tree_root;
    build_huffman_tree_from_q(frequency_queue, &huffman_tree_root);

    // build huffman codes
    int arr[MAX_TREE_HT];
    std::unordered_map<uint8_t, huffman_code> huffman_codes;
    build_huffman_codes(huffman_tree_root, arr, 0, huffman_codes);

    DeleteTree(huffman_tree_root);

    // encode and write back to a new file 
    std::ofstream out_stream(out_filename, std::ios::binary | std::ios::out | std::ios::trunc);

    unsigned long total_length_bits = total_number_of_bits(huffman_codes, symbol_frequencies);
    write_huffman_codes(out_stream, huffman_codes, total_length_bits);
    
    encode_append(out_stream, in_stream, huffman_codes, 0, in_file_size);

    in_stream.close();
    out_stream.close();
    return HUFFMAN_SUCCESS;
}


uint8_t huffman_decode_serial(std::string& in_filename, std::string& out_filename)
{
    // read huffman code table from file
    std::unordered_map<uint8_t, huffman_code> huffman_codes;
    std::ifstream in_stream(in_filename, std::ios::binary | std::ios::in);
    if(!in_stream.is_open()){
        std::cout << "Error opening file" << std::endl;
        return HUFFMAN_FAILURE;
    }

    // read huffman codes from file
    unsigned long total_length_bits;
    read_huffman_codes(in_stream, huffman_codes, total_length_bits);

    // build huffman tree
    Node* huffman_tree_root;
    build_huffman_tree_from_codes(huffman_codes, &huffman_tree_root);

    std::ofstream out_stream(out_filename, std::ios::binary | std::ios::out | std::ios::trunc);
    decode_file(in_stream, out_stream, huffman_tree_root, total_length_bits);

    in_stream.close();
    out_stream.close();

    return HUFFMAN_SUCCESS;
}

int main(int argv, char** argc){
    std::string in_filename(argc[1]);
    std::string encoded_filename = in_filename + ".huff";
    std::string decoded_filename(argc[2]);
    
    std::chrono::high_resolution_clock::time_point t_start, t_end;
    t_start = std::chrono::high_resolution_clock::now();
    if(huffman_encode_serial(in_filename, encoded_filename) == HUFFMAN_FAILURE)
    {
        std::cout << "File encoding failed!!" << std::endl;
        return -1;
    }
    t_end = std::chrono::high_resolution_clock::now();
    long encoding_time = std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count();
    std::cout << "Encoding time: " << (double)encoding_time/(double)1000000 << std::endl;

    // std::cout << "File encoding successful!!" << std::endl;
    t_start = std::chrono::high_resolution_clock::now();
    if(huffman_decode_serial(encoded_filename, decoded_filename) == HUFFMAN_FAILURE)
    {
        std::cout << "failed to decoded!!" << std::endl;
        return -1;
    }
    t_end = std::chrono::high_resolution_clock::now();
    long decoding_time = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count();
    std::cout << "Decoding time: " << (double)decoding_time/(double)1000000 << std::endl;
    
    // std::cout << "File decoded!!" << std::endl;
    unsigned long in_file_size = std::filesystem::file_size(in_filename);
    unsigned long out_file_size = std::filesystem::file_size(encoded_filename);
    unsigned long deocded_file_size = std::filesystem::file_size(decoded_filename);

    std::cout << in_file_size << " " << out_file_size << " Compression ratio: " << (double)in_file_size/(double)out_file_size << std::endl;


    return 0;
}