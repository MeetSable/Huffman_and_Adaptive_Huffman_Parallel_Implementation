#include <cstdint>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>
#include <unordered_map>
#include <chrono>

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
    return HUFFMAN_SUCCESS;
}

uint8_t huffman_decode_parallel(std::string& in_filename, std::string& out_filename, int my_rank, int p)
{
    std::unordered_map<uint8_t, huffman_code> huffman_codes;
    std::ifstream in_stream(in_filename, std::ios::binary | std::ios::in);
    if(!in_stream.is_open()){
        std::cout << "Error opening file" << std::endl;
        return HUFFMAN_FAILURE;
    }
    unsigned long total_length_bits;
    // read huffman codes from file
    read_huffman_codes(in_stream, huffman_codes, total_length_bits);
    
    // build huffman tree
    Node* huffman_tree_root;
    build_huffman_tree_from_codes(huffman_codes, &huffman_tree_root);
    if(my_rank == 0){
        std::ofstream out_stream(out_filename, std::ios::binary | std::ios::out | std::ios::trunc);
        decode_file(in_stream, out_stream, huffman_tree_root, total_length_bits);
        out_stream.close();
        // dummy_send to trigger next write process
        char a = 0;
        MPI_Send(&a, 1, MPI_CHAR, my_rank+1, 0, MPI_COMM_WORLD);
    }
    else{
        char a;
        // wait for previous process to complete
        MPI_Recv(&a, 1, MPI_CHAR, my_rank-1, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        std::ofstream out_stream(out_filename, std::ios::binary | std::ios::in | std::ios::ate);
        // long cur = out_stream.tellp();
        // out_stream.seekp(cur-3);
        decode_file(in_stream, out_stream, huffman_tree_root, total_length_bits);
        out_stream.close();
        // send to trigger next write
        if(my_rank < p-1){
            MPI_Send(&a, 1, MPI_CHAR, my_rank+1, 0, MPI_COMM_WORLD);
        }
    }
    in_stream.close();
    return HUFFMAN_SUCCESS;
}

int main(int argc, char** argv)
{
    int my_rank, p;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    std::string in_filename(argv[1]);
    std::string encoded_filename = in_filename + "." + std::to_string(my_rank) + ".huff";
    std::string decoded_filename(argv[2]);
    unsigned long in_file_size = std::filesystem::file_size(in_filename);
    long encoding_time, decoding_time;

    std::chrono::high_resolution_clock::time_point t_start, t_end;
    MPI_Barrier(MPI_COMM_WORLD);
    t_start = std::chrono::high_resolution_clock::now();
    if (huffman_encode_parallel(in_filename, encoded_filename, my_rank, p) == HUFFMAN_FAILURE)
    {
        std::cout << "File encoding failed!!" << std::endl;
        return -1;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    t_end = std::chrono::high_resolution_clock::now();
    if(my_rank == 0)
    {
        encoding_time = std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count();
        // std::cout << "Encoding time: " << (double)encoding_time/(double)1000000 << std::endl;
    }

    unsigned long chunk_encoded_file_size = std::filesystem::file_size(encoded_filename);
    unsigned long total_encoded_size = 0;
    MPI_Reduce(&chunk_encoded_file_size, &total_encoded_size, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    t_start = std::chrono::high_resolution_clock::now();
    if (huffman_decode_parallel(encoded_filename, decoded_filename, my_rank, p) == HUFFMAN_FAILURE)
    {
        std::cout << "File decoding failed!!" << std::endl;
        return -1;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    t_end = std::chrono::high_resolution_clock::now();
    if(my_rank == 0)
    {
        decoding_time = std::chrono::duration_cast<std::chrono::microseconds>(t_end-t_start).count();
        // std::cout << "Decoding time: " << (double)decoding_time/(double)1000000 << std::endl;
        // std::cout << in_file_size << " " << total_encoded_size << " Compression ratio: " << (double)in_file_size/(double)total_encoded_size << std::endl;
        std::cout << in_file_size << " " << (double)encoding_time/(double)1000000 << " " << (double)decoding_time/(double)1000000 << " " << (double)in_file_size/(double)total_encoded_size << std::endl;
    }

    MPI_Finalize();
    return 0;
}