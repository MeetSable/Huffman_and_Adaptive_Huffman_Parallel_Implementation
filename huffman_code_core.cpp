#include "huffman_code_core.h"

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cmath>
// #include <iostream>

#define MAX_TREE_HT 50



static unsigned char get_bit(std::vector<unsigned char> &bits, unsigned long i)
{
	return (bits[i / 8] >> i % 8) & 1;
}

void build_huffman_codes(Node *curr, int arr[], int top, std::unordered_map<uint8_t, huffman_code> &huffman_codes)
{
    if(curr->left){
        arr[top] = 0;
        build_huffman_codes(curr->left, arr, top+1, huffman_codes);
    }

    if(curr->right){
        arr[top] = 1;
        build_huffman_codes(curr->right, arr, top+1, huffman_codes);
    }

    if(curr->isLeaf()){
        std::vector<unsigned char> bits(std::ceil((float)top/(float)8), 0); // allocating enough memory
        unsigned long curr_byte = 0, curr_bit;
        for(int i = 0 ; i < top ; i++)
        {
            curr_byte = i / 8;
            curr_bit = i % 8;
            // change bit if one
            if(arr[i] == 1){
                bits[curr_byte] |= 1 << curr_bit;
            }
        }
        huffman_codes[curr->symbol].bits = bits;
        huffman_codes[curr->symbol].numbits = top;
        huffman_codes[curr->symbol].bits_arr_len = bits.size();
    }
}

void DeleteTree(Node *curr){
    if(curr == nullptr) return;
    DeleteTree(curr->right);
    DeleteTree(curr->left);
    delete curr;
}

pq generate_frequency_queue(std::ifstream &in_stream, unsigned long offset, unsigned long length, std::unordered_map<uint8_t, uint32_t> &symbol_frequencies)
{
    in_stream.clear(); // clearing error bits if exists
    in_stream.seekg(offset, in_stream.beg);
    // std::unordered_map<uint8_t, uint32_t> symbol_frequencies;
    char read_byte;
    for(unsigned long i = 0 ; i < length ; i++)
    {
        in_stream.read(&read_byte, 1);
        symbol_frequencies[read_byte]++;
    }
    pq frequency_queue;
    // update frequency_queue
    for(auto& item : symbol_frequencies){
        frequency_queue.push(new Node(item.first, item.second));
    }

    return frequency_queue;
}

void build_huffman_tree_from_q(pq &frequency_queue, Node** root)
{
    while(frequency_queue.size() != 1)
    {
        Node* left = frequency_queue.top();
        frequency_queue.pop();
        Node* right = frequency_queue.top();
        frequency_queue.pop();
        Node* top = new Node('$', left->freq + right->freq);
        top->left = left;
        top->right = right;
        frequency_queue.push(top);
    }
    (*root) = frequency_queue.top();
}

unsigned long total_number_of_bits(std::unordered_map<uint8_t, huffman_code> &huffman_codes, std::unordered_map<uint8_t, uint32_t> &symbol_frequencies)
{
    unsigned long total_bits = 0;
    for(auto& item : symbol_frequencies)
    {
        total_bits += huffman_codes[item.first].numbits * item.second;
    }
    return total_bits;
}

void write_huffman_codes(std::ofstream &out_stream, std::unordered_map<uint8_t, huffman_code> &huffman_codes, unsigned long total_length_bits)
{
    // write code table to file
    char write_line[512];
    int writer_cursor = 0;
    unsigned int num_codes = huffman_codes.size();
    memcpy(write_line, &total_length_bits, 8);
    out_stream.write(write_line, 8); // number of bits
    memcpy(write_line, &num_codes, 4);
    out_stream.write(write_line, 4); // number of codes in first 4 bytes

    for(auto& item : huffman_codes)
    {
        writer_cursor = 0;
        huffman_code curr_code = item.second;
        memcpy(&write_line[writer_cursor], &item.first, 1); // symbol of code
        writer_cursor++;
        memcpy(&write_line[writer_cursor], &curr_code.numbits, 4); // number of bits
        writer_cursor += 4;
        memcpy(&write_line[writer_cursor], &curr_code.bits_arr_len, 4); // length of bits array
        writer_cursor += 4;
        memcpy(&write_line[writer_cursor], curr_code.bits.data(), curr_code.bits_arr_len); // bits array
        writer_cursor += curr_code.bits_arr_len;
        out_stream.write(write_line, writer_cursor);
    }
}


void encode_append(std::ofstream &out_stream, std::ifstream &in_stream, std::unordered_map<uint8_t, huffman_code> &huffman_codes, unsigned long in_stream_offset, unsigned long in_stream_len)
{
    char write_byte = 0;
    char read_byte = 0;
    char curr_bit = 0;
    in_stream.clear();
    in_stream.seekg(in_stream_offset, in_stream.beg);
    while((in_stream_len-- > 0) && (in_stream.read(&read_byte, 1)))
    {
        huffman_code code = huffman_codes[read_byte];
        // std::cout << "Hello\n";
        for(int i = 0 ; i < code.numbits ; i++){
            write_byte |= get_bit(code.bits, i) << curr_bit;
            curr_bit++;
            // std::cout << (int)get_bit(code.bits, i);
            if(curr_bit == 8)
            {
                out_stream.write(&write_byte, 1);
                write_byte = 0;
                curr_bit = 0;
            }
        }
    }
    if(curr_bit > 0 && curr_bit < 8){
        out_stream.write(&write_byte, 1);
    }
}

void read_huffman_codes(std::ifstream &in_stream, std::unordered_map<uint8_t, huffman_code> &huffman_codes, unsigned long &total_length_bits)
{
    char read_line[512];
    int reader_cursor=0;
    unsigned int num_codes;
    in_stream.read(read_line, 8);
    total_length_bits = *((unsigned long*)read_line);
    in_stream.read(read_line, 4);
    num_codes = *((unsigned int*)read_line);
    for(int i = 0 ; i < num_codes ; i++)
    {
        huffman_code curr_code;
        char symbol;
        in_stream.read(&symbol, 1); // read symbol
        in_stream.read(read_line, 4); // read number of bits
        curr_code.numbits = *((unsigned int*)read_line);
        in_stream.read(read_line, 4); // bits array length
        curr_code.bits_arr_len = *((unsigned int*)read_line);
        in_stream.read(read_line, curr_code.bits_arr_len);
        for(int i = 0 ; i < curr_code.bits_arr_len ; i++)
            curr_code.bits.push_back(read_line[i]);
        huffman_codes[symbol].bits = curr_code.bits;
        huffman_codes[symbol].bits_arr_len = curr_code.bits_arr_len;
        huffman_codes[symbol].numbits = curr_code.numbits;
        curr_code.bits.clear();
    }
}

void build_huffman_tree_from_codes(std::unordered_map<uint8_t, huffman_code> &huffman_codes, Node** root)
{
    *root = new Node('$', 0);
    Node* curr = NULL;
    for(auto& item : huffman_codes)
    {
        curr = *root;
        for(int i = 0 ; i < item.second.numbits ; i++)
        {
            if(get_bit(item.second.bits, i) == 0)
            {
                if(curr->left == NULL) curr->left = new Node('$', 0);
                curr = curr->left;
            }
            else {
                if(curr->right == NULL) curr->right = new Node('$', 0);
                curr = curr->right;
            }
        }
        curr->symbol = item.first;
    }
}

void decode_file(std::ifstream &in_stream, std::ofstream &out_stream, Node* root, unsigned long total_length_bits)
{
    char read_byte;
    char write_byte;
    Node *curr = root;
    while((in_stream.read(&read_byte, 1)) && (total_length_bits > 0))
    {
        for(int i = 0 ; i < 8 ; i++)
        {
            if(((read_byte >> i) & 1) == 0)
                curr = curr->left;
            else
                curr = curr->right;
            if(curr->isLeaf()){
                write_byte = curr->symbol;
                out_stream.write(&write_byte, 1);
                curr = root;
            }
            total_length_bits-=1;
            if(total_length_bits == 0) break;
        }
    }
    // std::cout << total_length_bits << std::endl;
}

