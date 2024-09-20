#ifndef HUFFMAN_CODE_H
#define HUFFMAN_CODE_H

#include <string>
#include <stdint.h>
#include <vector>
#include <unordered_map>
#include <queue>

#define MAX_TREE_HT 50



struct Node {
    uint8_t symbol;
    uint32_t freq;
    Node *left, *right;
    Node(const uint8_t &item, const uint32_t &freq)
        : symbol(item), freq(freq), left(NULL), right(NULL) {}
    bool isLeaf(){
        return !(this->left) && !(this->right);
    }
};

struct huffman_code {
    unsigned int numbits;
    std::vector<unsigned char> bits;
    unsigned int bits_arr_len;
};

class Comparator{
public:
    bool operator()(const Node* a, const Node* b){
        return a->freq > b->freq;
    }
};

typedef std::priority_queue<Node*, std::vector<Node*>, Comparator> pq;

void build_huffman_tree_from_q(pq &frequency_queue, Node** root);
void build_huffman_codes(Node *curr, int arr[], int top, std::unordered_map<uint8_t, huffman_code> &huffman_codes);
void write_huffman_codes(std::ofstream &out_stream, std::unordered_map<uint8_t, huffman_code> &huffman_codes, unsigned long total_length_bits);
pq generate_frequency_queue(std::ifstream &in_stream, unsigned long offset, unsigned long length, std::unordered_map<uint8_t, uint32_t> &symbol_frequencies);
unsigned long total_number_of_bits(std::unordered_map<uint8_t, huffman_code> &huffman_codes, std::unordered_map<uint8_t, uint32_t> &symbol_frequencies);
void encode_append(std::ofstream &out_stream, std::ifstream &in_stream, std::unordered_map<uint8_t, huffman_code> &huffman_codes, unsigned long in_stream_offset, unsigned long in_stream_len);
void DeleteTree(Node *curr);
void read_huffman_codes(std::ifstream &in_stream, std::unordered_map<uint8_t, huffman_code> &huffman_codes, unsigned long &total_length_bits);
void build_huffman_tree_from_codes(std::unordered_map<uint8_t, huffman_code> &huffman_codes, Node** root);
void decode_file(std::ifstream &in_stream, std::ofstream &out_stream, Node* root, unsigned long total_length_bits);

#endif