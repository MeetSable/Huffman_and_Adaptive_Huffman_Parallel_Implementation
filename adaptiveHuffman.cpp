#include <iostream>
#include <map>
#include <sstream>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <bitset>
#include <time.h>
#include <vector>
#include <omp.h>
#include<bits/stdc++.h>

class AdaptiveHuffmanCoding
{
public:
	AdaptiveHuffmanCoding();
	void UpdateTreeModel(int);
	std::string Encode(int);
	std::string Encode2(int);
	std::string Decode(const std::string&);
	int Decode2(std::istream&);

	
	~AdaptiveHuffmanCoding();
	//Used in order to specify the end of the compressed text => on live decompressing the decoder will know when to stop
	static const int PSEUDO_EOF = 256;
private:
	AdaptiveHuffmanCoding(const AdaptiveHuffmanCoding&); //both copy constructor and assignment operator are useless in this case
	AdaptiveHuffmanCoding& operator= (const AdaptiveHuffmanCoding&);
	int e = 4;
	int r = 10;
	struct HuffmanNode
	{
		int value;      //8-bit character contained in tree node
		int weight;   //number of times val has occured in file so far
		int order;    //ordering system to track weights

		HuffmanNode* left;
		HuffmanNode* right;
		HuffmanNode* parent;

		bool isNYT;

		HuffmanNode(int value, int weight, int order, HuffmanNode* left, HuffmanNode* right,
			HuffmanNode* parent, bool isNYT = false)
			:value(value), weight(weight), order(order), left(left), right(right), parent(parent), isNYT(isNYT)
		{}
		HuffmanNode()
			: value(0), weight(0), order(0), left(nullptr), right(nullptr), parent(nullptr), isNYT(false)
		{}

		HuffmanNode(bool isNYT, int value)
			: HuffmanNode()
		{
			this->isNYT = isNYT;
			this->value = value;
		}
	};

	HuffmanNode* GetSymbolNode(int, HuffmanNode*) const;
	void FindBlockLeader(HuffmanNode*, HuffmanNode*&) const;
	void SwapNodes(HuffmanNode*, HuffmanNode*);
	std::string GetPathToSymbol(HuffmanNode*, HuffmanNode*, std::string) const;
	void SlideAndIncrement(HuffmanNode*);

	void DeleteHuffmanTree(HuffmanNode*);

	HuffmanNode	*NYTNode;
	HuffmanNode *root;
};

AdaptiveHuffmanCoding::AdaptiveHuffmanCoding()
	: root(new HuffmanNode(-1, 0, 512, nullptr, nullptr, nullptr, true))
{
	NYTNode = root;
}

void AdaptiveHuffmanCoding::UpdateTreeModel(int symbol)
{
	HuffmanNode *leafToIncrement = nullptr;
	HuffmanNode *symbolNode = GetSymbolNode(symbol, root);

	if (symbolNode == nullptr)
	{
		NYTNode->isNYT = false;
		int currentMinOrder = NYTNode->order;
		NYTNode->left = new HuffmanNode(-1, 0, currentMinOrder - 2, nullptr, nullptr, NYTNode, true);
		NYTNode->right = new HuffmanNode(symbol, 0, currentMinOrder - 1, nullptr, nullptr, NYTNode, false);
		symbolNode = NYTNode->right;
		NYTNode = NYTNode->left;
	}
			
	SlideAndIncrement(symbolNode);	
}

std::string AdaptiveHuffmanCoding::Encode(int symbol)
{
	HuffmanNode *symbolNode = GetSymbolNode(symbol, root);
	if (symbolNode != nullptr)
	{
		std::string result = GetPathToSymbol(root, symbolNode, "");
		UpdateTreeModel(symbol);
		return result;
	}
	
	std::stringstream ss;
	if(symbol <= 2*r)
		ss << GetPathToSymbol(root, NYTNode, "") << std::bitset<5>(symbol-1);
	else
		ss << GetPathToSymbol(root, NYTNode, "") << std::bitset<4>(symbol-r-1);
	UpdateTreeModel(symbol);
	return ss.str();
}

std::string AdaptiveHuffmanCoding::Encode2(int symbol)
{
	HuffmanNode *symbolNode = GetSymbolNode(symbol, root);
	if (symbolNode != nullptr)
	{
		std::string result = GetPathToSymbol(root, symbolNode, "");
		UpdateTreeModel(symbol);
		return result;
	}
	
	std::stringstream ss;
	ss << GetPathToSymbol(root, NYTNode, "") << std::bitset<9>(symbol);

	UpdateTreeModel(symbol);
	return ss.str();
}


std::string AdaptiveHuffmanCoding::Decode(const std::string& inputStr)
{
	std::string decoded;
	char *Arr1 = new char[4];
	char *Arr2 = new char[5];
	
	int val, index=0;
	bool traversing = false;
	HuffmanNode* currNode;


	while(index!=inputStr.size())
	{
		if(traversing)
		{
			if(inputStr[index]=='1')
			{
				currNode = currNode->right;
			}
			else
			{
				currNode = currNode->left;
			}
			if(currNode->value != -1)
			{
				val = currNode->value;
				decoded += val+'a'-1;
				UpdateTreeModel(val);
				currNode = root;
			}
			if(currNode->isNYT)
				traversing = false;
			index++;
		}
		else
		{
			inputStr.copy(Arr1,e,index);
			val = std::stoi(Arr1, nullptr, 2);
			if(val < r)
			{
				inputStr.copy(Arr2,e+1,index);
				val = std::stoi(Arr2, nullptr, 2) + 1;
				decoded += val + 'a' - 1;
				index += e+1;
				
			}
			else
			{
				inputStr.copy(Arr1,e,index);
				val = std::stoi(Arr1, nullptr, 2) + r + 1;
				decoded += val + 'a' - 1;
				index += e;
			}
			UpdateTreeModel(val);
			traversing = true;
			currNode = root;
		}
	}
	delete Arr1;
	delete Arr2;
	return decoded;
	
}

int AdaptiveHuffmanCoding::Decode2(std::istream &inputStr)
{
	int result = -1;
	HuffmanNode *crr = root;
	while (result == -1)
	{
		if (crr->value != -1)
		{
			result = crr->value;
			UpdateTreeModel(crr->value);
			crr = root;
		}
		else if (crr->isNYT)
		{
			char numberArr[10];
			inputStr.getline(numberArr, 10);
			inputStr.clear();
			int number = std::stoi(numberArr, nullptr, 2);
			result = number;
			UpdateTreeModel(number);
			crr = root;
		}
		else if(inputStr.get() == '0')
		{
			crr = crr->left;
		}
		else crr = crr->right;
	}

	return result;
}

AdaptiveHuffmanCoding::~AdaptiveHuffmanCoding()
{
	DeleteHuffmanTree(root);
}

AdaptiveHuffmanCoding::HuffmanNode * AdaptiveHuffmanCoding::GetSymbolNode(int symbol, HuffmanNode *crr) const
{
	if (crr == nullptr || crr->value == symbol) 
	{
		return crr;
	}

	HuffmanNode *leftResult = GetSymbolNode(symbol, crr->left);
	
	return leftResult == nullptr ? GetSymbolNode(symbol, crr->right) : leftResult;
}

void AdaptiveHuffmanCoding::FindBlockLeader(HuffmanNode *crr, HuffmanNode *&crrMax) const
{
	if (crr == nullptr || crrMax == root)
	{
		return;
	}

	if (crr->weight == crrMax->weight && crr != crrMax->parent  && crr->order > crrMax->order)
	{
		crrMax = crr;
	}

	FindBlockLeader(crr->left, crrMax);
	FindBlockLeader(crr->right, crrMax);
}

void AdaptiveHuffmanCoding::SwapNodes(HuffmanNode *first, HuffmanNode *second)
{
	if (first->parent == nullptr || second->parent == nullptr) return;
	
	if (first->parent == second || second->parent == first) return;

	HuffmanNode *&firstRef = first->parent->left == first ? first->parent->left : first->parent->right;
	HuffmanNode *&secondRef = second->parent->left == second ? second->parent->left : second->parent->right;

	std::swap(firstRef, secondRef);
	std::swap(firstRef->parent, secondRef->parent);
	std::swap(firstRef->order, secondRef->order);
}

std::string AdaptiveHuffmanCoding::GetPathToSymbol(HuffmanNode *crr, HuffmanNode *result, std::string currentPath) const
{
	if (crr == result)
	{
		return currentPath;
	}
	if (crr == nullptr)
	{
		return "";
	}
	
	std::string left = GetPathToSymbol(crr->left, result, currentPath + "0");
	return left != "" ? left : GetPathToSymbol(crr->right, result, currentPath + "1");
}

void AdaptiveHuffmanCoding::SlideAndIncrement(HuffmanNode *node)
{
	if (node == nullptr)
	{
		return;
	}

	HuffmanNode *blockLeader = node;
	FindBlockLeader(root, blockLeader);
	if (blockLeader != node)
	{
		SwapNodes(blockLeader, node);
	}

	++node->weight;
	SlideAndIncrement(node->parent);
}

void AdaptiveHuffmanCoding::DeleteHuffmanTree(HuffmanNode *crrNode)
{
	if (crrNode == nullptr) return;

	DeleteHuffmanTree(crrNode->left);
	DeleteHuffmanTree(crrNode->right);

	delete crrNode;
}

const int AdaptiveHuffmanCoding::PSEUDO_EOF;


void HuffTest(unsigned noOfChar)
{
	// std::cout << "Paragraph count: " + filename << std::endl;
    std::ifstream inputFile("testData");
	if (!inputFile.is_open())
    {
        std::cerr << "Failed to open input file"<<std::endl;
        exit(EXIT_FAILURE);
    }
    AdaptiveHuffmanCoding encoder, decoder;
	std::string encoded, decoded, message;
    std::stringstream encodedStream,decodedStream, messageStream;
	char symbol = inputFile.get();
    for(int i = 0 ; i < noOfChar; i++)
	{
		messageStream << symbol;
        encodedStream << encoder.Encode2(symbol);
		symbol = inputFile.get();
	}
	encodedStream << encoder.Encode2(AdaptiveHuffmanCoding::PSEUDO_EOF);
	encoded = encodedStream.str();

	int symbolAsc = decoder.Decode2(encodedStream);
	while(symbolAsc != AdaptiveHuffmanCoding::PSEUDO_EOF)
	{
		decodedStream << (char)symbolAsc;
		symbolAsc = decoder.Decode2(encodedStream);
	}
	decoded=decodedStream.str();
	message = messageStream.str();
	// assert(message == decoded);
	// std::cout << "Original Message: " << message << std::endl;
	// std::cout << "Encoded Message: " << encoded << std::endl;
	// std::cout << "Decoded Message: " << decoded << std::endl;
	inputFile.close();
}

void parallelHuffTest(unsigned noOfChar, int nthreads)
{
	std::ifstream inputFile("testData");
	if (!inputFile.is_open())
    {
        std::cerr << "Failed to open input file"<<std::endl;
        exit(EXIT_FAILURE);
    }
	std::vector<AdaptiveHuffmanCoding> encoder(nthreads), decoder(nthreads);
	std::vector<std::string> encoded(nthreads), decoded(nthreads), messageArr(nthreads);
	std::stringstream messageStream, encoderStream, decodedStream;
	char symbol = inputFile.get();
    for(int i = 0 ; i < noOfChar; i++)
	{
		messageStream << symbol;
		symbol = inputFile.get();
	}	
	inputFile.close();
	std::string message = messageStream.str();
	for(int i = 0 ; i < message.size();i++)
		messageArr[i%nthreads] += message[i];
	omp_set_num_threads(nthreads);
	#pragma omp parallel private(encoderStream, decodedStream)
	{
		int n = omp_get_thread_num();
		for(int i = 0 ; i < messageArr[n].size() ; i++)
		{
			encoderStream << encoder[n].Encode2(messageArr[n][i]);
		}
		encoderStream << encoder[n].Encode2(AdaptiveHuffmanCoding::PSEUDO_EOF);
		encoded[n] = encoderStream.str();

		int symbolAsc = decoder[n].Decode2(encoderStream);
		while(symbolAsc != AdaptiveHuffmanCoding::PSEUDO_EOF)
		{
			decodedStream << (char)symbolAsc;
			symbolAsc = decoder[n].Decode2(encoderStream);
		}
		decoded[n]=decodedStream.str();

	}
}

void ReadData(int& size, std::vector<uint16_t>& dataArr)
{
	int ch;
    std::string filename;

    std::cout <<  
        "\nChoose dataset : \n" 
        "1. bonsai CT scan (256x256x256) \n"
        "2. aneurism CT scan (256x256x256) \n"
        "3. head aneurism CT scan (512x512x512) \n";
    std::cin >> ch;
    switch(ch)
    {   
        case 1:
            filename = "bonsai_256x256x256_uint8.raw";
            size = 1;
            break;
        case 2:
            filename = "aneurism_256x256x256_uint8.raw";
            size = 1;
            break;
        case 3:
            filename = "vertebra_512x512x512_uint16.raw";
            size = 2;
            break;
    }

    std::ifstream dataBuff(filename, std::ios::binary | std::ios::in);
    if(!dataBuff.is_open())
    {
        std::cout << "Error opening file " << std::endl;
        return;
    }   

    char temp[size];
    int ind = 0;
    while(dataBuff.read(temp, sizeof(temp)))
    {
        dataArr.push_back((uint16_t)temp[0]);
    }
}

void serialImp()
{
	int size;
	std::vector<uint16_t> dataArr, decoded;
	std::cout << "Read Data \n";
	ReadData(size, dataArr);
	#pragma omp parallel num_threads(5)
	{
		AdaptiveHuffmanCoding encoder, decoder;
		std::string encoded, message;
		std::stringstream encodedStream;

		std::cout << "Encode \n";
		#pragma omp for
		for(int i = 0 ; i < dataArr.size() ; i++)
		{
			encodedStream << encoder.Encode2(dataArr[i]);
		}
		std::cout << "decode \n";
		#pragma omp for
		for(int i = 0 ; i < dataArr.size() ; i++)
		{
			decoded.push_back(decoder.Decode2(encodedStream));
		}
	}
	std::cout << "Verify \n";
	int err = 0;
	for(int i = 0 ; i < dataArr.size(); i++)
		if( dataArr[i] != decoded[i])
			err++;
	std::cout << err << std::endl;
}

void parallelImp()
{

}

enum Code{
	serial=1,
	parallel
};

int main()
{
	int ch;
    std::cout << "Serial[1] or Parallel[2] ? ";
    std::cin >> ch ;
    switch (ch)
    {
    case Code::serial:
        serialImp();
        break;
    
    case Code::parallel:
        parallelImp();
        break;

    default:
        break;
    }
	return 0;
}
