#ifndef DATATOOL_H
#define DATATOOL_H
#define ori first
#define triplet second
#include <vector>
#include <string>
#include <iostream>
#include "NumType.h"
#include <ctype.h>
#include <unordered_map>

//typedef double NumType; // To represent Probability numerically;
typedef std::vector<int8_t> proSeqType; //Encoding the ProSeqType, index start from 1
typedef std::pair<std::vector<int8_t>, std::vector<int8_t> > dnaSeqType;

// .ori .triplet 

class DataTool {
public:
DataTool();
proSeqType encodePro(const std::string & pro);
dnaSeqType encodeDNA(const std::string & dna);
//private:
int8_t encodeBase(const char & base);
int8_t encodeTriplet(const std::string & dnaTriplet);
int8_t encodeAA(const char & aa);
char decodeAA(int8_t num);
char decodeBase(int8_t num);
std::string decodeTriplet(int num);
std::unordered_map<char, int8_t> check;
std::vector<char> bases;
std::vector<char> aas;
bool checkDNA(std::string & dna);
bool checkPro(std::string & pro);
};

#endif