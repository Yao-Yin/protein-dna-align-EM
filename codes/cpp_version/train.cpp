#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include "PairHMM.h"
#include "DataTool.h"
#include "State.h"
#include "Transition.h"
#define ori first
#define triplet second
typedef double ProbType; // To represent Probability numerically;
typedef std::vector<int8_t> proSeqType; //Encoding the ProSeqType, index start from 1
typedef std::pair<std::vector<int8_t>, std::vector<int8_t> > dnaSeqType;

int main() 
{
    PairHMM tiny;
    DataTool dt;
    std::string dna = "ATCGCTCG";
    std::string pro = "ACNDSVNSD";
    proSeqType p = dt.encodePro(pro);
    dnaSeqType d = dt.encodeDNA(dna);
    //std::cout << p.size() << " " << d.ori.size() << " " << std::endl;
    tiny.BaumWelchSingleStep(p, d, 0);
    std::cout << tiny.finishFwd << std::endl;
    return 0;
}

