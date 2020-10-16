#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include "PairHMM.h"
#include "DataTool.h"
#include "State.h"
#include "Transition.h"
#include "NumType.h"
#include <random>
#define ori first
#define triplet second
typedef double ProbType; // To represent Probability numerically;
typedef std::vector<int8_t> proSeqType; //Encoding the ProSeqType, index start from 1
typedef std::pair<std::vector<int8_t>, std::vector<int8_t> > dnaSeqType;

bool testForDt(int n) {
    DataTool dt;
    for(int i = 0; i < n; i ++) {
        int currb = rand() % 4;
        int curra = rand() % 20;
        if(dt.encodeAA(dt.decodeAA(curra)) != curra) return false;
        if(dt.encodeBase(dt.decodeBase(currb)) != currb) return false;
    }
    return true;
}

std::string genDNA(int n) {
    DataTool dt;
    std::string res;
    std::vector<int> cnts(4, 0);
    for (int i = 0; i < n; i ++) {
        int curr = rand()%4;
        res.push_back(dt.decodeBase(curr));
        cnts[curr] ++;
    }
    for (auto c:cnts) std::cout << c << " ";
    std::cout << std::endl;
    return res;
}

std::string genPro(int n) {
    DataTool dt;
    std::string res;
    for (int i = 0; i < n; i ++) {
        res.push_back(dt.decodeAA(rand()%20));
    }
    return res;
}

int main() 
{
    std::cout << "hello" << std::endl;
    PairHMM tiny;
    DataTool dt;
    std::string dna = genDNA(4000);
    std::string pro = genPro(1000);
    proSeqType p = dt.encodePro(pro);
    dnaSeqType d = dt.encodeDNA(dna);
    //std::cout << LogSumExp(log(0), log(0)) << std::endl;
    //std::cout << p.size() << " " << d.ori.size() << " " << std::endl;
    //tiny.BaumWelchSingleStep(p, d, 0);
    //std::cout << tiny.finishFwd << std::endl;
    tiny.displayEmissionCnts();
    tiny.BaumWelchSingleStep(p, d, 1);
    std::cout << tiny.logFinishFwd <<" "<<tiny.logStartBwd<< std::endl;
    //std::cout << tiny.finishFwd <<" "<<tiny.startBwd<< std::endl;
    tiny.displayEmissionCnts();
    std::cout << dna.size() << std::endl;
    return 0;
}

