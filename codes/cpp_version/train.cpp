#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
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
        int curra = rand() % 21;
        if(dt.encodeAA(dt.decodeAA(curra)) != curra) return false;
        if(dt.encodeBase(dt.decodeBase(currb)) != currb) return false;
    }
    return true;
}

std::string genDNA(int n) {
    DataTool dt;
    std::string res;
    std::vector<int> cnts(4, 0);
    srand((unsigned)time(NULL)); 
    for (int i = 0; i < n; i ++) {
        int curr = rand()%4;
        res.push_back(dt.decodeBase(curr));
        cnts[curr] ++;
    }
    //for (auto c:cnts) std::cout << c << " ";
    //std::cout << std::endl;
    return res;
}

std::string genPro(int n) {
    DataTool dt;
    std::string res;
    std::vector<int> cnts(21, 0);
    srand((unsigned)time(NULL)); 
    for (int i = 0; i < n; i ++) {
        int curr = rand() % 21;
        res.push_back(dt.decodeAA(curr));
        cnts[curr] ++;
    }
    //for (auto c:cnts) std::cout << c << " ";
    //std::cout << std::endl;
    return res;
}

std::vector<std::vector<NumType>> readPi(const std::string & filePath) {
    std::ifstream in(filePath, std::ios::in);
    std::vector<std::vector<NumType>> curr(21, std::vector<NumType>(64, 0));
    std::vector<std::vector<NumType>> prob(21, std::vector<NumType>(64, 0));
    for(int i = 0; i < 21; i ++) {
        for (int j = 0; j < 64; j ++) {
            in >> curr[i][j];
        }
    }
    NumType total = 0.0;
    for (auto i: curr) {
        for (auto j : i) total += exp(j);
    }
    for(int i = 0; i < 21; i ++) {
        for (int j = 0; j < 64; j ++) {
            prob[i][j] = exp(curr[i][j]) / total;
        }
    }
    in.close();
    return prob;
}

std::vector<std::vector<NumType>> directReadPi(const std::string & filePath) {
    std::ifstream in(filePath, std::ios::in);
    std::vector<std::vector<NumType>> curr(21, std::vector<NumType>(64, 0));
    for(int i = 0; i < 21; i ++) {
        for (int j = 0; j < 64; j ++) {
            in >> curr[i][j];
        }
    }
    in.close();
    return curr;
}

int main1() 
{
    //std::cout << "hello" << std::endl;
    PairHMM tiny;
    tiny.eps = 1e-5;
    //tiny.default_filepath = "C:\\Users\\InYuo\\Documents\\GitHub\\protein-dna-align-EM\\codes\\cpp_version\\test_para_pg_450_nocon.txt";
    //tiny.default_filepath = "C:\\Users\\InYuo\\Documents\\GitHub\\protein-dna-align-EM\\codes\\cpp_version\\test_for_set_output.txt";
    //tiny.default_filepath = "C:\\Users\\InYuo\\Documents\\GitHub\\protein-dna-align-EM\\codes\\cpp_version\\pg500_500_training.txt";
    //tiny.default_filepath = "C:\\Users\\InYuo\\Documents\\GitHub\\protein-dna-align-EM\\codes\\cpp_version\\pg500_377.txt";
    tiny.default_filepath = "C:\\Users\\InYuo\\Documents\\GitHub\\protein-dna-align-EM\\codes\\cpp_version\\balanced.txt";
    //tiny.setParameters("C:\\Users\\InYuo\\Documents\\GitHub\\protein-dna-align-EM\\codes\\cpp_version\\test_for_set_input.txt");
    //tiny.reNormalize();
    //tiny.naiveTolog();
    //tiny.default_filepath = "C:\\Users\\InYuo\\Documents\\GitHub\\protein-dna-align-EM\\codes\\cpp_version\\parameter_log.txt";
    //tiny.error_filepath = "C:\\Users\\InYuo\\Documents\\GitHub\\protein-dna-align-EM\\codes\\cpp_version\\error_log.txt";
    //tiny.setPi(directReadPi("C:\\Users\\InYuo\\Documents\\GitHub\\protein-dna-align-EM\\codes\\cpp_version\\initProb\\piProb.txt"));
    //for(int i = 0; i < 20; i ++) tiny.phi[i] = 0.995/20.0;
    //tiny.phi[20] = 0.005;
    //tiny.reNormalize();
    tiny.setValid("C:\\Users\\InYuo\\Documents\\GitHub\\protein-dna-align-EM\\data\\pg_500_hg19_presuf10_valid.txt");
    //std::cout << "hello" << std::endl;
    //std::cout << tiny.pi.size() << tiny.pi[0].size() << std::endl;
    //tiny.setInsertion();
    //tiny.displayParameters("start", tiny.default_filepath);
    tiny.testTraining("C:\\Users\\InYuo\\Documents\\GitHub\\protein-dna-align-EM\\data\\pg_500_hg19_presuf10_train.txt", 30, 1);
    //tiny.get_total();
    //tiny.testTraining("C:\\Users\\InYuo\\Documents\\GitHub\\protein-dna-align-EM\\codes\\cpp_version\\test_cnt.txt", 1, 0);
    //DataTool dt;
    //std::string dna = genDNA(100);
    //std::string pro = genPro(50);
    //std::cout << dna << " " << pro << std::endl;
    //proSeqType p = dt.encodePro(pro);
    //dnaSeqType d = dt.encodeDNA(dna);
    //std::vector<proSeqType> testpro {p};
    //std::vector<dnaSeqType> testdna {d};
    //tiny.BaumWelchSingleStep(p, d, 0);
    //tiny.naiveBaumWelch(testpro, testdna, 1, 1);
    // tiny.displayTransitionCnts();
    // std::cout << tiny.overAllProb << std::endl;
    //tiny.BaumWelchSingleStep(p, d, 1);
    //std::cout << tiny.logFinishFwd <<" "<<tiny.logStartBwd<< std::endl;
    //tiny.displayTransitionCnts();
    //tiny.get_total();
    //std::cout << tiny.omega_d << " " << tiny.omega_i << " " << std::endl;
    //tiny.naiveBaumWelch(testpro, testdna, 3, 1);
    //tiny.displayEmissionCnts();
    //std::cout << tiny.omega_d << " " << tiny.omega_i << " " << std::endl;
    //std::cout << tiny.logFinishFwd <<" "<<tiny.logStartBwd<< std::endl;
    //tiny.displayEmissionCnts();
    //tiny.displayTransitionCnts();
    //tiny.checkEmissions();
    //tiny.checkTransitionParameters();
    //tiny.insertionSolver();
    /*for(double i = 0.1; i < 1.0; i += 0.1) {
        tiny.deltaItoObject(i);
    }*/
    //std::cout << dna.size() << std::endl;*/
    return 0;
}

