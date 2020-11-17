#include "DataTool.h"
#define ori first
#define tripletSeq second
DataTool::DataTool() {
    int curr = 0;
    aas = std::vector<char> {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','*'};
    /*for (int i = 0; i < 26; i ++) {
        if (i == ('B' - 'A') || i == ('J' - 'A') || i == ('O' - 'A') 
         || i == ('U' - 'A') || i == ('X' - 'A') || i == ('Z' - 'A')) check[i] = -1;
        else {
            aas[curr] = ('A' + i);
            check[i] = curr ++;
        }
    }*/
    for (int i = 0; i < aas.size(); i ++) {
        check[aas[i]] = i;
    }
    bases = std::vector<char> {'T', 'C', 'A', 'G'};
    
}

int8_t DataTool::encodeBase(const char & base) {
    /* 
    Encode the base, might be better to do sth to control invalid inputs.
    */
    char newbase = toupper(base); 
    if (newbase == 'T') return 0;
    else if (newbase == 'C') return 1;
    else if (newbase == 'A') return 2;
    else if (newbase == 'G') return 3;
    return -1;
}

int8_t DataTool::encodeAA(const char & AA) {
    /* 
    Encode the amino acid, might be better to do sth to control invalid inputs.
    */
    return check[AA];
}

int8_t DataTool::encodeTriplet(const std::string & dnaTriplet) {
    /* 
    Encode the dna triplets, regards them as 4-based.
    */
    return (encodeBase(dnaTriplet[0]) << 4) | (encodeBase(dnaTriplet[1]) << 2) | (encodeBase(dnaTriplet[2]));
}

std::string DataTool::decodeTriplet(int num) {
    std::string res;
    res.push_back(decodeBase((num >> 4) & 3));
    res.push_back(decodeBase((num >> 2) & 3));
    res.push_back(decodeBase((num) & 3));
    return res;
}

proSeqType DataTool::encodePro(const std::string & proSeq) {
    proSeqType res;
    res.push_back(-1);
    for (int i = 0; i < proSeq.size(); i ++) res.push_back(encodeAA(proSeq[i]));
    return res;
}

dnaSeqType DataTool::encodeDNA(const std::string & dna) {
    // MUST BE LONGER THAN 2
    dnaSeqType res;
    res.ori.push_back(-1);
    res.tripletSeq.push_back(-1);
    for (int i = 0; i < dna.size(); i ++) {
        res.ori.push_back(encodeBase(dna[i]));
    }
    for (int i = 0; i < dna.size(); i ++) {
        if (i < 2) {
            res.triplet.push_back(-1);
            continue;
        }
        std::string curr_triplet;
        for (int j = i - 2; j <= i; j ++) {
            curr_triplet.push_back(dna[j]);
        }
        res.tripletSeq.push_back(encodeTriplet(curr_triplet));
    }
    return res;
}

char DataTool::decodeAA (int8_t aa) {
    return aas[aa];
}

char DataTool::decodeBase (int8_t base) {
    return bases[base];
}

bool DataTool::checkDNA (std::string & str) {
    for(auto & c: str) {
        if(!isalpha(c)) return false;
        c = toupper(c);
        if(c != 'T' && c != 'A' && c != 'G' && c != 'C') return false;
    }
    return true;
}

bool DataTool::checkPro (std::string & str) {
    for(auto & c: str) {
        if (c == '*') continue;
        if(!isalpha(c)) return false;
        c = toupper(c);
        int i = c - 'A';
        if (i == ('B' - 'A') || i == ('J' - 'A') || i == ('O' - 'A') 
         || i == ('U' - 'A') || i == ('X' - 'A') || i == ('Z' - 'A'))  return false;
    }
    return true;
}