#include "DataTool.h"
#define ori first
#define tripletSeq second
DataTool::DataTool() {
    check = std::vector<int8_t> (26, 0);
    int curr = 0;
    for (int i = 0; i < 26; i ++) {
        if (i == ('B' - 'A') || i == ('J' - 'A') || i == ('O' - 'A') 
         || i == ('U' - 'A') || i == ('X' - 'A') || i == ('Z' - 'A')) check[i] = -1;
        else check[i] = curr ++;
    }
}

int8_t DataTool::encodeBase(const char & base) {
    /* 
    Encode the base, might be better to do sth to control invalid inputs.
    */
    if (base == 'T') return 0;
    else if (base == 'C') return 1;
    else if (base == 'A') return 2;
    else if (base == 'G') return 3;
    return -1;
}

int8_t DataTool::encodeAA(const char & AA) {
    /* 
    Encode the amino acid, might be better to do sth to control invalid inputs.
    */
    return check[AA - 'A'];
}

int8_t DataTool::encodeTriplet(const std::string & dnaTriplet) {
    /* 
    Encode the dna triplets, regards them as 4-based.
    */
    return (encodeBase(dnaTriplet[0]) << 4) & (encodeBase(dnaTriplet[1]) << 2) & (encodeBase(dnaTriplet[2]));
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

