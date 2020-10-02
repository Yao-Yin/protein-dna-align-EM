#include "Transition.h"

void Transition::add_cnt(NumType reverseP) {
    NumType curr(0);
    int N = from->f.size();
    int M = from->f[0].size();
    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < M; j ++) {
            curr += (from->f[i][j])*(to->b[i][j]);
        }
    }
    cnt += curr*reverseP;
}

void Transition::set(State* f, State* t) {
    this->from = f;
    this->to = t;
}