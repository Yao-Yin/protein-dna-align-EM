#include "Transition.h"

Transition::Transition(){
    cnt = NumType(0.0);
    logCnt = LogNumType(0.0);
    from = nullptr;
    to = nullptr;
}

Transition::Transition(State* f, State* t){
    cnt = NumType(0.0);
    logCnt = LogNumType(0.0);
    from = nullptr;
    to = nullptr;
    set(f, t);
}

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

void Transition::add_log_cnt(LogNumType LogProb) {
    NumType curr(0);
    int N = from->logf.size();
    int M = from->logf[0].size();
    std::vector<LogNumType> curr_list;
    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < M; j ++) {
            curr_list.emplace_back((from->logf[i][j])+(to->logb[i][j])-LogProb);
        }
    }
    cnt += exp(log_sum_exp(curr_list.begin(), curr_list.end()));
}

void Transition::set(State* f, State* t) {
    this->from = f;
    this->to = t;
}