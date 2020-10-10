#ifndef TRANSITION_H
#define TRANSITION_H
#include <vector>
#include "State.h"
#include "NumType.h"

class Transition {
public:
    Transition();
    Transition::Transition(State* f, State* t);
    NumType cnt;
    LogNumType logCnt;    
    State* from;
    State* to;
    void add_cnt(NumType Prob);
    void add_log_cnt(LogNumType LogProb);
    void set(State* f, State* t);
};


#endif