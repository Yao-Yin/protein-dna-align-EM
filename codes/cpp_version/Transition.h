#ifndef TRANSITION_H
#define TRANSITION_H
#include <vector>
#include "State.h"
typedef double NumType;

class Transition {
public:
    Transition();
    NumType cnt;
    State* from;
    State* to;
    void add_cnt(NumType reverseP);
    void set(State* f, State* t);
};


#endif