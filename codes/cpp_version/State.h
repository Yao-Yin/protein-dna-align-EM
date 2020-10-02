#ifndef STATE_H
#define STATE_H
#include <vector>

typedef double NumType;

class State {
public:
State();
std::vector<std::vector<NumType> > f, b; // forward[i][j], backward[i][j]
};



#endif