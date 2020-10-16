#ifndef STATE_H
#define STATE_H
#include <vector>
#include "NumType.h"

class State {
public:
    State();
    std::vector<std::vector<NumType> > f, b; // forward[i][j], backward[i][j]
    std::vector<std::vector<LogNumType> > logf, logb; // forward[i][j], backward[i][j]
    void setSize(int n, int m);
    void setLogSize(int n, int m);
    void setNaiveSize(int n, int m);
};



#endif