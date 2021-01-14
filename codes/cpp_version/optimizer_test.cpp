#include "Optimizer.h"
#include <vector>
#include <chrono>
#include <iostream>
#include <random>
using namespace std;

int main() {
    Optimizer opt;
    vector<double> init{0.94, 0.94, 0.86, 0.04, 0.04, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.0};
    opt.setValue(init);
    opt.setCnts(12868.0365843772, 31085.6461251031, 128868.861587414, 1232, 12868.0365843772, 31085.6461251031, 5726.9888229679, 793.076420863528, 786.899273318932, 3322.85247813292, 1617.23707151605, 3197.21276569851, 1579.97569418246, 885.765314250953, 1726.92943881095, 3365.82219152276, 5978.51694458466, 3544.09606901925, 1817.1666302083, 931.401315957345);
    cout << opt.getObject() << endl;
    vector<double> t = opt.gradientDescent(0.00001);
    cout << opt.getObject() << endl;
    
    for (auto i : t) {
        cout << i << endl;
    }
    cout << "hello: "<<opt.check() << endl;
    return 0;
}
