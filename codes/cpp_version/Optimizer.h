#ifndef OPTIMIZER_H
#define OPTIMIZER_H
#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>
#include <algorithm>

class Optimizer {
private:
    double lr, omega_i, omega_d, 
    Gamma, alpha_i, alpha_d, 
    delta_i, delta_d, epsilon_i, epsilon_d, 
    beta_i, beta_d, lambda;
    double M, A;
    double B_i, D_i, Z_i, F_i, G_i, H_i, X_i, J_i, K_i;
    double B_d, D_d, Z_d, F_d, G_d, H_d, X_d, J_d, K_d;
    double global;
    double pdomega_i();
    double pdomega_d();
    double pdGamma();
    double pdalpha_i();
    double pdalpha_d();
    double pddelta_i(); 
    double pddelta_d(); 
    double pdepsilon_i();
    double pdepsilon_d();
    double pdbeta_i();
    double pdbeta_d();
    double pdlambda();
public:
    Optimizer();
    void setCnts(double J_d, double J_i, double M, double A, double K_d, double K_i, double 
        F_d, double X_d, double B_d, double D_d, double E_d, double G_d, double H_d, double
        B_i, double E_i, double D_i, double F_i, double G_i, double H_i, double X_i);
    void setValue(double omega_i, double omega_d, double Gamma, double alpha_i, double alpha_d, double delta_i, double delta_d, 
                double epsilon_i, double epsilon_d, double beta_i, double beta_d, double lambda);
    void setValue(const std::vector<double> & probs);
    void setLearningRate(double alpha);
    std::vector<double> gradientDescent(double eps);
    double getObject();
    double check();
    bool checkValid();
    std::vector<double> grad();
};

#endif