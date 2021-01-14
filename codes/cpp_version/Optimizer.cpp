#include "Optimizer.h"

double Optimizer::pdomega_i(){
    double res = -((2*A)/(1 - (omega_i))) + (J_i)/(omega_i) + (K_i)/(omega_i) + 
  lambda*(3*Gamma*((-(beta_d))*(delta_d)*(epsilon_d) + (omega_d))*
     pow(omega_i,  2) - 3*(omega_d)*((-(beta_d))*(delta_d)*(epsilon_d) + 
      (omega_d))*pow(omega_i, 5) + (alpha_i)*(omega_d)*
     ((-(beta_d))*(delta_d)*(epsilon_d) + (omega_d))*pow(omega_i,  3)*
     ((delta_i)*(1 - (epsilon_i)) + 2*(1 - (beta_i))*(omega_i)) + 
    3*(alpha_i)*(omega_d)*((-(beta_d))*(delta_d)*(epsilon_d) + (omega_d))*
     pow(omega_i,  2)*((1 - (beta_i))*(delta_i)*(epsilon_i) + 
      (delta_i)*(1 - (epsilon_i))*(omega_i) + (1 - (beta_i))*pow(omega_i,  2)) + 
    3*(alpha_d)*(omega_d)*pow(omega_i,  3)*(1 - (beta_d) + 
      (delta_d)*(1 - (epsilon_d))*(omega_i) + (1 - (beta_d))*(delta_d)*
       (epsilon_d)*pow(omega_i,  2)) - 3*(omega_d)*
     ((-(beta_d))*(delta_d)*(epsilon_d) + (omega_d))*pow(omega_i,  2)*
     ((-(beta_i))*(delta_i)*(epsilon_i) + pow(omega_i,  3)) + 
    (alpha_d)*(omega_d)*(omega_i)*((delta_d)*(1 - (epsilon_d)) + 
      2*(1 - (beta_d))*(delta_d)*(epsilon_d)*(omega_i))*
     ((-(beta_i))*(delta_i)*(epsilon_i) + pow(omega_i,  3)) + 
    (alpha_d)*(omega_d)*(1 - (beta_d) + (delta_d)*(1 - (epsilon_d))*
       (omega_i) + (1 - (beta_d))*(delta_d)*(epsilon_d)*pow(omega_i,  2))*
     ((-(beta_i))*(delta_i)*(epsilon_i) + pow(omega_i,  3))) ;
    return res;
    }

double Optimizer::pdomega_d(){
    double res = -((2*A)/(1 - (omega_d))) + (J_d)/(omega_d) + (K_d)/(omega_d) + 
  lambda*((alpha_i)*(omega_d)*pow(omega_i,3)*
     ((1 - (beta_i))*(delta_i)*(epsilon_i) + (delta_i)*(1 - (epsilon_i))*
       (omega_i) + (1 - (beta_i))*pow(omega_i,2)) + 
    (alpha_i)*((-(beta_d))*(delta_d)*(epsilon_d) + (omega_d))*pow(omega_i,3)*
     ((1 - (beta_i))*(delta_i)*(epsilon_i) + (delta_i)*(1 - (epsilon_i))*
       (omega_i) + (1 - (beta_i))*pow(omega_i,2)) + 
    Gamma*((-(beta_i))*(delta_i)*(epsilon_i) + pow(omega_i,3)) - 
    (omega_d)*pow(omega_i,3)*((-(beta_i))*(delta_i)*(epsilon_i) + pow(omega_i,3)) - 
    ((-(beta_d))*(delta_d)*(epsilon_d) + (omega_d))*pow(omega_i, 3)*
     ((-(beta_i))*(delta_i)*(epsilon_i) + pow(omega_i, 3)) + 
    (alpha_d)*(omega_i)*(1 - (beta_d) + (delta_d)*(1 - (epsilon_d))*
       (omega_i) + (1 - (beta_d))*(delta_d)*(epsilon_d)*pow(omega_i,2))*
     ((-(beta_i))*(delta_i)*(epsilon_i) + pow(omega_i,3))); 
    return res;
    }

double Optimizer::pdGamma(){
    double res = M/Gamma - A/(1 - Gamma - alpha_d - alpha_i) + 
    lambda* (-beta_d *delta_d * epsilon_d + 
    omega_d) *(-beta_i *delta_i* epsilon_i + pow(omega_i,3)); 
    return res;
    }

double Optimizer::pdalpha_i(){
    double res = -(A/(1-Gamma-(alpha_d)-(alpha_i)))+
    lambda*(omega_d)*((-(beta_d))*(delta_d)*(epsilon_d)+(omega_d))*pow(omega_i, 3)*((1-(beta_i))*(delta_i)*(epsilon_i)+
    (delta_i)*(1-(epsilon_i))*(omega_i)+(1-(beta_i))*pow(omega_i, 2))+((B_i)+(D_i)+(Z_i))/(alpha_i);
    return res;
    }

double Optimizer::pdalpha_d(){
    double res = -(A/(1-Gamma-(alpha_d)-(alpha_i)))+lambda*(omega_d)*(omega_i)*(1-(beta_d)+(delta_d)*(1-(epsilon_d))*(omega_i)+
    (1-(beta_d))*(delta_d)*(epsilon_d)*pow(omega_i, 2))*((-(beta_i))*(delta_i)*(epsilon_i)+pow(omega_i, 3))+((B_d)+(D_d)+(Z_d))/(alpha_d);
    return res;}

double Optimizer::pddelta_i(){
    double res = 
        -((D_i)/(1-(delta_i)))+lambda*((-Gamma)*(beta_i)*(epsilon_i)*((-(beta_d))*(delta_d)*(epsilon_d)+(omega_d))+
        (beta_i)*(epsilon_i)*(omega_d)*((-(beta_d))*(delta_d)*(epsilon_d)+(omega_d))*pow(omega_i, 3)+(alpha_i)*(omega_d)*((-(beta_d))*(delta_d)*(epsilon_d)+
        (omega_d))*pow(omega_i, 3)*((1-(beta_i))*(epsilon_i)+(1-(epsilon_i))*(omega_i))-(alpha_d)*(beta_i)*(epsilon_i)*(omega_d)*(omega_i)*(1-(beta_d)+(delta_d)*(1-(epsilon_d))*(omega_i)
        +(1-(beta_d))*(delta_d)*(epsilon_d)*pow(omega_i, 2)))+((B_i)+(X_i)+(Z_i))/(delta_i);
    return res;} 

double Optimizer::pddelta_d(){
    double res = -((D_d)/(1-(delta_d)))+lambda*((-(alpha_i))*(beta_d)*(epsilon_d)*(omega_d)*pow(omega_i, 3)*((1-(beta_i))*(delta_i)*(epsilon_i)+(delta_i)*(1-(epsilon_i))*(omega_i)+(1-(beta_i))*pow(omega_i, 2))
    -Gamma*(beta_d)*(epsilon_d)*((-(beta_i))*(delta_i)*(epsilon_i)+pow(omega_i, 3))+
    (beta_d)*(epsilon_d)*(omega_d)*pow(omega_i, 3)*((-(beta_i))*(delta_i)*(epsilon_i)+pow(omega_i, 3))+
    (alpha_d)*(omega_d)*(omega_i)*((1-(epsilon_d))*(omega_i)+(1-(beta_d))*(epsilon_d)*pow(omega_i, 2))*((-(beta_i))*(delta_i)*(epsilon_i)
    +pow(omega_i, 3)))+((B_d)+(X_d)+(Z_d))/(delta_d);
    return res;} 

double Optimizer::pdepsilon_i(){
    double res = lambda*((-Gamma)*(beta_i)*(delta_i)*((-(beta_d))*(delta_d)*(epsilon_d)+(omega_d))+
    (beta_i)*(delta_i)*(omega_d)*((-(beta_d))*(delta_d)*(epsilon_d)+(omega_d))*pow(omega_i, 3)+
    (alpha_i)*(omega_d)*((-(beta_d))*(delta_d)*(epsilon_d)+
    (omega_d))*pow(omega_i, 3)*((1-(beta_i))*(delta_i)-(delta_i)*(omega_i))-
    (alpha_d)*(beta_i)*(delta_i)*(omega_d)*(omega_i)*(1-(beta_d)+(delta_d)*(1-(epsilon_d))*(omega_i)+
    (1-(beta_d))*(delta_d)*(epsilon_d)*pow(omega_i, 2)))+((B_i)+(X_i))/(epsilon_i)-(Z_i)/(1-(epsilon_i));
    return res;}

double Optimizer::pdepsilon_d(){
    double res = lambda*((-(alpha_i))*(beta_d)*(delta_d)*(omega_d)*pow(omega_i, 3)*((1-(beta_i))*(delta_i)*(epsilon_i)+(delta_i)*(1-(epsilon_i))*(omega_i)+
    (1-(beta_i))*pow(omega_i, 2))-Gamma*(beta_d)*(delta_d)*((-(beta_i))*(delta_i)*(epsilon_i)+pow(omega_i, 3))+(beta_d)*(delta_d)*(omega_d)*pow(omega_i, 3)*((-(beta_i))*(delta_i)*(epsilon_i)+
    pow(omega_i, 3))+(alpha_d)*(omega_d)*(omega_i)*((-(delta_d))*(omega_i)+(1-(beta_d))*(delta_d)*pow(omega_i, 2))*((-(beta_i))*(delta_i)*(epsilon_i)+pow(omega_i, 3)))+((B_d)+(X_d))/(epsilon_d)-(Z_d)/(1-(epsilon_d));
    return res;}

double Optimizer::pdbeta_i(){
    double res = -((B_i)/(1-(beta_i)))+lambda*((-Gamma)*(delta_i)*(epsilon_i)*((-(beta_d))*(delta_d)*(epsilon_d)+(omega_d))+(delta_i)*(epsilon_i)*(omega_d)*((-(beta_d))*(delta_d)*(epsilon_d)+(omega_d))*pow(omega_i, 3)+(alpha_i)*(omega_d)*((-(beta_d))*(delta_d)*(epsilon_d)+(omega_d))*pow(omega_i, 3)*((-(delta_i))*(epsilon_i)-pow(omega_i, 2))-(alpha_d)*(delta_i)*(epsilon_i)*(omega_d)*(omega_i)*(1-(beta_d)+(delta_d)*(1-(epsilon_d))*(omega_i)+(1-(beta_d))*(delta_d)*(epsilon_d)*pow(omega_i, 2)))+(X_i)/(beta_i);
    return res;}

double Optimizer::pdbeta_d(){
    double res = -((B_d)/(1-(beta_d)))+lambda*((-(alpha_i))*(delta_d)*(epsilon_d)*(omega_d)*pow(omega_i, 3)*((1-(beta_i))*(delta_i)*(epsilon_i)+(delta_i)*(1-(epsilon_i))*(omega_i)+(1-(beta_i))*pow(omega_i, 2))-Gamma*(delta_d)*(epsilon_d)*((-(beta_i))*(delta_i)*(epsilon_i)+pow(omega_i, 3))+(delta_d)*(epsilon_d)*(omega_d)*pow(omega_i, 3)*((-(beta_i))*(delta_i)*(epsilon_i)+pow(omega_i, 3))+(alpha_d)*(omega_d)*(omega_i)*(-1-(delta_d)*(epsilon_d)*pow(omega_i, 2))*((-(beta_i))*(delta_i)*(epsilon_i)+pow(omega_i, 3)))+(X_d)/(beta_d);
    return res;}

double Optimizer::pdlambda(){
    double res = (alpha_i)*(omega_d)*((-(beta_d))*(delta_d)*(epsilon_d)+(omega_d))*pow(omega_i, 3)*((1-(beta_i))*(delta_i)*(epsilon_i)+(delta_i)*(1-(epsilon_i))*(omega_i)+(1-(beta_i))*pow(omega_i, 2))+Gamma*((-(beta_d))*(delta_d)*(epsilon_d)+(omega_d))*((-(beta_i))*(delta_i)*(epsilon_i)+pow(omega_i, 3))-(omega_d)*((-(beta_d))*(delta_d)*(epsilon_d)+(omega_d))*pow(omega_i, 3)*((-(beta_i))*(delta_i)*(epsilon_i)+pow(omega_i, 3))+(alpha_d)*(omega_d)*(omega_i)*(1-(beta_d)+(delta_d)*(1-(epsilon_d))*(omega_i)+(1-(beta_d))*(delta_d)*(epsilon_d)*pow(omega_i, 2))*((-(beta_i))*(delta_i)*(epsilon_i)+pow(omega_i, 3));
    return res;
    }

std::vector<double> Optimizer::gradientDescent(double realeps) {
    double beta1 = 0.9;
    double beta2 = 0.999;
    double eps = 1e-6;
    std::vector<double> res {omega_i, omega_d, Gamma, alpha_i, alpha_d, 
    delta_i, delta_d, epsilon_i, epsilon_d, beta_i, beta_d, lambda};
    int t = 1;
    std::vector<double> v(12, 0.0), s(12, 0.0), v_bias_error(12, 0.0), s_bias_error(12, 0.0);
    std::vector<double> delta (12, 0.0);
    while (abs(getObject() - global) > realeps) {
        // 
        global = getObject();
        delta = grad();
        //std::cout << global << std::endl;
        for (int i = 0; i < 12; i ++) {
            v[i] = beta1*v[i] + (1 - beta1)*delta[i];
            s[i] = beta2*s[i] + (1 - beta2)*(delta[i]*delta[i]);
            v_bias_error[i] =  v[i] / (1 - pow(beta1, t));
            s_bias_error[i] =  s[i] / (1 - pow(beta2, t));
            res[i] += this->lr*v_bias_error[i] / (std::sqrt(s_bias_error[i]) + eps);
        }
        setValue(res);
        t ++;
    }
    return res;
}

std::vector<double> Optimizer::grad() {
    std::vector<double> res(12, 0);
    res[0] = pdomega_i();
    res[1] = pdomega_d(); 
    res[2] = pdGamma();
    res[3] = pdalpha_i(); 
    res[4] = pdalpha_d(); 
    res[5] = pddelta_i();
    res[6] = pddelta_d();
    res[7] = pdepsilon_i();
    res[8] = pdepsilon_d();
    res[9] = pdbeta_i();
    res[10] = pdbeta_d();
    res[11] = pdlambda();
    return res;
}

Optimizer::Optimizer() {
    lr = 3e-4;
    global = 0.0;
}

void Optimizer::setValue(const std::vector<double> & paras) {
    omega_i = paras[0];
    omega_d = paras[1];
    Gamma = paras[2];
    alpha_i = paras[3];
    alpha_d = paras[4]; 
    delta_i = paras[5];
    delta_d = paras[6];
    epsilon_i = paras[7];
    epsilon_d = paras[8];
    beta_i = paras[9];
    beta_d = paras[10];
    lambda = paras[11];
}

double Optimizer::getObject(){
    double alignTransition = M*log(Gamma) + A*log(1 - Gamma - alpha_d - alpha_i) + (B_i + D_i + Z_i)*log(alpha_i) + (B_d + D_d + Z_d)*log(alpha_d) 
            + B_i*log(1 -beta_i) + X_i*log(beta_i) + Z_i*log(1 -epsilon_i) + (X_i + B_i) * log(epsilon_i) 
            + D_i*(log(1 -delta_i)) + (X_i + B_i + Z_i)*log(delta_i) + B_d*(log(1 -beta_d)) + X_d*log(beta_d)
            + Z_d*log(1 -epsilon_d) + (X_d + B_d)*log(epsilon_d) + D_d*log(1 -delta_d) + (X_d + B_d + Z_d)*log(delta_d);
    double otherTransition = J_i*log(omega_i) + J_d*log(omega_d) + K_i*log(omega_i) + K_d*log(omega_d) + A*(2*log(1 -omega_i) + 2*log(1-omega_d));
    return alignTransition + otherTransition;
}

double Optimizer::check(){
    return Gamma*(pow(omega_i, 3) - delta_i * epsilon_i * beta_i)*(omega_d - delta_d*epsilon_d*beta_d) + alpha_i * ((1-beta_i)*pow(omega_i, 2) + delta_i*(1-epsilon_i)*omega_i + delta_i*epsilon_i*(1-beta_i))*omega_d*pow(omega_i, 3)*(omega_d-delta_d*epsilon_d*beta_d) + alpha_d*((1-beta_d)+delta_d*omega_i*(1-epsilon_d)+pow(omega_i, 2)*epsilon_d*delta_d*(1-beta_d))*omega_d*omega_i*(pow(omega_i, 3) - epsilon_i*delta_i*beta_i) - omega_d*pow(omega_i, 3)*(pow(omega_i, 3) - delta_i *epsilon_i*beta_i)*(omega_d - delta_d*epsilon_d*beta_d);
}

bool Optimizer::checkValid() {
    std::vector<double> para {omega_i, omega_d, Gamma, alpha_i, alpha_d, 
    delta_i, delta_d, epsilon_i, epsilon_d, beta_i, beta_d, lambda};
    for (auto t: para) {
        if (t <= 0 || t >= 1) return false;
    }
    return alpha_i + alpha_d + Gamma < 1;
}

void Optimizer::setCnts(double J_d, double J_i, double M, double A, double K_d, double K_i, double 
        F_d, double X_d, double B_d, double D_d, double E_d, double G_d, double H_d, double
        B_i, double E_i, double D_i, double F_i, double G_i, double H_i, double X_i){
                    this->B_i = B_i;
                    this->D_i = D_i;
                    this->Z_i = E_i;
                    this->F_i = F_i;
                    this->G_i = G_i;
                    this->H_i = H_i;
                    this->X_i = X_i;
                    this->J_i = J_i;
                    this->K_i = K_i;
                    this->B_d = B_d;
                    this->D_d = D_d;
                    this->Z_d = E_d;
                    this->F_d = F_d;
                    this->G_d = G_d;
                    this->H_d = H_d;
                    this->X_d = X_d;
                    this->J_d = J_d;
                    this->K_d = K_d;
                    this->M = M;
                    this->A = A;
                }

