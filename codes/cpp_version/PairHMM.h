#ifndef PAIRHMM_H
#define PAIRHMM_H

#include <vector>
#include <string>
#include <iostream>
#include "NumType.h"
#include "DataTool.h"
#include "State.h"
#include "Transition.h"

//typedef double NumType; // To represent Probability numerically;
//typedef double LogNumType;
typedef std::vector<int8_t> proSeqType; //Encoding the ProSeqType, index start from 1
typedef std::pair<std::vector<int8_t>, std::vector<int8_t> > dnaSeqType;
class PairHMM {
public:
    PairHMM();
    void naiveForward(const proSeqType & proSeq, const dnaSeqType & dnaSeq);
    void logForward(const proSeqType & proSeq, const dnaSeqType & dnaSeq);
    void naiveBackward(const proSeqType & proSeq, const dnaSeqType & dnaSeq);
    void logBackward(const proSeqType & proSeq, const dnaSeqType & dnaSeq);
    void forward(const proSeqType & proSeq, const dnaSeqType & dnaSeq, int option);
    void backward(const proSeqType & proSeq, const dnaSeqType & dnaSeq, int option);
    void BaumWelchSingleStep(const proSeqType & proSeq, const dnaSeqType & dnaSeq, int option);
    void BaumWelch(const std::vector<std::string> & rawProSeq, const std::vector<std::string> & rawdnaSeq, int iterTimes, int option);
    void naiveBaumWelch(const std::vector<proSeqType> & proSeqs, const std::vector<dnaSeqType> & dnaSeqs, int iterTimes, int option);
    void initialize();
    void updateTransitions(int option);
    void logUpdateTransitions();
    void naiveUpdateTransitions();
    void updateEmissions(const proSeqType & proSeq, const dnaSeqType & dnaSeq, int option);
    void naiveUpdateEmissions(const proSeqType & proSeq, const dnaSeqType & dnaSeq);
    void logUpdateEmissions(const proSeqType & proSeq, const dnaSeqType & dnaSeq);
    void logUpdatePossibilities();
    void displayParameters(int option);
    NumType startFwd, startBwd, finishFwd, finishBwd, reversep;
    LogNumType logStartFwd, logStartBwd, logFinishFwd, logFinishBwd, logReversep;
    NumType Prob; // To measure the converge of BW algorithm
    LogNumType logProb;
    State start, finish; 
    State D_1, D_2, D_3;
    State I_1, I_2, I_3, I_4, I_5, I_6, I_7;
    State H_1, H_2, H_3, H_4, H_5, H_6, H_7; // H for central hidden node
    State Match;
    NumType omega_d, omega_i, gamma, alpha_d, beta_d, epsilon_d, delta_d, alpha_i, beta_i, epsilon_i, delta_i;
    LogNumType log_omega_d, log_omega_i, log_gamma, log_alpha_d, log_beta_d, log_epsilon_d, log_delta_d, log_alpha_i, log_beta_i, log_epsilon_i, log_delta_i;
    Transition M, A;
    Transition B_i, C_i, D_i, E_i, F_i, G_i, H_i, X_i, J_i, K_i;
    Transition B_d, D_d, E_d, F_d, G_d, H_d, X_d, J_d, K_d;  
    std::vector<NumType> psi, phi, psi_cnt, phi_cnt;  //psi: insertion, phi: deletion
    std::vector<std::vector<NumType> > pi, pi_cnt;
    std::vector<LogNumType> log_psi, log_phi, log_psi_cnt, log_phi_cnt;  //psi: insertion, phi: deletion
    std::vector<std::vector<LogNumType> > log_pi, log_pi_cnt;
    void resetForward(int n, int m);
    void resetBackward(int n, int m);
    void updatePossibilities(int option);
    void naiveUpdatePossibilities();
    void emissionInitialize();
    void transitionInitialize();
    void parameterInitialize();
    void shapeInitialize();
    void BaumWelchSingleStepInitialize(int n, int m, int option);
private:

};

#endif