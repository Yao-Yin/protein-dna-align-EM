#include "PairHMM.h"

PairHMM::PairHMM() {
    initialize();
}

void PairHMM::forward(const proSeqType & proSeq, const dnaSeqType & dnaSeq, int option=0) {
    switch(option) {
        case 0:
            naiveForward(proSeq, dnaSeq);
            break;
    }
    return;
}

void PairHMM::backward(const proSeqType & proSeq, const dnaSeqType & dnaSeq, int option=0) {
    switch(option) {
        case 0:
            naiveBackward(proSeq, dnaSeq);
            break;
    }
    return;
}

void PairHMM::naiveForward(const proSeqType & proSeq, const dnaSeqType & dnaSeq) {
    int N = proSeq.size(); // N = size of protein + 1
    int M = dnaSeq.ori.size(); // M = size of dna + 1
    int n = N - 1;
    int m = M - 1;
    std::cout << N << " " << M << " " << n << " " << m << std::endl;
    //psi: insertion, phi: deletion
    startFwd = NumType(1.0);
    for (int i = 0; i <= n; i ++) {
        for (int j = 0; j <= m; j ++) {
            start.f[i][j] = (i || j) ? NumType(0) : NumType(1);
            D_1.f[i][j] = i ? phi[proSeq[i]]*H_1.f[i-1][j]*omega_d : NumType(0); 
            H_1.f[i][j] = start.f[i][j] + D_1.f[i][j];
            I_1.f[i][j] = j ? psi[dnaSeq.ori[j]]*H_2.f[i][j-1]*omega_i : NumType(0);
            H_2.f[i][j] = H_1.f[i][j]*(NumType(1)-omega_d) + I_1.f[i][j];
            Match.f[i][j] = (i > 0 && j > 2) ? pi[proSeq[i]][dnaSeq.triplet[j]]*H_3.f[i-1][j-3]*gamma : NumType(0); 
            D_3.f[i][j] = i ? phi[proSeq[i]]*H_4.f[i-1][j]*omega_d : NumType(0);
            D_2.f[i][j] =  i ? phi[proSeq[i]]*(H_3.f[i-1][j]*alpha_d + H_6.f[i-1][j]*beta_d) : NumType(0);
            I_2.f[i][j] = j ? psi[dnaSeq.ori[j]]*(D_2.f[i][j-1]*(NumType(1)-delta_d)) : NumType(0);
            I_3.f[i][j] = j ? psi[dnaSeq.ori[j]]*(I_2.f[i][j-1] + H_7.f[i][j-1]*(NumType(1)-epsilon_d)) : NumType(0);
            H_7.f[i][j] = D_2.f[i][j]*delta_d;
            H_6.f[i][j] = H_7.f[i][j]*epsilon_d;
            I_5.f[i][j] = j ? psi[dnaSeq.ori[j]]*I_4.f[i][j-1]*delta_i : NumType(0);
            I_4.f[i][j] = j ? psi[dnaSeq.ori[j]]*(I_6.f[i][j-1]*beta_i + H_3.f[i][j-1]*alpha_i): NumType(0);
            I_6.f[i][j] = j ? psi[dnaSeq.ori[j]]*(I_5.f[i][j-1]*epsilon_i) : NumType(0);
            I_7.f[i][j] = j ? psi[dnaSeq.ori[j]]*H_5.f[i][j-1]*omega_i : NumType(0);
            H_3.f[i][j] = Match.f[i][j] + H_2.f[i][j]*(NumType(1)-omega_i) + H_6.f[i][j]*(NumType(1)-beta_d)
                        + I_3.f[i][j] + I_4.f[i][j]*(NumType(1)-delta_i) + I_5.f[i][j]*(1-epsilon_i) + I_6.f[i][j]*(NumType(1)-beta_i);
            H_4.f[i][j] = H_3.f[i][j]*(NumType(1)-alpha_i-alpha_d-gamma) + D_3.f[i][j];
            H_5.f[i][j] = H_4.f[i][j]*(NumType(1) - omega_d) + I_7.f[i][j];
            finish.f[i][j] = H_5.f[i][j]*(NumType(1)-omega_i);
            //std::cout << i << " " << j << ": " << H_1.f[i][j]<<" "<<H_2.f[i][j]<<" "<<H_3.f[i][j]<<" "<<H_4.f[i][j]<<" "<<H_5.f[i][j]<<std::endl;
        }
    }
    finishFwd = finish.f[n][m];
    reversep = NumType(1.0)/finishFwd;
}

void PairHMM::naiveBackward(const proSeqType & proSeq, const dnaSeqType & dnaSeq) {
    int N = proSeq.size(); // N = size of protein + 1
    int M = dnaSeq.ori.size(); // M = size of dna + 1
    int n = N - 1;
    int m = M - 1;
    finishBwd = NumType(1);
    //psi: insertion, phi: deletion
    for (int i = n; i >= 0; i --) {
        for (int j = m; j >= 0; j --) {
            /* def: just leave
            finish.b[i][j] = (i == n && j == m) ? NumType(1) : NumType(0);
            H_5.b[i][j] = ((j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_7.b[i][j+1]*omega_i) + finish.b[i][j]*(NumType(1)-omega_i);
            I_7.b[i][j] = H_5.b[i][j];
            H_4.b[i][j] = H_5.b[i][j]*(NumType(1)-omega_d) + (i == n ? NumType(0) : phi[proSeq[i+1]]*D_3.b[i+1][j]*omega_d);
            D_3.b[i][j] = H_4.b[i][j];
            H_3.b[i][j] = ((i == n || j >= m - 2) ? NumType(0) : pi[proSeq[i+1]][dnaSeq.triplet[j+3]]*Match.b[i+1][j+3]*gamma) + H_4.b[i][j]*(NumType(1)-gamma-alpha_d-alpha_i)
                        + ((i == n) ? NumType(0) : phi[proSeq[i+1]]*D_2.b[i+1][j]*alpha_d) + ((j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_4.b[i][j+1]*alpha_i);
            Match.b[i][j] = H_3.b[i][j];
            H_6.b[i][j] = H_3.b[i][j]*(NumType(1)-beta_d) + ((i == n) ? NumType(0) : phi[proSeq[i+1]]*D_2.b[i+1][j]*beta_d);
            H_7.b[i][j] = H_6.b[i][j]*epsilon_d + ((j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_3.b[i][j+1]*(NumType(1)-epsilon_d));
            D_2.b[i][j] = H_7.b[i][j]*delta_d + ((j == m) ? NumType(0) : I_2.b[i][j]*(NumType(1)-delta_d));
            I_2.b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_3.b[i][j+1];
            I_3.b[i][j] = H_3.b[i][j];
            I_5.b[i][j] = H_3.b[i][j]*(NumType(1)-epsilon_i) + ((j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_6.b[i][j+1]*epsilon_i);
            I_4.b[i][j] = H_3.b[i][j]*(NumType(1)-delta_i) + ((j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_5.b[i][j+1]*epsilon_i);
            I_6.b[i][j] = H_3.b[i][j]*(NumType(1)-beta_i) + ((j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_4.b[i][j+1]*epsilon_i);
            H_2.b[i][j] = ((j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_1.b[i][j+1]*omega_i) + H_3.b[i][j]*(NumType(1)-omega_i);
            I_1.b[i][j] = H_2.b[i][j];
            H_1.b[i][j] = H_2.b[i][j]*(NumType(1)-omega_d) + (i == n ? NumType(0) : phi[proSeq[i+1]]*D_1.b[i+1][j]*omega_d);
            D_1.b[i][j] = H_1.b[i][j];
            start.b[i][j] = H_1.b[i][j];
            */
           // def: just arrive
           //psi: insertion, phi: deletion
            finish.b[i][j] = (i == n && j == m) ? NumType(1) : NumType(0);
            I_7.b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*H_5.b[i][j+1];
            H_5.b[i][j] = I_7.b[i][j]*omega_i + finish.b[i][j]*(NumType(1)-omega_i);
            D_3.b[i][j] = (i == n) ? NumType(0) : phi[proSeq[i+1]]*H_4.b[i+1][j];
            H_4.b[i][j] = H_5.b[i][j]*(NumType(1)-omega_d) + D_3.b[i][j]*omega_d;
            I_2.b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_3.b[i][j+1];
            Match.b[i][j] = (i == n || j >= m - 2) ? NumType(0) : pi[proSeq[i+1]][dnaSeq.triplet[j+3]]*H_3.b[i+1][j+3];
            D_2.b[i][j] = (i == n) ? NumType(0) : phi[proSeq[i+1]]*(I_2.b[i+1][j]*(NumType(1)-delta_d) + H_7.b[i+1][j]*delta_d);
            I_4.b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*(H_3.b[i][j+1]*(NumType(1)-delta_i)+I_5.b[i][j+1]*epsilon_i);
            H_3.b[i][j] = Match.b[i][j]*gamma + H_4.b[i][j]*(NumType(1)-gamma-alpha_d-alpha_i)
                        + D_2.b[i][j]*alpha_d + I_4.b[i][j]*alpha_i;
            H_6.b[i][j] = H_3.b[i][j]*(NumType(1)-beta_d) + D_2.b[i][j]*beta_d;
            I_3.b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*H_3.b[i][j+1];
            H_7.b[i][j] = H_6.b[i][j]*epsilon_d + I_3.b[i][j]*(NumType(1)-epsilon_d);
            I_5.b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*(H_3.b[i][j+1]*(NumType(1)-epsilon_i)+I_6.b[i][j+1]*epsilon_i);
            I_1.b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*H_2.b[i][j+1];
            I_6.b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*(H_3.b[i][j+1]*(NumType(1)-beta_i)+I_4.b[i][j+1]*epsilon_i);
            H_2.b[i][j] = I_1.b[i][j]*omega_i + H_3.b[i][j]*(NumType(1)-omega_i);
            D_1.b[i][j] = (i == n) ? NumType(0) : phi[proSeq[i+1]]*H_1.b[i+1][j];
            H_1.b[i][j] = D_1.b[i][j]*omega_d + H_2.b[i][j]*(NumType(1)-omega_d);
            start.b[i][j] = H_1.b[i][j];
            //std::cout << i << " " << j << ": " << I_1.b[i][j]<<" "<<I_2.b[i][j]<<" "<<I_3.b[i][j]<<" "<<I_4.b[i][j]<<" "<<I_5.b[i][j]<<" "<<I_6.b[i][j]<<" "<<I_7.b[i][j]<<std::endl;
        }
    }
    startBwd = start.b[0][0];
}

void PairHMM::BaumWelchSingleStep(const proSeqType & proSeq, const dnaSeqType & dnaSeq, int option=0) {
    int n = proSeq.size();
    int m = dnaSeq.ori.size();
    BaumWelchSingleStepInitialize(n, m);
    forward(proSeq, dnaSeq, option);
    backward(proSeq, dnaSeq, option);
    std::cout << finishFwd << " " << startBwd << std::endl;
    //updateTransitions();
    //updateEmissions(proSeq, dnaSeq);
}

void PairHMM::updateTransitions() {
    std::vector<Transition*> vt {&B_i, &C_i, &D_i, &E_i, &F_i, &G_i, &H_i, &X_i, &J_i, &K_i, 
    &B_d, &D_d, &E_d, &F_d, &G_d, &F_d, &H_d, &X_d, &J_d, &K_d, &M, &A};
    for (auto & i: vt) {
        i->add_cnt(reversep);
    }
}

void PairHMM::updateEmissions(const proSeqType & proSeq, const dnaSeqType & dnaSeq) {
    int n = start.f.size() - 1;
    int m = start.f[0].size() - 1; 
    //psi: insertion phi: deletion
    std::vector<NumType> curr_psi_cnt (4, NumType(0));
    std::vector<NumType> curr_phi_cnt (20, NumType(0));
    std::vector<std::vector<NumType> > curr_pi_cnt (20, std::vector<NumType> (64, NumType(0)));
    for (int i = 0; i <= n; i ++) {
        for (int j = 0; j <= m; j ++) {
            for (int k = 0; k < 4; k ++) {
                curr_psi_cnt[k] += dnaSeq.ori[j] == k && j ? I_1.f[i][j]*I_1.b[i][j-1]*psi[k] : NumType(0);
                curr_psi_cnt[k] += dnaSeq.ori[j] == k && j ? I_2.f[i][j]*I_2.b[i][j-1]*psi[k] : NumType(0);
                curr_psi_cnt[k] += dnaSeq.ori[j] == k && j ? I_3.f[i][j]*I_3.b[i][j-1]*psi[k] : NumType(0);
                curr_psi_cnt[k] += dnaSeq.ori[j] == k && j ? I_4.f[i][j]*I_4.b[i][j-1]*psi[k] : NumType(0);
                curr_psi_cnt[k] += dnaSeq.ori[j] == k && j ? I_5.f[i][j]*I_5.b[i][j-1]*psi[k] : NumType(0);
                curr_psi_cnt[k] += dnaSeq.ori[j] == k && j ? I_6.f[i][j]*I_6.b[i][j-1]*psi[k] : NumType(0);
                curr_psi_cnt[k] += dnaSeq.ori[j] == k && j ? I_7.f[i][j]*I_7.b[i][j-1]*psi[k] : NumType(0);
            } 
            for (int k = 0; k < 20; k ++) {
                curr_phi_cnt[k] += proSeq[i] == k && i ? D_1.f[i][j]*D_1.b[i-1][j]*phi[k] : NumType(0);
                curr_phi_cnt[k] += proSeq[i] == k && i ? D_2.f[i][j]*D_2.b[i-1][j]*phi[k] : NumType(0);
                curr_phi_cnt[k] += proSeq[i] == k && i ? D_3.f[i][j]*D_3.b[i-1][j]*phi[k] : NumType(0);
            }
            for (int s = 0; s < 20; s ++) {
                for (int t = 0; t < 64; t ++) {
                    curr_pi_cnt[s][t] += proSeq[i] == s && dnaSeq.triplet[j] == t && i && j > 2 ? 
                                        Match.f[i][j] * Match.b[i-1][j-3] * pi[s][t] : NumType(0);
                }
            }
        }
    }
    for (int k = 0; k < 4; k ++) {
        psi[k] += curr_psi_cnt[k]*reversep;
    }
    for (int k = 0; k < 20; k ++) {
        phi[k] += curr_phi_cnt[k]*reversep;
    }
    for (int s = 0; s < 20; s ++) {
        for (int t = 0; t < 64; t ++) {
            pi[s][t] += curr_pi_cnt[s][t]*reversep;
        }
    }
}

void PairHMM::naiveBaumWelch(const std::vector<proSeqType> & proSeqs, const std::vector<dnaSeqType> & dnaSeqs, int iterTimes) {
    int nums = proSeqs.size();
    for (int iter = 0; iter < iterTimes; iter ++) {
        for (int i = 0; i < nums; i ++) {
            BaumWelchSingleStep(proSeqs[i], dnaSeqs[i]);
        }
        updatePossibilities();
    }
    displayParameters();
}

void PairHMM::naiveUpdatePossibilities() {
    // transition part
    NumType temp = A.cnt + M.cnt + B_d.cnt + D_d.cnt + E_d.cnt + B_i.cnt + D_i.cnt;
    gamma = M.cnt / temp;
    alpha_i = (B_i.cnt + D_i.cnt + E_i.cnt) / temp;
    alpha_d = (B_d.cnt + D_d.cnt + E_d.cnt) / temp;
    beta_i = X_i.cnt / (X_i.cnt + B_i.cnt);
    epsilon_i = (X_i.cnt + B_i.cnt) / (X_i.cnt + B_i.cnt + E_i.cnt);
    delta_i = (X_i.cnt + B_i.cnt + E_i.cnt) / (X_i.cnt + B_i.cnt + E_i.cnt + D_i.cnt);
    beta_d = X_d.cnt / (X_d.cnt + B_d.cnt);
    epsilon_d = (X_d.cnt + B_d.cnt) / (X_d.cnt + B_d.cnt + E_d.cnt);
    delta_d = (X_d.cnt + B_d.cnt + E_d.cnt) / (X_d.cnt + B_d.cnt + E_d.cnt + D_d.cnt);
    omega_d = (J_d.cnt + K_d.cnt) / (temp + temp);
    omega_i = (J_i.cnt + K_i.cnt) / (temp + temp);
    // emmission part
    NumType total_psi(0);
    NumType total_phi(0);
    NumType total_match(0);
    for (int i = 0; i < phi_cnt.size(); i ++) { total_phi += phi_cnt[i]; }
    for (int i = 0; i < psi_cnt.size() ; i ++) { total_psi += psi_cnt[i]; }
    for (int i = 0; i < phi_cnt.size(); i ++) { phi[i] = phi_cnt[i] / total_phi; }
    for (int i = 0; i < psi_cnt.size() ; i ++) { psi[i] = psi_cnt[i] / total_psi; }
    for (int i = 0; i < pi_cnt.size(); i ++) {
        for (int j = 0; j < pi_cnt[0].size(); j ++) {
            total_match += pi_cnt[i][j];
        }
    }
    for (int i = 0; i < 20; i ++) {
        for (int j = 0; j < 64; j ++) {
            pi[i][j] = pi_cnt[i][j] / total_match;
        }
    }
}

void PairHMM::initialize() {
    shapeInitialize();
    transitionInitialize();
    emissionInitialize();
    parameterInitialize();
}

void PairHMM::emissionInitialize() {
    std::fill(psi.begin(), psi.end(), NumType(0.25));
    for (int i = 0; i < 20; i ++) {
        phi[i] = NumType(1.0/20);
    }
    for (int i = 0; i < 20; i ++) {
        for (int j = 0; j < 64; j ++) {
            pi[i][j] = NumType(1.0 / 20 / 64);
        }
    }
}

void PairHMM::transitionInitialize() {
    /*
    &B_i, &C_i, &D_i, &E_i, &F_i, &G_i, &H_i, &X_i, &J_i, &K_i, 
    &B_d, &D_d, &E_d, &F_d, &G_d, &F_d, &H_d, &X_d, &J_d, &K_d, &M, &A
    */
    J_d.set(&H_1, &D_1);
    J_i.set(&H_2, &I_1);
    M.set(&H_3, &Match);
    A.set(&H_3, &H_4);
    K_d.set(&H_4, &D_3);
    K_i.set(&H_5, &I_7);
    F_d.set(&H_3, &D_2);
    X_d.set(&H_6, &D_2);
    B_d.set(&H_6, &H_3);
    D_d.set(&D_2, &I_2);
    E_d.set(&H_7, &I_3);
    G_d.set(&D_2, &H_7);
    H_d.set(&H_7, &H_6);
    B_i.set(&I_6, &H_3);
    E_i.set(&I_5, &H_3);
    D_i.set(&I_4, &H_3);
    F_i.set(&H_3, &I_4);
    G_i.set(&I_4, &I_5);
    H_i.set(&I_5, &I_6);
    X_i.set(&I_6, &I_4);
}

void PairHMM::parameterInitialize() {
    omega_i = NumType(0.9);
    omega_d = NumType(0.9);
    gamma = NumType(0.25);
    //deletion
    alpha_d = NumType(0.25);
    delta_d = NumType(0.5);
    epsilon_d = NumType(0.5);
    beta_d = NumType(0.5);
    //insertion
    alpha_i = NumType(0.25);
    delta_i = NumType(0.5);
    beta_i = NumType(0.5);
    epsilon_i = NumType(0.5);
    //phi, psi
    for (int i = 0; i < phi.size(); i ++) phi[i] = NumType(1.0/phi.size());
    for (int i = 0; i < psi.size(); i ++) psi[i] = NumType(1.0/psi.size());
    for (int i = 0; i < pi.size() ; i ++) {
        for (int j = 0; j < pi[0].size(); j ++) pi[i][j] = NumType(1.0/20/64);
    }
}

void PairHMM::shapeInitialize() {
    phi = std::vector<NumType> (20, NumType(0));
    phi_cnt = std::vector<NumType> (20, NumType(0));
    psi = std::vector<NumType> (4, NumType(0));
    psi_cnt = std::vector<NumType> (4, NumType(0));
    pi = std::vector<std::vector<NumType> > (20, std::vector<NumType>(64, NumType(0)));
    pi_cnt = std::vector<std::vector<NumType> > (20, std::vector<NumType>(64, NumType(0)));
}

void PairHMM::BaumWelchSingleStepInitialize(int n, int m) {
    std::vector<State*> list { &start, &finish, &D_1, &D_2, &D_3, &I_1, &I_2, &I_3, &I_4, &I_5, &I_6, &I_7,
                          &H_1, &H_2, &H_3, &H_4, &H_5, &H_6, &H_7, &Match};
    for(auto & ptr: list) {
        ptr->f = std::vector<std::vector<NumType> >(n, std::vector<NumType> (m, NumType(0)));
        ptr->b = std::vector<std::vector<NumType> >(n, std::vector<NumType> (m, NumType(0)));
    }
}

void PairHMM::displayParameters(){}

void PairHMM::updatePossibilities(){}