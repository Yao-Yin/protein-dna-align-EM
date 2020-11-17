#include "PairHMM.h"

PairHMM::PairHMM() {
    initialize();
}

PairHMM::~PairHMM() {
    std::vector<State*> list {start, finish, D_1, D_2, D_3, I_1, I_2, I_3, I_4, I_5, I_6, I_7,
                        H_1, H_2, H_3, H_4, H_5, H_6, H_7, Match};
    for (auto & ptr: list) {
        delete ptr;
    }
}

void PairHMM::initialize() {
    shapeInitialize();
    statesBuild();
    transitionInitialize();
    emissionInitialize();
    parameterInitialize();
    naiveTolog();
}

void PairHMM::shapeInitialize() {
    phi = std::vector<NumType> ((int)21, NumType(0));
    phi_cnt = std::vector<NumType> (21, NumType(0));
    psi = std::vector<NumType> (4, NumType(0));
    psi_cnt = std::vector<NumType> (4, NumType(0));
    pi = std::vector<std::vector<NumType> > (21, std::vector<NumType>(64, NumType(0)));
    pi_cnt = std::vector<std::vector<NumType> > (21, std::vector<NumType>(64, NumType(0)));

    log_phi = std::vector<LogNumType> (21, LogNumType(0));
    log_phi_cnt = std::vector<LogNumType> (21, LogNumType(0));
    log_psi = std::vector<LogNumType> (4, LogNumType(0));
    log_psi_cnt = std::vector<LogNumType> (4, LogNumType(0));
    log_pi = std::vector<std::vector<LogNumType> > ((int)21, std::vector<LogNumType>(64, LogNumType(0)));
    log_pi_cnt = std::vector<std::vector<LogNumType> > ((int)21, std::vector<LogNumType>(64, LogNumType(0)));
}

void PairHMM::statesBuild() {
    std::vector<State**> list {&start, &finish, &D_1, &D_2, &D_3, &I_1, &I_2, &I_3, &I_4, &I_5, &I_6, &I_7,
                          &H_1, &H_2, &H_3, &H_4, &H_5, &H_6, &H_7, &Match};
    for (auto & pptr: list) {
        *pptr = new State; 
    }
}

void PairHMM::transitionInitialize() {
    /*
    &B_i, &C_i, &D_i, &E_i, &F_i, &G_i, &H_i, &X_i, &J_i, &K_i, 
    &B_d, &D_d, &E_d, &F_d, &G_d, &F_d, &H_d, &X_d, &J_d, &K_d, &M, &A
    */
    J_d.set(H_1, D_1);
    J_i.set(H_2, I_1);
    M.set(H_3, Match);
    A.set(H_3, H_4);
    K_d.set(H_4, D_3);
    K_i.set(H_5, I_7);
    F_d.set(H_3, D_2);
    X_d.set(H_6, D_2);
    B_d.set(H_6, H_3);
    D_d.set(D_2, I_2);
    E_d.set(H_7, I_3);
    G_d.set(D_2, H_7);
    H_d.set(H_7, H_6);
    B_i.set(I_6, H_3);
    E_i.set(I_5, H_3);
    D_i.set(I_4, H_3);
    F_i.set(H_3, I_4);
    G_i.set(I_4, I_5);
    H_i.set(I_5, I_6);
    X_i.set(I_6, I_4);
}

void PairHMM::emissionInitialize() {
    std::fill(psi.begin(), psi.end(), NumType(0.25));
    //psi = std::vector<NumType> {0.1, 0.2, 0.3, 0.4};
    for (int i = 0; i < phi.size(); i ++) {
        phi[i] = NumType(1.0/21);
    }
    for (int i = 0; i < 21; i ++) {
        for (int j = 0; j < 64; j ++) {
            pi[i][j] = NumType(1.0/21/64);
        }
    }
}

void PairHMM::parameterInitialize() {
    omega_i = NumType(0.9);
    omega_d = NumType(0.9);
    gamma = NumType(0.7);
    //deletion
    alpha_d = NumType(0.1);
    delta_d = NumType(0.1);
    epsilon_d = NumType(0.1);
    beta_d = NumType(0.1);
    //insertion
    alpha_i = NumType(0.1);
    delta_i = NumType(0.1); 
    beta_i = NumType(0.1);
    epsilon_i = NumType(0.1);
}

void PairHMM::forward(const proSeqType & proSeq, const dnaSeqType & dnaSeq, int option) {
    switch(option) {
        case 0:
            naiveForward(proSeq, dnaSeq);
            break;
        case 1:
            logForward(proSeq, dnaSeq);
            break;
    }
    return;
}

void PairHMM::naiveForward(const proSeqType & proSeq, const dnaSeqType & dnaSeq) {
    int N = proSeq.size(); // N = size of protein + 1
    int M = dnaSeq.ori.size(); // M = size of dna + 1
    int n = N - 1;
    int m = M - 1;
    //std::cout << N << " " << M << " " << n << " " << m << std::endl;
    //psi: insertion, phi: deletion
    startFwd = NumType(1.0);
    for (int i = 0; i <= n; i ++) {
        for (int j = 0; j <= m; j ++) {
            start->f[i][j] = (i || j) ? NumType(0) : NumType(1);
            D_1->f[i][j] = i ? phi[proSeq[i]]*H_1->f[i-1][j]*omega_d : NumType(0); 
            H_1->f[i][j] = start->f[i][j] + D_1->f[i][j];
            I_1->f[i][j] = j ? psi[dnaSeq.ori[j]]*H_2->f[i][j-1]*omega_i : NumType(0);
            H_2->f[i][j] = H_1->f[i][j]*(NumType(1)-omega_d) + I_1->f[i][j];
            Match->f[i][j] = (i > 0 && j > 2) ? pi[proSeq[i]][dnaSeq.triplet[j]]*H_3->f[i-1][j-3]*gamma : NumType(0); 
            D_3->f[i][j] = i ? phi[proSeq[i]]*H_4->f[i-1][j]*omega_d : NumType(0);
            D_2->f[i][j] =  i ? phi[proSeq[i]]*(H_3->f[i-1][j]*alpha_d + H_6->f[i-1][j]*beta_d) : NumType(0);
            I_2->f[i][j] = j ? psi[dnaSeq.ori[j]]*(D_2->f[i][j-1]*(NumType(1)-delta_d)) : NumType(0);
            I_3->f[i][j] = j ? psi[dnaSeq.ori[j]]*(I_2->f[i][j-1] + H_7->f[i][j-1]*(NumType(1)-epsilon_d)) : NumType(0);
            H_7->f[i][j] = D_2->f[i][j]*delta_d;
            H_6->f[i][j] = H_7->f[i][j]*epsilon_d;
            I_5->f[i][j] = j ? psi[dnaSeq.ori[j]]*I_4->f[i][j-1]*delta_i : NumType(0);
            I_4->f[i][j] = j ? psi[dnaSeq.ori[j]]*(I_6->f[i][j-1]*beta_i + H_3->f[i][j-1]*alpha_i): NumType(0);
            I_6->f[i][j] = j ? psi[dnaSeq.ori[j]]*(I_5->f[i][j-1]*epsilon_i) : NumType(0);
            I_7->f[i][j] = j ? psi[dnaSeq.ori[j]]*H_5->f[i][j-1]*omega_i : NumType(0);
            H_3->f[i][j] = Match->f[i][j] + H_2->f[i][j]*(NumType(1)-omega_i) + H_6->f[i][j]*(NumType(1)-beta_d)
                        + I_3->f[i][j] + I_4->f[i][j]*(NumType(1)-delta_i) + I_5->f[i][j]*(1-epsilon_i) + I_6->f[i][j]*(NumType(1)-beta_i);
            H_4->f[i][j] = H_3->f[i][j]*(NumType(1.0)-alpha_i-alpha_d-gamma) + D_3->f[i][j];
            H_5->f[i][j] = H_4->f[i][j]*(NumType(1) - omega_d) + I_7->f[i][j];
            finish->f[i][j] = H_5->f[i][j]*(NumType(1)-omega_i);
            //std::cout << i << " " << j << " (naiveForward): " << I_1->f[i][j]<<" "<<I_2->f[i][j]<<" "<<I_3->f[i][j]<<" "<<I_4->f[i][j]<<" "<<I_5->f[i][j]<<" "<<I_6->f[i][j]<<" "<<I_7->f[i][j]<<std::endl;
            //std::cout << i << " " << j << " (naiveForward): " << H_1->f[i][j]<<" "<<H_2->f[i][j]<<" "<<H_3->f[i][j]<<" "<<H_4->f[i][j]<<" "<<H_5->f[i][j]<<" "<<H_6->f[i][j]<<" "<<H_7->f[i][j]<<std::endl;
        }
    }
    finishFwd = finish->f[n][m];
    reversep = NumType(1.0)/finishFwd;
}

void PairHMM::logForward(const proSeqType & proSeq, const dnaSeqType & dnaSeq) {
    int N = proSeq.size(); // N = size of protein + 1
    int M = dnaSeq.ori.size(); // M = size of dna + 1
    int n = N - 1;
    int m = M - 1;
    //std::cout << N << " " << M << " " << n << " " << m << std::endl;
    //psi: insertion, phi: deletion
    LogNumType lzero = log(0);
    LogNumType lone (0);
    logStartFwd = lone;
    for (int i = 0; i <= n; i ++) {
        for (int j = 0; j <= m; j ++) {
            //start->f[i][j] = (i || j) ? NumType(0) : NumType(1);
            start->logf[i][j] = (i || j) ? lzero : lone;
            //D_1->f[i][j] = i ? phi[proSeq[i]]*H_1->f[i-1][j]*omega_d : NumType(0); 
            D_1->logf[i][j] = i ? H_1->logf[i-1][j]+log_omega_d : lzero; 
            //H_1->f[i][j] = start->f[i][j] + D_1->f[i][j];
            //std::cout << i << " " << j << ":: " <<start->logf[i][j]<<" "<<D_1->logf[i][j] << std::endl;
            H_1->logf[i][j] = LogSumExp(start->logf[i][j], D_1->logf[i][j]);
            //std::cout << i << " " << j << ":: " <<H_1->logf[i][j]<<" "<<D_1->logf[i][j] << std::endl;
            //I_1->f[i][j] = j ? psi[dnaSeq.ori[j]]*H_2->f[i][j-1]*omega_i : NumType(0);
            I_1->logf[i][j] = j ? H_2->logf[i][j-1] + log_omega_i : lzero;
            //H_2->f[i][j] = H_1->f[i][j]*(NumType(1)-omega_d) + I_1->f[i][j];
            H_2->logf[i][j] = LogSumExp(H_1->logf[i][j] + Log1m(log_omega_d), I_1->logf[i][j]);
            //Match->f[i][j] = (i > 0 && j > 2) ? pi[proSeq[i]][dnaSeq.triplet[j]]*H_3->f[i-1][j-3]*gamma : NumType(0); 
            Match->logf[i][j] = (i > 0 && j > 2) ? log_pi[proSeq[i]][dnaSeq.triplet[j]]+H_3->logf[i-1][j-3]+log_gamma
                                                  -log_phi[proSeq[i]]-log_psi[dnaSeq.ori[j-2]]-log_psi[dnaSeq.ori[j-1]]-log_psi[dnaSeq.ori[j]]:lzero;
            //D_3->f[i][j] = i ? phi[proSeq[i]]*H_4->f[i-1][j]*omega_d : NumType(0);
            D_3->logf[i][j] = i ? H_4->logf[i-1][j] + log_omega_d : lzero;//2020/10/10
            //D_2->f[i][j] =  i ? phi[proSeq[i]]*(H_3->f[i-1][j]*alpha_d + H_6->f[i-1][j]*beta_d) : NumType(0);
            D_2->logf[i][j] = i ? LogSumExp(H_3->logf[i-1][j]+log_alpha_d, H_6->logf[i-1][j] + log_beta_d) : lzero;
            //I_2->f[i][j] = j ? psi[dnaSeq.ori[j]]*(D_2->f[i][j-1]*(NumType(1)-delta_d)) : NumType(0);
            I_2->logf[i][j] = j ? D_2->logf[i][j-1] + Log1m(log_delta_d) : lzero; 
            //I_3->f[i][j] = j ? psi[dnaSeq.ori[j]]*(I_2->f[i][j-1] + H_7->f[i][j-1]*(NumType(1)-epsilon_d)) : NumType(0);
            I_3->logf[i][j] = j ? LogSumExp(I_2->logf[i][j-1], H_7->logf[i][j-1]+Log1m(log_epsilon_d)) : lzero;
            //H_7->f[i][j] = D_2->f[i][j]*delta_d;
            H_7->logf[i][j] = D_2->logf[i][j] + log_delta_d;
            //H_6->f[i][j] = H_7->f[i][j]*epsilon_d;
            H_6->logf[i][j] = H_7->logf[i][j] + log_epsilon_d;
            //I_5->f[i][j] = j ? psi[dnaSeq.ori[j]]*I_4->f[i][j-1]*delta_i : NumType(0);
            I_5->logf[i][j] = j ? I_4->logf[i][j-1] + log_delta_i : lzero;
            //I_4->f[i][j] = j ? psi[dnaSeq.ori[j]]*(I_6->f[i][j-1]*beta_i + H_3->f[i][j-1]*alpha_i): NumType(0);
            I_4->logf[i][j] = j ? LogSumExp(I_6->logf[i][j-1]+log_beta_i, H_3->logf[i][j-1]+log_alpha_i) : lzero;
            //I_6->f[i][j] = j ? psi[dnaSeq.ori[j]]*(I_5->f[i][j-1]*epsilon_i) : NumType(0);
            I_6->logf[i][j] = j ? I_5->logf[i][j-1]+log_epsilon_i : lzero;
            //I_7->f[i][j] = j ? psi[dnaSeq.ori[j]]*H_5->f[i][j-1]*omega_i : NumType(0);
            I_7->logf[i][j] = j ? H_5->logf[i][j-1]+log_omega_i : lzero;
            //H_3->f[i][j] = Match->f[i][j] + H_2->f[i][j]*(NumType(1)-omega_i) + H_6->f[i][j]*(NumType(1)-beta_d)
            //            + I_3->f[i][j] + I_4->f[i][j]*(NumType(1)-delta_i) + I_5->f[i][j]*(1-epsilon_i) + I_6->f[i][j]*(NumType(1)-beta_i);
            H_3->logf[i][j] = LogSumExp(Match->logf[i][j], H_2->logf[i][j]+Log1m(log_omega_i), H_6->logf[i][j]+Log1m(log_beta_d), I_3->logf[i][j],
                                    I_4->logf[i][j]+Log1m(log_delta_i), I_5->logf[i][j]+Log1m(log_epsilon_i), I_6->logf[i][j]+Log1m(log_beta_i));
            //H_4->f[i][j] = H_3->f[i][j]*(NumType(1)-alpha_i-alpha_d-gamma) + D_3->f[i][j];
            H_4->logf[i][j] = LogSumExp(H_3->logf[i][j]+Log1m(LogSumExp(log_alpha_i, log_alpha_d, log_gamma)), D_3->logf[i][j]);
            //H_5->f[i][j] = H_4->f[i][j]*(NumType(1) - omega_d) + I_7->f[i][j];
            H_5->logf[i][j] = LogSumExp(H_4->logf[i][j]+Log1m(log_omega_d),I_7->logf[i][j]);
            //finish->f[i][j] = H_5->f[i][j]*(NumType(1)-omega_i);
            finish->logf[i][j] = H_5->logf[i][j]+Log1m(log_omega_i);
            //checkforward(i, j);
        }
    }
    logFinishFwd = finish->logf[n][m];
    logProb = logFinishFwd;
    //std::cout << logProb <<std::endl;
}

void PairHMM::backward(const proSeqType & proSeq, const dnaSeqType & dnaSeq, int option) {
    switch(option) {
        case 0:
            naiveBackward(proSeq, dnaSeq);
            break;
        case 1:
            logBackward(proSeq, dnaSeq);
            break;
    }
    return;
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
            finish->b[i][j] = (i == n && j == m) ? NumType(1) : NumType(0);
            H_5->b[i][j] = ((j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_7->b[i][j+1]*omega_i) + finish->b[i][j]*(NumType(1)-omega_i);
            I_7->b[i][j] = H_5->b[i][j];
            H_4->b[i][j] = H_5->b[i][j]*(NumType(1)-omega_d) + (i == n ? NumType(0) : phi[proSeq[i+1]]*D_3->b[i+1][j]*omega_d);
            D_3->b[i][j] = H_4->b[i][j];
            H_3->b[i][j] = ((i == n || j >= m - 2) ? NumType(0) : pi[proSeq[i+1]][dnaSeq.triplet[j+3]]*Match->b[i+1][j+3]*gamma) + H_4->b[i][j]*(NumType(1)-gamma-alpha_d-alpha_i)
                        + ((i == n) ? NumType(0) : phi[proSeq[i+1]]*D_2->b[i+1][j]*alpha_d) + ((j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_4->b[i][j+1]*alpha_i);
            Match->b[i][j] = H_3->b[i][j];
            H_6->b[i][j] = H_3->b[i][j]*(NumType(1)-beta_d) + ((i == n) ? NumType(0) : phi[proSeq[i+1]]*D_2->b[i+1][j]*beta_d);
            H_7->b[i][j] = H_6->b[i][j]*epsilon_d + ((j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_3->b[i][j+1]*(NumType(1)-epsilon_d));
            D_2->b[i][j] = H_7->b[i][j]*delta_d + ((j == m) ? NumType(0) : I_2->b[i][j]*(NumType(1)-delta_d));
            I_2->b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_3->b[i][j+1];
            I_3->b[i][j] = H_3->b[i][j];
            I_5->b[i][j] = H_3->b[i][j]*(NumType(1)-epsilon_i) + ((j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_6->b[i][j+1]*epsilon_i);
            I_4->b[i][j] = H_3->b[i][j]*(NumType(1)-delta_i) + ((j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_5->b[i][j+1]*epsilon_i);
            I_6->b[i][j] = H_3->b[i][j]*(NumType(1)-beta_i) + ((j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_4->b[i][j+1]*epsilon_i);
            H_2->b[i][j] = ((j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_1->b[i][j+1]*omega_i) + H_3->b[i][j]*(NumType(1)-omega_i);
            I_1->b[i][j] = H_2->b[i][j];
            H_1->b[i][j] = H_2->b[i][j]*(NumType(1)-omega_d) + (i == n ? NumType(0) : phi[proSeq[i+1]]*D_1->b[i+1][j]*omega_d);
            D_1->b[i][j] = H_1->b[i][j];
            start->b[i][j] = H_1->b[i][j];
            */
           // def: just arrive
           //psi: insertion, phi: deletion
            finish->b[i][j] = (i == n && j == m) ? NumType(1) : NumType(0);
            I_7->b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*H_5->b[i][j+1];
            H_5->b[i][j] = I_7->b[i][j]*omega_i + finish->b[i][j]*(NumType(1)-omega_i);
            D_3->b[i][j] = (i == n) ? NumType(0) : phi[proSeq[i+1]]*H_4->b[i+1][j];
            H_4->b[i][j] = H_5->b[i][j]*(NumType(1)-omega_d) + D_3->b[i][j]*omega_d;
            I_2->b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_3->b[i][j+1];
            Match->b[i][j] = (i == n || j >= m - 2) ? NumType(0) : pi[proSeq[i+1]][dnaSeq.triplet[j+3]]*H_3->b[i+1][j+3];
            D_2->b[i][j] = (i == n) ? NumType(0) : phi[proSeq[i+1]]*(I_2->b[i+1][j]*(NumType(1)-delta_d) + H_7->b[i+1][j]*delta_d);
            I_4->b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*(H_3->b[i][j+1]*(NumType(1)-delta_i)+I_5->b[i][j+1]*epsilon_i);
            H_3->b[i][j] = Match->b[i][j]*gamma + H_4->b[i][j]*(NumType(1)-gamma-alpha_d-alpha_i)
                        + D_2->b[i][j]*alpha_d + I_4->b[i][j]*alpha_i;
            H_6->b[i][j] = H_3->b[i][j]*(NumType(1)-beta_d) + D_2->b[i][j]*beta_d;
            I_3->b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*H_3->b[i][j+1];
            H_7->b[i][j] = H_6->b[i][j]*epsilon_d + I_3->b[i][j]*(NumType(1)-epsilon_d);
            I_5->b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*(H_3->b[i][j+1]*(NumType(1)-epsilon_i)+I_6->b[i][j+1]*epsilon_i);
            I_1->b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*H_2->b[i][j+1];
            I_6->b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*(H_3->b[i][j+1]*(NumType(1)-beta_i)+I_4->b[i][j+1]*epsilon_i);
            H_2->b[i][j] = I_1->b[i][j]*omega_i + H_3->b[i][j]*(NumType(1)-omega_i);
            D_1->b[i][j] = (i == n) ? NumType(0) : phi[proSeq[i+1]]*H_1->b[i+1][j];
            H_1->b[i][j] = D_1->b[i][j]*omega_d + H_2->b[i][j]*(NumType(1)-omega_d);
            start->b[i][j] = H_1->b[i][j];
            //std::cout << i << " " << j << " (naiveBack): " << I_1->b[i][j]<<" "<<I_2->b[i][j]<<" "<<I_3->b[i][j]<<" "<<I_4->b[i][j]<<" "<<I_5->b[i][j]<<" "<<I_6->b[i][j]<<" "<<I_7->b[i][j]<<std::endl;
            //std::cout << i << " " << j << " (naiveBack): " << H_1->b[i][j]<<" "<<H_2->b[i][j]<<" "<<H_3->b[i][j]<<" "<<H_4->b[i][j]<<" "<<H_5->b[i][j]<<" "<<H_6->b[i][j]<<" "<<H_7->b[i][j]<<std::endl;
            //std::cout << i << " " << j << ": " << I_1->b[i][j]<<" "<<I_2->b[i][j]<<" "<<I_3->b[i][j]<<" "<<I_4->b[i][j]<<" "<<I_5->b[i][j]<<" "<<I_6->b[i][j]<<" "<<I_7->b[i][j]<<std::endl;
        }
    }
    startBwd = start->b[0][0];
}

void PairHMM::logBackward(const proSeqType & proSeq, const dnaSeqType & dnaSeq) {
    int N = proSeq.size(); // N = size of protein + 1
    int M = dnaSeq.ori.size(); // M = size of dna + 1
    int n = N - 1;
    int m = M - 1;
    //std::cout << N << " " << M << " " << n << " " << m << std::endl;
    //psi: insertion, phi: deletion
    LogNumType lzero = log(0);
    LogNumType lone (0);
    logFinishBwd = lone;
    //psi: insertion, phi: deletion
    for (int i = n; i >= 0; i --) {
        for (int j = m; j >= 0; j --) {
            // def: just arrive
            //psi: insertion, phi: deletion
            //finish->b[i][j] = (i == n && j == m) ? NumType(1) : NumType(0);
            finish->logb[i][j] = (i == n && j == m) ? lone : lzero;
            //I_7->b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*H_5->b[i][j+1];
            I_7->logb[i][j] = (j == m) ? lzero : H_5->logb[i][j+1];
            //H_5->b[i][j] = I_7->b[i][j]*omega_i + finish->b[i][j]*(NumType(1)-omega_i);
            H_5->logb[i][j] = LogSumExp(I_7->logb[i][j]+log_omega_i, finish->logb[i][j]+Log1m(log_omega_i));
            //D_3->b[i][j] = (i == n) ? NumType(0) : phi[proSeq[i+1]]*H_4->b[i+1][j];
            D_3->logb[i][j] = (i == n) ? lzero : H_4->logb[i+1][j];
            //H_4->b[i][j] = H_5->b[i][j]*(NumType(1)-omega_d) + D_3->b[i][j]*omega_d;
            H_4->logb[i][j] = LogSumExp(H_5->logb[i][j]+Log1m(log_omega_d), D_3->logb[i][j]+log_omega_d);
            //I_2->b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*I_3->b[i][j+1];
            I_2->logb[i][j] = (j == m) ? lzero : I_3->logb[i][j+1];
            //Match->b[i][j] = (i == n || j >= m - 2) ? NumType(0) : pi[proSeq[i+1]][dnaSeq.triplet[j+3]]*H_3->b[i+1][j+3];
            Match->logb[i][j] = (i == n || j >= m - 2) ? lzero : log_pi[proSeq[i+1]][dnaSeq.triplet[j+3]]+H_3->logb[i+1][j+3]
                                                                - log_phi[proSeq[i+1]] - log_psi[dnaSeq.ori[j+1]] - log_psi[dnaSeq.ori[j+2]] - log_psi[dnaSeq.ori[j+3]];
            //D_2->b[i][j] = (i == n) ? NumType(0) : phi[proSeq[i+1]]*(I_2->b[i+1][j]*(NumType(1)-delta_d) + H_7->b[i+1][j]*delta_d);
            D_2->logb[i][j] = (i == n) ? lzero : LogSumExp(I_2->logb[i+1][j]+Log1m(log_delta_d), H_7->logb[i+1][j]+log_delta_d);
            //I_4->b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*(H_3->b[i][j+1]*(NumType(1)-delta_i)+I_5->b[i][j+1]*epsilon_i);
            I_4->logb[i][j] = (j == m) ? lzero : LogSumExp(H_3->logb[i][j+1]+Log1m(log_delta_i), I_5->logb[i][j+1]+log_epsilon_i);
            //H_3->b[i][j] = Match->b[i][j]*gamma + H_4->b[i][j]*(NumType(1)-gamma-alpha_d-alpha_i)
            //            + D_2->b[i][j]*alpha_d + I_4->b[i][j]*alpha_i;
            H_3->logb[i][j] = LogSumExp(Match->logb[i][j]+log_gamma, H_4->logb[i][j]+log(1.0-gamma-alpha_d-alpha_i),/*Log1m(LogSumExp(log_gamma, log_alpha_d, log_alpha_i))*/
                                       D_2->logb[i][j]+log_alpha_d, I_4->logb[i][j]+log_alpha_i);
            //H_6->b[i][j] = H_3->b[i][j]*(NumType(1)-beta_d) + D_2->b[i][j]*beta_d;
            H_6->logb[i][j] = LogSumExp(H_3->logb[i][j]+Log1m(log_beta_d), D_2->logb[i][j]+log_beta_d);
            //I_3->b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*H_3->b[i][j+1];
            I_3->logb[i][j] = (j == m) ? lzero : H_3->logb[i][j+1];
            //H_7->b[i][j] = H_6->b[i][j]*epsilon_d + I_3->b[i][j]*(NumType(1)-epsilon_d);
            H_7->logb[i][j] = LogSumExp(H_6->logb[i][j]+log_epsilon_d, I_3->logb[i][j]+Log1m(log_epsilon_d));
            //I_5->b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*(H_3->b[i][j+1]*(NumType(1)-epsilon_i)+I_6->b[i][j+1]*epsilon_i);
            I_5->logb[i][j] = (j == m) ? lzero : LogSumExp(H_3->logb[i][j+1]+Log1m(log_epsilon_i), I_6->logb[i][j+1]+log_epsilon_i);
            //I_1->b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*H_2->b[i][j+1];
            I_1->logb[i][j] = (j == m) ? lzero : H_2->logb[i][j+1];
            //I_6->b[i][j] = (j == m) ? NumType(0) : psi[dnaSeq.ori[j+1]]*(H_3->b[i][j+1]*(NumType(1)-beta_i)+I_4->b[i][j+1]*epsilon_i);
            I_6->logb[i][j] = (j == m) ? lzero : LogSumExp(H_3->logb[i][j+1]+Log1m(log_beta_i), I_4->logb[i][j+1]+log_epsilon_i);
            //H_2->b[i][j] = I_1->b[i][j]*omega_i + H_3->b[i][j]*(NumType(1)-omega_i);
            H_2->logb[i][j] = LogSumExp(I_1->logb[i][j]+log_omega_i, H_3->logb[i][j]+Log1m(log_omega_i));
            //D_1->b[i][j] = (i == n) ? NumType(0) : phi[proSeq[i+1]]*H_1->b[i+1][j];
            D_1->logb[i][j] = (i == n) ? lzero : H_1->logb[i+1][j];
            //H_1->b[i][j] = D_1->b[i][j]*omega_d + H_2->b[i][j]*(NumType(1)-omega_d);
            H_1->logb[i][j] = LogSumExp(D_1->logb[i][j]+log_omega_d, H_2->logb[i][j]+Log1m(log_omega_d));
            //start->b[i][j] = H_1->b[i][j];
            start->logb[i][j] = H_1->logb[i][j];
            //std::cout << i << " " << j << " (logBack): " << I_1->logb[i][j]<<" "<<I_2->logb[i][j]<<" "<<I_3->logb[i][j]<<" "<<I_4->logb[i][j]<<" "<<I_5->logb[i][j]<<" "<<I_6->logb[i][j]<<" "<<I_7->logb[i][j]<<std::endl;
            //std::cout << i << " " << j << " (logBack): " << H_1->logb[i][j]<<" "<<H_2->logb[i][j]<<" "<<H_3->logb[i][j]<<" "<<H_4->logb[i][j]<<" "<<H_5->logb[i][j]<<" "<<H_6->logb[i][j]<<" "<<H_7->logb[i][j]<<std::endl;
            //checkback(i, j, n, m);
        }
    }
    logStartBwd = start->logb[0][0];
}

void PairHMM::naiveBaumWelch(const std::vector<proSeqType> & proSeqs, const std::vector<dnaSeqType> & dnaSeqs, int iterTimes, int option) {
    int nums = proSeqs.size();
    for (int iter = 0; iter < iterTimes; iter ++) {
        std::cout << "This is the " << iter <<" epoch. " << std::endl;
        naiveTolog();
        std::vector<Transition*> vt {
            &J_d, &J_i, &M, &A, &K_d, &K_i, 
            &F_d, &X_d, &B_d, &D_d, &E_d, &G_d, &H_d,
            &B_i, &E_i, &D_i, &F_i, &G_i, &H_i, &X_i
        };
        for (auto & ptr: vt) ptr->cnt = 0;
        for (int i = 0; i < phi_cnt.size(); i ++) phi_cnt[i] = 0;
        for (int i = 0; i < psi_cnt.size(); i ++) psi_cnt[i] = 0;
        for (int s = 0; s < pi_cnt.size(); s ++) {
            for (int t = 0; t < pi_cnt[0].size(); t ++) pi_cnt[s][t] = 0;
        }   
        displayParameters("Before the " + std::to_string(iter) + " epoch: ", default_filepath);
        for (int i = 0; i < nums; i ++) {
            BaumWelchSingleStep(proSeqs[i], dnaSeqs[i], 1);
        }
        get_total();
        int n = 1;
        //pseudocount(n);
        updateProbabilities(option);
        double tot_psi = 0;
        for (int i = 0; i < psi.size(); i ++) tot_psi += psi[i];
        double tot_phi = 0;
        for (int i = 0; i < phi.size(); i ++) tot_phi += phi[i];
        double tot_pi = 0;
        for (int i = 0; i < pi.size(); i ++) {
            for (int j = 0; j < pi[0].size(); j ++) {
                tot_pi += pi[i][j];
            }
        }
        std::cout <<"check: "<<tot_psi<<" "<<tot_phi<<" "<<tot_pi  << std::endl;
        displayParameters("This is the " + std::to_string(iter) + " epoch: ", default_filepath);
    }
}

void PairHMM::BaumWelchSingleStep(const proSeqType & proSeq, const dnaSeqType & dnaSeq, int option) {
    int n = proSeq.size();
    int m = dnaSeq.ori.size();
    //std::cout << "Start !" << std::endl;
    //time_t curr = clock();
    BaumWelchSingleStepInitialize(n, m, option);
    //std::cout << "Initialize step takes " << (clock() - curr )/CLOCKS_PER_SEC << "sec"<<std::endl;
    //curr = clock();
    forward(proSeq, dnaSeq, option);
    //std::cout << "Forward step takes " << (clock() - curr )/CLOCKS_PER_SEC << "sec"<<std::endl;
    //curr = clock();
    backward(proSeq, dnaSeq, option);
    //std::cout << "Backward step takes " << (clock() - curr )/CLOCKS_PER_SEC << "sec"<<std::endl;
    //curr = clock();
    updateTransitions(option);
    //std::cout << "update transitions step takes " << (clock() - curr )/CLOCKS_PER_SEC << "sec"<<std::endl;
    //curr = clock();
    updateEmissions(proSeq, dnaSeq, option);
    //std::cout << "updatee emissions step takes " << (clock() - curr )/CLOCKS_PER_SEC << "sec"<<std::endl;
    //curr = clock();
}

void PairHMM::BaumWelchSingleStepInitialize(int n, int m, int option) {
    std::vector<State*> list { start, finish, D_1, D_2, D_3, I_1, I_2, I_3, I_4, I_5, I_6, I_7,
                          H_1, H_2, H_3, H_4, H_5, H_6, H_7, Match};
    switch (option) {
        case (0):
            for(auto & ptr: list) {
                //ptr->f.resize(n, std::vector<NumType> (m, NumType(0)));
                ptr->f= std::vector<std::vector<NumType> >(n, std::vector<NumType> (m, NumType(0)));
                //ptr->b.resize(n, std::vector<NumType> (m, NumType(0)));
                ptr->b = std::vector<std::vector<NumType> >(n, std::vector<NumType> (m, NumType(0)));
            }
            break;
        case (1):
                LogNumType lzero (log(0.0));
                for(auto & ptr: list) {
                    //ptr->logf.resize(n, std::vector<LogNumType> (m, LogNumType(0)));
                    //ptr->logb.resize(n, std::vector<LogNumType> (m, LogNumType(0)));
                    ptr->logf = std::vector<std::vector<LogNumType> >(n, std::vector<LogNumType> (m, lzero));
                    ptr->logb = std::vector<std::vector<LogNumType> >(n, std::vector<LogNumType> (m, lzero));
                }
            break;
    }
}

void PairHMM::updateTransitions(int option) {
    int idx = 0;
    std::vector<Transition*> vt {
        &J_d, &J_i, &M, &A, &K_d, &K_i, 
        &F_d, &X_d, &B_d, &D_d, &E_d, &G_d, &H_d,
        &F_i, &X_i, &B_i, &D_i, &E_i, &G_i, &H_i
    };
    switch (option)
    {
    case 0:
        for (auto & i: vt) {
            //std::cout << idx ++ << std::endl;
            i->add_cnt(reversep);
        }
        break;
    
    case 1:
        /*
        J_d.set(H_1, D_1);
        J_i.set(H_2, I_1);
        M.set(H_3, Match);
        A.set(H_3, H_4);
        K_d.set(H_4, D_3);
        K_i.set(H_5, I_7);
        F_d.set(H_3, D_2);
        X_d.set(H_6, D_2);
        B_d.set(H_6, H_3);
        D_d.set(D_2, I_2);
        E_d.set(H_7, I_3);
        G_d.set(D_2, H_7);
        H_d.set(H_7, H_6);
        B_i.set(I_6, H_3);
        E_i.set(I_5, H_3);
        D_i.set(I_4, H_3);
        F_i.set(H_3, I_4);
        G_i.set(I_4, I_5);
        H_i.set(I_5, I_6);
        X_i.set(I_6, I_4);
            break;
        }*/
        J_d.add_log_cnt(logProb, log_omega_d);
        J_i.add_log_cnt(logProb, log_omega_i);
        M.add_log_cnt(logProb, log_gamma);
        A.add_log_cnt(logProb, Log1m(LogSumExp(log_gamma, log_alpha_d, log_alpha_i)));
        K_d.add_log_cnt(logProb, log_omega_d);
        K_i.add_log_cnt(logProb, log_omega_i);
        F_d.add_log_cnt(logProb, log_alpha_d);
        X_d.add_log_cnt(logProb, log_beta_d);
        B_d.add_log_cnt(logProb, Log1m(log_beta_d));
        D_d.add_log_cnt(logProb, Log1m(log_delta_d));
        E_d.add_log_cnt(logProb, Log1m(log_epsilon_d));
        G_d.add_log_cnt(logProb, Log1m(log_delta_d));
        H_d.add_log_cnt(logProb, log_epsilon_d);
        F_i.add_log_cnt(logProb, log_alpha_i);
        X_i.add_log_cnt(logProb, log_beta_i);
        B_i.add_log_cnt(logProb, Log1m(log_beta_i));
        D_i.add_log_cnt(logProb, Log1m(log_delta_i));
        E_i.add_log_cnt(logProb, Log1m(log_epsilon_i));
        G_i.add_log_cnt(logProb, log_delta_i);
        H_i.add_log_cnt(logProb, log_epsilon_i);
        break;
    }
}

void PairHMM::updateEmissions(const proSeqType & proSeq, const dnaSeqType & dnaSeq, int option) {
    switch (option) {
    case 0:
        naiveUpdateEmissions(proSeq, dnaSeq);
        break;
    
    case 1:
        logUpdateEmissions(proSeq, dnaSeq);
        break;
    }
}

void PairHMM::naiveUpdateEmissions(const proSeqType & proSeq, const dnaSeqType & dnaSeq) {
    int n = start->f.size() - 1;
    int m = start->f[0].size() - 1; 
    //psi: insertion phi: deletion
    std::vector<NumType> curr_psi_cnt (4, NumType(0));
    std::vector<NumType> curr_phi_cnt (21, NumType(0));
    std::vector<std::vector<NumType> > curr_pi_cnt (21, std::vector<NumType> (64, NumType(0)));
    for (int i = 0; i <= n; i ++) {
        for (int j = 0; j <= m; j ++) {
            for (int k = 0; k < 4; k ++) {
                curr_psi_cnt[k] += dnaSeq.ori[j] == k && j ? I_1->f[i][j]*I_1->b[i][j-1]*psi[k] : NumType(0);
                curr_psi_cnt[k] += dnaSeq.ori[j] == k && j ? I_2->f[i][j]*I_2->b[i][j-1]*psi[k] : NumType(0);
                curr_psi_cnt[k] += dnaSeq.ori[j] == k && j ? I_3->f[i][j]*I_3->b[i][j-1]*psi[k] : NumType(0);
                curr_psi_cnt[k] += dnaSeq.ori[j] == k && j ? I_4->f[i][j]*I_4->b[i][j-1]*psi[k] : NumType(0);
                curr_psi_cnt[k] += dnaSeq.ori[j] == k && j ? I_5->f[i][j]*I_5->b[i][j-1]*psi[k] : NumType(0);
                curr_psi_cnt[k] += dnaSeq.ori[j] == k && j ? I_6->f[i][j]*I_6->b[i][j-1]*psi[k] : NumType(0);
                curr_psi_cnt[k] += dnaSeq.ori[j] == k && j ? I_7->f[i][j]*I_7->b[i][j-1]*psi[k] : NumType(0);
            } 
            for (int k = 0; k < 21; k ++) {
                curr_phi_cnt[k] += proSeq[i] == k && i ? D_1->f[i][j]*D_1->b[i-1][j]*phi[k] : NumType(0);
                curr_phi_cnt[k] += proSeq[i] == k && i ? D_2->f[i][j]*D_2->b[i-1][j]*phi[k] : NumType(0);
                curr_phi_cnt[k] += proSeq[i] == k && i ? D_3->f[i][j]*D_3->b[i-1][j]*phi[k] : NumType(0);
            }
            for (int s = 0; s < 21; s ++) {
                for (int t = 0; t < 64; t ++) {
                    curr_pi_cnt[s][t] += proSeq[i] == s && dnaSeq.triplet[j] == t && i && j > 2 ? 
                                        Match->f[i][j] * Match->b[i-1][j-3] * pi[s][t] : NumType(0);
                }
            }
        }
    }
    for (int k = 0; k < 4; k ++) {
        psi_cnt[k] += curr_psi_cnt[k]*reversep;
    }
    for (int k = 0; k < 21; k ++) {
        phi_cnt[k] += curr_phi_cnt[k]*reversep;
    }
    for (int s = 0; s < 21; s ++) {
        for (int t = 0; t < 64; t ++) {
            pi_cnt[s][t] += curr_pi_cnt[s][t]*reversep;
        }
    }
}

void PairHMM::logUpdateEmissions(const proSeqType & proSeq, const dnaSeqType & dnaSeq) {
    int n = start->logf.size() - 1;
    int m = start->logf[0].size() - 1; 
    //psi: insertion phi: deletion
    std::vector<std::vector<LogNumType> > curr_log_psi_cnt (4, std::vector<LogNumType> ());
    std::vector<std::vector<LogNumType> > curr_log_phi_cnt (21, std::vector<LogNumType> ());
    std::vector<std::vector<std::vector<LogNumType>>> curr_log_pi_cnt (21, std::vector<std::vector<LogNumType>> (64, std::vector<LogNumType> ()));
    for (int i = 0; i <= n; i ++) {
        for (int j = 1; j <= m; j ++) {
            int k = dnaSeq.ori[j];
            curr_log_psi_cnt[k].push_back(I_1->logf[i][j]+I_1->logb[i][j-1]);
        }
    }
    for (int i = 0; i <= n; i ++) {
        for (int j = 1; j <= m; j ++) {
            int k = dnaSeq.ori[j];
            curr_log_psi_cnt[k].push_back(I_2->logf[i][j]+I_2->logb[i][j-1]);
        }
    }
    for (int i = 0; i <= n; i ++) {
        for (int j = 1; j <= m; j ++) {
            int k = dnaSeq.ori[j];
            curr_log_psi_cnt[k].push_back(I_3->logf[i][j]+I_3->logb[i][j-1]);
        }
    }
    for (int i = 0; i <= n; i ++) {
        for (int j = 1; j <= m; j ++) {
            int k = dnaSeq.ori[j];
            curr_log_psi_cnt[k].push_back(I_4->logf[i][j]+I_4->logb[i][j-1]);
        }
    }
    for (int i = 0; i <= n; i ++) {
        for (int j = 1; j <= m; j ++) {
            int k = dnaSeq.ori[j];
            curr_log_psi_cnt[k].push_back(I_5->logf[i][j]+I_5->logb[i][j-1]);
        }
    }
    for (int i = 0; i <= n; i ++) {
        for (int j = 1; j <= m; j ++) {
            int k = dnaSeq.ori[j];
            curr_log_psi_cnt[k].push_back(I_6->logf[i][j]+I_6->logb[i][j-1]);
        }
    }
    for (int i = 0; i <= n; i ++) {
        for (int j = 1; j <= m; j ++) {
            int k = dnaSeq.ori[j];
            curr_log_psi_cnt[k].push_back(I_7->logf[i][j]+I_7->logb[i][j-1]);
        }
    }
    for (int k = 0; k < 4; k ++) psi_cnt[k] += exp(log_sum_exp(curr_log_psi_cnt[k].begin(), curr_log_psi_cnt[k].end())-logProb);
        //std::cout << k << std::endl;
    for (int i = 1; i <= n; i ++) {
        int k = proSeq[i];
        for (int j = 0; j <= m; j ++) {
            curr_log_phi_cnt[k].push_back(D_1->logf[i][j]+D_1->logb[i-1][j]);
        }
    }
    for (int i = 1; i <= n; i ++) {
        int k = proSeq[i];
        for (int j = 0; j <= m; j ++) {
            curr_log_phi_cnt[k].push_back(D_2->logf[i][j]+D_2->logb[i-1][j]);
        }
    }
    for (int i = 1; i <= n; i ++) {
        int k = proSeq[i];
        for (int j = 0; j <= m; j ++) {
            curr_log_phi_cnt[k].push_back(D_3->logf[i][j]+D_3->logb[i-1][j]);
        }
    }

    for (int k = 0; k < 21; k ++) phi_cnt[k] += exp(log_sum_exp(curr_log_phi_cnt[k].begin(), curr_log_phi_cnt[k].end())-logProb);

    for (int i = 1; i <= n; i ++) {
        for (int j = 3; j <= m; j ++) {
            int s = proSeq[i];
            int t = dnaSeq.triplet[j];
            curr_log_pi_cnt[s][t].push_back(Match->logf[i][j] + Match->logb[i-1][j-3]);
        }
    }
    //std::cout << "value check: " << Match->logf[1][3] << " " <<Match->logb[0][0] << " " << logProb << std::endl;
    for (int s = 0; s < 21; s ++) {
        for (int t = 0; t < 64; t ++) {
            double curr = exp(log_sum_exp(curr_log_pi_cnt[s][t].begin(), curr_log_pi_cnt[s][t].end())-log_pi[s][t]+log_psi[t&3]+log_psi[(t>>2)&3]+log_psi[(t>>4)&3]+log_phi[s]-logProb); 
            pi_cnt[s][t] += isnan(curr) ? 0 : curr;
            if(isnan(pi_cnt[s][t])) {
                std::cout << "NaN! " << s<<" " << t<<": "<<log_sum_exp(curr_log_pi_cnt[s][t].begin(), curr_log_pi_cnt[s][t].end()) <<" "
                <<log_pi[s][t]<<" "<<log_psi[t&3]<<" "<< log_psi[(t>>2)&3]<<" "<< log_psi[(t>>4)&3]<<" "<<log_phi[s]<< std::endl;
                exit(0);
            }
        }
    }
}

void PairHMM::updateEmissionProbabilities() {
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
    for (int i = 0; i < 21; i ++) {
        for (int j = 0; j < 64; j ++) {
            pi[i][j] = pi_cnt[i][j] / total_match;
        }
    }
}

void PairHMM::naiveUpdateProbabilities() {
    // transition part
    updateAlignProbabilities();
    naiveUpdateInsertionProbabilities();
    naiveUpdateDeletionProbabilities();
    // emmission part
    updateEmissionProbabilities();
    naiveTolog();
}

void PairHMM::naiveUpdateInsertionProbabilities() {
    beta_i = X_i.cnt / (X_i.cnt + B_i.cnt);
    epsilon_i = (X_i.cnt + B_i.cnt) / (X_i.cnt + B_i.cnt + E_i.cnt);
    delta_i = (X_i.cnt + B_i.cnt + E_i.cnt) / (X_i.cnt + B_i.cnt + E_i.cnt + D_i.cnt);
}

void PairHMM::naiveUpdateDeletionProbabilities() {
    beta_d = X_d.cnt / (X_d.cnt + B_d.cnt);
    epsilon_d = (X_d.cnt + B_d.cnt) / (X_d.cnt + B_d.cnt + E_d.cnt);
    delta_d = (X_d.cnt + B_d.cnt + E_d.cnt) / (X_d.cnt + B_d.cnt + E_d.cnt + D_d.cnt);
}

LogNumType PairHMM::calculateOverallLogProb(const std::vector<LogNumType> & parameters) const {
    LogNumType ln_gamma = parameters[0];
    LogNumType ln_alpha_d = parameters[1];
    LogNumType ln_alpha_i = parameters[2];
    LogNumType ln_beta_i = parameters[3];
    LogNumType ln_epsilon_i = parameters[4];
    LogNumType ln_delta_i = parameters[5];
    LogNumType ln_beta_d = parameters[6];
    LogNumType ln_epsilon_d = parameters[7];
    LogNumType ln_delta_d = parameters[8];
    return M.cnt*ln_gamma + A.cnt*Log1m(LogSumExp(ln_gamma, ln_alpha_d, ln_alpha_i)) + (B_d.cnt + D_d.cnt + E_d.cnt)*ln_alpha_d 
            + B_i.cnt*Log1m(ln_beta_i) + X_i.cnt*ln_beta_i + E_i.cnt*Log1m(ln_epsilon_i) + (X_i.cnt + B_i.cnt) * ln_epsilon_i 
            + D_i.cnt*(Log1m(ln_beta_i)) + (X_i.cnt + B_i.cnt + E_i.cnt)*ln_delta_i + B_d.cnt*(Log1m(ln_beta_d)) + X_d.cnt*ln_beta_d
            + E_d.cnt*Log1m(ln_epsilon_d) + (X_d.cnt + D_d.cnt)*ln_epsilon_d + D_d.cnt*Log1m(ln_epsilon_d) + (X_d.cnt + B_d.cnt + E_d.cnt)*ln_delta_d;
}

std::vector<LogNumType> PairHMM::deltaItoParameters(const LogNumType & deltai) const {
    /*
    LogNumType ln_gamma = parameters[0];
    LogNumType ln_alpha_d = parameters[1];
    LogNumType ln_alpha_i = parameters[2];
    LogNumType ln_beta_i = parameters[3];
    LogNumType ln_epsilon_i = parameters[4];
    LogNumType ln_delta_i = parameters[5];
    LogNumType ln_beta_d = parameters[6];
    LogNumType ln_epsilon_d = parameters[7];
    LogNumType ln_delta_d = parameters[8];
    */
    std::vector<LogNumType> parameters (9, 0.0);
    return parameters;
}

int PairHMM::insertionSolver() {
    double a = (2*D_i.cnt + 3*E_i.cnt + 3*X_i.cnt + 3*B_i.cnt);
    a *= a;
    double b = E_i.cnt + 3*X_i.cnt + B_i.cnt;
    double c = D_i.cnt + E_i.cnt + B_i.cnt;
    double d = 3*D_i.cnt + 4*E_i.cnt + 3*X_i.cnt + 3*B_i.cnt;
    double e = E_i.cnt + 3*X_i.cnt + 3*B_i.cnt;
    DComplex* roots = solve_quartic((-a*b - 3*a*c - d)/(a*c), (3*a*b + 3*a*c +3*e*d)/(a*c), (-3*a*b-3*e*e*d-a*c)/(a*c), (a*b + e*e*e)/(a*c));
    /*
    1. check real number in [0, 1]
    2. check other parameters in [0, 1]
    3. check objects
    4. check max
    5. update parameters
    6. what if failed => set error
    */
    std::vector<NumType> validDeltaI;
    std::vector<std::pair<NumType,NumType>> objects ;
    for (int i = 0; i < 4; i ++) {
        //std::cout << roots[i].real() << " + " << roots[i].imag() << "i" << std::endl;
        if(roots[i].imag() || roots[i].real() >= 1.0 || roots[i].real() <= 0) continue;
        //validDeltaI.push_back(roots[i].real());
        if(checkValidInsertionParameters(roots[i].real())){
            objects.push_back({deltaItoObject(roots[i].real()), roots[i].real()});
        }
    }
    delete roots;
    if(objects.size() == 0) {
        displayParameters("can not find optimal insertion parameters: ", error_filepath);
        return -1;
    }
    auto max_elem = max_element(objects.begin(), objects.end());
    setInsertionParameters(max_elem->second);
    return 0;
}

bool PairHMM::checkValidInsertionParameters(NumType DeltaI) const {
    // use deltaI to calculate beta_i, epsilon_i, check them in range [0, 1]
    if(DeltaI < 0 || DeltaI >= 1) return false;
    NumType EpsilonI = ((E_i.cnt + 3*X_i.cnt + 3*B_i.cnt))/((2*D_i.cnt + 3*E_i.cnt + 3*X_i.cnt + 3*B_i.cnt)*DeltaI) - (D_i.cnt + E_i.cnt)/(2*D_i.cnt + 3*E_i.cnt + 3*X_i.cnt + 3*B_i.cnt);
    NumType BetaI = pow(DeltaI, 2) * pow(1-EpsilonI, 3) / (EpsilonI*pow(1-DeltaI, 3));
    return EpsilonI > 0 && EpsilonI < 1 && BetaI > 0 && BetaI < 1;
}

NumType PairHMM::deltaItoObject(LogNumType DeltaI) {    
    NumType EpsilonI = ((E_i.cnt + 3*X_i.cnt + 3*B_i.cnt))/((2*D_i.cnt + 3*E_i.cnt + 3*X_i.cnt + 3*B_i.cnt)*DeltaI) - (D_i.cnt + E_i.cnt)/(2*D_i.cnt + 3*E_i.cnt + 3*X_i.cnt + 3*B_i.cnt);
    //std::cout << EpsilonI << " " << DeltaI << std::endl;
    return B_i.cnt * (log(EpsilonI*pow((1-DeltaI), 3)-DeltaI*DeltaI*pow(1-EpsilonI, 3))) + (E_i.cnt + 3*X_i.cnt + B_i.cnt)*log(DeltaI) 
        + (D_i.cnt - 3*X_i.cnt - 3*B_i.cnt)*log(1-DeltaI) + (E_i.cnt + 3*X_i.cnt)*log(1-EpsilonI);
}

void PairHMM::setInsertionParameters(NumType DeltaI) {
    NumType EpsilonI = ((E_i.cnt + 3*X_i.cnt + 3*B_i.cnt))/((2*D_i.cnt + 3*E_i.cnt + 3*X_i.cnt + 3*B_i.cnt)*DeltaI) - (D_i.cnt + E_i.cnt)/(2*D_i.cnt + 3*E_i.cnt + 3*X_i.cnt + 3*B_i.cnt);
    NumType BetaI = pow(DeltaI, 2) * pow(1-EpsilonI, 3) / (EpsilonI*pow(1-DeltaI, 3));
    setInsertion(BetaI, EpsilonI, DeltaI);
}

int PairHMM::deletionSolver() {
    double a = (2*D_d.cnt + 3*E_d.cnt + 3*X_d.cnt + 3*B_d.cnt);
    a *= a;
    double b = E_d.cnt + 3*X_d.cnt + B_d.cnt;
    double c = D_d.cnt + E_d.cnt + B_d.cnt;
    double d = 3*D_d.cnt + 4*E_d.cnt + 3*X_d.cnt + 3*B_d.cnt;
    double e = E_d.cnt + 3*X_d.cnt + 3*B_d.cnt;
    DComplex* roots = solve_quartic((-a*b - 3*a*c - d)/(a*c), (3*a*b + 3*a*c +3*e*d)/(a*c), (-3*a*b-3*e*e*d-a*c)/(a*c), (a*b + e*e*e)/(a*c));
    /*
    1. check real number in [0, 1]
    2. check other parameters in [0, 1]
    3. check objects
    4. check max
    5. update parameters
    6. what if failed => set error
    */
    std::vector<NumType> validDeltaD;
    std::vector<std::pair<NumType,NumType>> objects ;
    for (int i = 0; i < 4; i ++) {
        //std::cout << roots[i].real() << " + " << roots[i].imag() << "i" << std::endl;
        if(roots[i].imag() || roots[i].real() >= 1.0 || roots[i].real() <= 0) continue;
        //validDeltaD.push_back(roots[i].real());
        if(checkValidDeletionParameters(roots[i].real())){
            objects.push_back({deltaDtoObject(roots[i].real()), roots[i].real()});
        }
    }
    delete roots;
    if(objects.size() == 0) {
        displayParameters("can not find optimal deletion parameters: ", error_filepath);
        return -1;
    }
    auto max_elem = max_element(objects.begin(), objects.end());
    setDeletionParameters(max_elem->second);
    return 0;
}

bool PairHMM::checkValidDeletionParameters(NumType DeltaD) const {
    // use deltaD to calculate beta_d, epsilon_d, check them in range [0, 1]
    if(DeltaD < 0 || DeltaD >= 1) return false;
    NumType EpsilonD = ((E_d.cnt + 3*X_d.cnt + 3*B_d.cnt))/((2*D_d.cnt + 3*E_d.cnt + 3*X_d.cnt + 3*B_d.cnt)*DeltaD) - (D_d.cnt + E_d.cnt)/(2*D_d.cnt + 3*E_d.cnt + 3*X_d.cnt + 3*B_d.cnt);
    NumType BetaD = pow(DeltaD, 2) * pow(1-EpsilonD, 3) / (EpsilonD*pow(1-DeltaD, 3));
    return EpsilonD >= 0 && EpsilonD <= 1 && BetaD >= 0 && BetaD <= 1;
}

NumType PairHMM::deltaDtoObject(LogNumType DeltaD) {    
    NumType EpsilonD = ((E_d.cnt + 3*X_d.cnt + 3*B_d.cnt))/((2*D_d.cnt + 3*E_d.cnt + 3*X_d.cnt + 3*B_d.cnt)*DeltaD) - (D_d.cnt + E_d.cnt)/(2*D_d.cnt + 3*E_d.cnt + 3*X_d.cnt + 3*B_d.cnt);
    //std::cout << EpsilonD << " " << DeltaD << std::endl;
    return B_d.cnt * (log(EpsilonD*pow((1-DeltaD), 3)-DeltaD*DeltaD*pow(1-EpsilonD, 3))) + (E_d.cnt + 3*X_d.cnt + B_d.cnt)*log(DeltaD) 
        + (D_d.cnt - 3*X_d.cnt - 3*B_d.cnt)*log(1-DeltaD) + (E_d.cnt + 3*X_d.cnt)*log(1-EpsilonD);
}

void PairHMM::setDeletionParameters(NumType DeltaD) {
    NumType EpsilonD = ((E_d.cnt + 3*X_d.cnt + 3*B_d.cnt))/((2*D_d.cnt + 3*E_d.cnt + 3*X_d.cnt + 3*B_d.cnt)*DeltaD) - (D_d.cnt + E_d.cnt)/(2*D_d.cnt + 3*E_d.cnt + 3*X_d.cnt + 3*B_d.cnt);
    NumType BetaD = pow(DeltaD, 2) * pow(1-EpsilonD, 3) / (EpsilonD*pow(1-DeltaD, 3));
    setDeletion(BetaD, EpsilonD, DeltaD);
}

void PairHMM::setPsi(std::vector<NumType> & Psi) {
    for (int i = 0; i < Psi.size(); i ++) {
        psi[i] = Psi[i];
    }
}

void PairHMM::setPhi(std::vector<NumType> & Phi) {
    for (int i = 0; i < Phi.size(); i ++) {
        phi[i] = Phi[i];
    }
}

void PairHMM::setPi(std::vector<std::vector<NumType>> mat) {
    pi = mat;
    for (int i = 0; i < pi.size(); i ++) {
        for (int j = 0; j < pi[0].size(); j ++) {
            log_pi[i][j] = log(pi[i][j]);
        }
    }
}

void PairHMM::setInsertion(NumType betaI, NumType epsilonI, NumType deltaI) {
    beta_i = betaI;
    epsilon_i = epsilonI;
    delta_i = deltaI;
} 

void PairHMM::setDeletion(NumType betaD, NumType epsilonD, NumType deltaD) {
    beta_d = betaD;
    epsilon_d = epsilonD;
    delta_d = deltaD;
}

void PairHMM::setAlign(NumType omegaI, NumType omegaD, NumType Gamma, NumType alphaI, NumType alphaD) {
    omega_i = omegaI;
    omega_d = omegaD;
    gamma = Gamma;
    alpha_i = alphaI;
    alpha_d = alphaD;
}

void PairHMM::updateProbabilities(int option){
    switch (option)
    {
    case 0:
        naiveUpdateProbabilities();
        break;
    
    case 1:
        optimizedUpdateProbabilities();
        break;
    }
}

void PairHMM::optimizedUpdateProbabilities(){
    updateAlignProbabilities();
    if (insertionSolver() == -1) naiveUpdateInsertionProbabilities();
    if (deletionSolver() == -1) naiveUpdateDeletionProbabilities();
    updateEmissionProbabilities();
    naiveTolog();
}

void PairHMM::updateAlignProbabilities() {
    NumType temp = A.cnt + M.cnt + B_d.cnt + D_d.cnt + E_d.cnt + B_i.cnt + D_i.cnt + E_i.cnt;
    gamma = M.cnt / temp;
    alpha_i = (B_i.cnt + D_i.cnt + E_i.cnt) / temp;
    alpha_d = (B_d.cnt + D_d.cnt + E_d.cnt) / temp;
    omega_d = (J_d.cnt + K_d.cnt) / (J_d.cnt + K_d.cnt + 2*A.cnt);
    omega_i = (J_i.cnt + K_i.cnt) / (J_i.cnt + K_i.cnt + 2*A.cnt);
}

void PairHMM::displayParameters(const std::string msg, const std::string & filename) const{
    std::ofstream out;
    out.open(filename, std::ios::out|std::ios::app);
    std::vector<NumType> vt {
        J_d.cnt, J_i.cnt, M.cnt, A.cnt, K_d.cnt, K_i.cnt, 
        F_d.cnt, X_d.cnt, B_d.cnt, D_d.cnt, E_d.cnt, G_d.cnt, H_d.cnt,
        B_i.cnt, E_i.cnt, D_i.cnt, F_i.cnt, G_i.cnt, H_i.cnt, X_i.cnt,
    };
    std::vector<NumType> tp {
        omega_i, omega_d, gamma, alpha_i, alpha_d, gamma, delta_i, beta_i, epsilon_i, delta_d, beta_d, epsilon_d,
    };
    if(msg.size()) {
        out << msg << std::endl;
    }
    /*for (int i = 0; i < vt.size(); i ++) {
        if(i) out << '\t';
        out << vt[i];
    }
    out << std::endl;*/
    out << "cnts: ";
    for (auto num: vt) out << num << "\t";
    out << std::endl;
    out << logProb << std::endl;
    for (int i = 0; i < tp.size(); i ++) {
        if(i) out << '\t';
        out << tp[i];
    }
    out << std::endl;
    /*
    for (int i = 0; i < phi_cnt.size(); i ++) {
        //if(i) out << '\t';
        out <<'\t' << phi_cnt[i];
    }
    out << std::endl;
    for (int i = 0; i < psi_cnt.size(); i ++) {
        //if(i) out << '\t';
        out << '\t'<<psi_cnt[i];
    }
    */
    for (int i = 0; i < phi.size(); i ++) out << phi[i] << '\t';
    out << std::endl;
    for (int i = 0; i < psi.size(); i ++) out << psi[i] << '\t';
    out << std::endl;
    for (int i = 0; i < pi_cnt.size(); i ++) {
        for (int j = 0; j < pi_cnt[0].size(); j ++) {
            out << pi[i][j] << '\t';
        }
    } 
    out << std::endl;
}

void PairHMM::displayEmissionCnts() {
    std::cout << "Emission counts:" << std::endl;
    std::cout <<  "T: " << psi_cnt[0] << " C: " << psi_cnt[1] << ", A: " << psi_cnt[2] << ", G: " << psi_cnt[3] << std::endl;
    DataTool dt;
    for (int i = 0; i < 21; i ++) {
        if(i) std::cout << " ";
        std::cout << dt.decodeAA(i) << " " << phi_cnt[i] ;
    }
    std::cout << std::endl;
    for (int i = 0; i < 21; i ++) {
        for (int j = 0; j < 64; j ++) {
            if(j) std::cout << " ";
            std::cout << dt.decodeAA(i) << " " << dt.decodeTriplet(j) << ": " <<pi_cnt[i][j]; 
        }
        std::cout << std::endl;
    }
}

void PairHMM::displayTransitionCnts() {
    std::cout << "Transition counts" << std::endl;
    std::vector<Transition*> vt {
        &J_d, &J_i, &M, &A, &K_d, &K_i, 
        &F_d, &X_d, &B_d, &D_d, &E_d, &G_d, &H_d,
        &B_i, &E_i, &D_i, &F_i, &G_i, &H_i, &X_i
    };
    std::cout << "J_d: " << J_d.cnt << std::endl; 
    std::cout << "J_i: " << J_i.cnt << std::endl; 
    std::cout << "K_d: " << K_d.cnt << std::endl; 
    std::cout << "K_i: " << K_i.cnt << std::endl; 
    std::cout << "M: " << M.cnt << std::endl;
    std::cout << "A: " << A.cnt << std::endl;  
    std::cout << "F_d: " << F_d.cnt << std::endl; 
    std::cout << "X_d: " << X_d.cnt << std::endl;
    std::cout << "B_d: " << B_d.cnt << std::endl;
    std::cout << "D_d: " << D_d.cnt << std::endl;
    std::cout << "E_d: " << E_d.cnt << std::endl;
    std::cout << "G_d: " << G_d.cnt << std::endl;
    std::cout << "H_d: " << H_d.cnt << std::endl;
    std::cout << "F_i: " << F_i.cnt << std::endl; 
    std::cout << "X_i: " << X_i.cnt << std::endl;
    std::cout << "B_i: " << B_i.cnt << std::endl;
    std::cout << "D_i: " << D_i.cnt << std::endl;
    std::cout << "E_i: " << E_i.cnt << std::endl;
    std::cout << "G_i: " << G_i.cnt << std::endl;
    std::cout << "H_i: " << H_i.cnt << std::endl;
}

void PairHMM::naiveTolog() {
    //omega_d, omega_i, gamma, alpha_d, beta_d, epsilon_d, delta_d, alpha_i, beta_i, epsilon_i, delta_i
    //log_omega_d, log_omega_i, log_gamma, log_alpha_d, log_beta_d, log_epsilon_d, log_delta_d, log_alpha_i, log_beta_i, log_epsilon_i, log_delta_i
    log_omega_d = log(omega_d);
    log_omega_i = log(omega_i);
    log_gamma = log(gamma);
    log_alpha_d = log(alpha_d);
    log_beta_d = log(beta_d);
    log_epsilon_d = log(epsilon_d);
    log_delta_d = log(delta_d);
    log_alpha_i = log(alpha_i);
    log_beta_i = log(beta_i);
    log_epsilon_i = log(epsilon_i);
    log_delta_i = log(delta_i);
    for (int i = 0; i < phi.size(); i ++) log_phi[i] = log(phi[i]);
    for (int i = 0; i < psi.size(); i ++) log_psi[i] = log(psi[i]);
    for (int i = 0; i < pi.size(); i ++) {
        for (int j = 0; j < pi[0].size(); j ++) {
            log_pi[i][j] = log(pi[i][j]);
        }
    }
}

void PairHMM::logToNaive() {
    omega_d = exp(log_omega_d);
    omega_i = exp(log_omega_i);
    gamma = exp(log_gamma);
    alpha_d = exp(log_alpha_d);
    beta_d = exp(log_beta_d);
    epsilon_d = exp(log_epsilon_d);
    delta_d = exp(log_delta_d);
    alpha_i = exp(log_alpha_i);
    beta_i = exp(log_beta_i);
    epsilon_i = exp(log_epsilon_i);
    delta_i = exp(log_delta_i);
    for (int i = 0; i < phi.size(); i ++) phi[i] = exp(log_phi[i]);
    for (int i = 0; i < psi.size(); i ++) psi[i] = exp(log_psi[i]);
    for (int i = 0; i < pi.size(); i ++) {
        for (int j = 0; j <  pi[0].size(); j ++) {
            pi[i][j] = exp(log_pi[i][j]);
        }
    }
}

void PairHMM::checkback(int i, int j, int n, int m) {
    double suffix = (n-i) *log_phi[0] + (m-j)*log_psi[0];
    std::vector<State*> veci {I_1, I_2, I_3, I_4, I_5, I_6, I_7};
    std::vector<State*> vech {H_1, H_2, H_3, H_4, H_5, H_6, H_7};
    std::vector<State*> vecd {D_1, D_2, D_3, Match};
    std::cout << "Backward check: "<<i<<" & "<< j << ": " << std::endl;
    for(auto & ptr: veci) {
        std::cout << "ori: " << ptr->b[i][j] <<" log: "<<ptr->logb[i][j]<<", error: " <<exp(ptr->logb[i][j] + suffix) - ptr->b[i][j] << std::endl; 
    }
    std::cout << std::endl;
    for(auto & ptr: vech) {
        std::cout << "ori: " << ptr->b[i][j] <<" log: "<<ptr->logb[i][j]<<", error: " <<exp(ptr->logb[i][j] + suffix) - ptr->b[i][j] << std::endl; 
    }
    std::cout << std::endl;
    for(auto & ptr: vecd) {
        std::cout << "ori: " << ptr->b[i][j] <<" log: "<<ptr->logb[i][j]<<", error: " <<exp(ptr->logb[i][j] + suffix) - ptr->b[i][j] << std::endl; 
    }
    std::cout << std::endl;
}

void PairHMM::checkforward(int i, int j) {
    double prefix = (i) *log_phi[0] + (j)*log_psi[0];
    std::vector<State*> veci {I_1, I_2, I_3, I_4, I_5, I_6, I_7};
    std::vector<State*> vech {H_1, H_2, H_3, H_4, H_5, H_6, H_7};
    std::vector<State*> vecd {D_1, D_2, D_3, Match};
    std::cout << "Forward check: "<<i<<" & "<< j << ": " << std::endl;
    for(auto & ptr: veci) {
        std::cout << "ori: " << ptr->f[i][j] <<" log: "<<ptr->logf[i][j]<<", error: " <<exp(ptr->logf[i][j] + prefix) - ptr->f[i][j] << std::endl; 
    }
    std::cout << std::endl;
    for(auto & ptr: vech) {
        std::cout << "ori: " << ptr->f[i][j] <<" log: "<<ptr->logf[i][j]<<", error: " <<exp(ptr->logf[i][j] + prefix) - ptr->f[i][j] << std::endl; 
    }
    std::cout << std::endl;
    for(auto & ptr: vecd) {
        std::cout << "ori: " << ptr->f[i][j] <<" log: "<<ptr->logf[i][j]<<", error: " <<exp(ptr->logf[i][j] + prefix) - ptr->f[i][j] << std::endl; 
    }
    std::cout << std::endl;
}

void PairHMM::checkEmissions() {
    DataTool dt;
    for(int i = 0; i < phi_cnt.size(); i ++) {
        double curr = phi_cnt[i];
        for(int j = 0; j < pi_cnt[0].size(); j ++) {
            curr += pi_cnt[i][j];
        }
        std::cout << dt.decodeAA(i) << ": " << curr << " ";
    }
    std::cout << std::endl;

    for(int i = 0; i < 4; i ++) {
        double curr = psi_cnt[i];
        for(int s = 0; s < 21; s ++) {
            for (int j = 0; j < 64; j ++) {
                if((j & 3) == i) curr += pi_cnt[s][j];
                if(((j >> 2) & 3) == i) curr += pi_cnt[s][j];
                if(((j >> 4) & 3) == i) curr += pi_cnt[s][j];
            }
        }
        std::cout << dt.decodeBase(i) << ": " << curr << " ";
    }
    std::cout << std::endl;
}

void PairHMM::checkTransitionParameters() {
    // equal frameshift cost
    std::cout << beta_i*epsilon_i*pow(1.0-delta_i, 3) - pow(delta_i, 2)*pow(1-epsilon_i, 3) << std::endl;
    std::cout << beta_d*epsilon_d*pow(1.0-delta_d, 3) - pow(delta_d, 2)*pow(1-epsilon_d, 3) << std::endl;
    std::cout << (2*D_i.cnt + 3*E_i.cnt + 3*X_i.cnt + 3*B_i.cnt)*delta_i*epsilon_i + (D_i.cnt + E_i.cnt)*delta_i - (E_i.cnt + 3*X_i.cnt + 3*B_i.cnt) << std::endl;
    std::cout << (2*D_d.cnt + 3*E_d.cnt + 3*X_d.cnt + 3*B_d.cnt)*delta_d*epsilon_d + (D_d.cnt + E_d.cnt)*delta_d - (E_d.cnt + 3*X_d.cnt + 3*B_d.cnt) << std::endl;
}

void PairHMM::setInsertion(NumType deltaI, NumType epsilonI) {
    NumType BetaI = pow(deltaI, 2)*pow(1-epsilonI, 3) / (pow(1-deltaI, 3)*epsilonI);
    setInsertion(BetaI, deltaI, epsilonI);
}

void PairHMM::setDeletion(NumType deltaD, NumType epsilonD) {
    NumType BetaD = pow(deltaD, 2)*pow(1-epsilonD, 3) / (pow(1-deltaD, 3)*epsilonD);
    setDeletion(BetaD, deltaD, epsilonD);
}

void PairHMM::testTraining(const std::string & filename) {
    std::ifstream in(filename, std::ios::in);
    std::string dna_position, dna_seq, protein_id, protein_seq;
    std::vector<proSeqType> pros;
    std::vector<dnaSeqType> dnas;
    DataTool dt; 
    while(in >> protein_id >> protein_seq >> dna_position >> dna_seq) {
        if(dt.checkDNA(dna_seq) && dt.checkPro(protein_seq)) {
            pros.push_back(dt.encodePro(protein_seq+"*"));
            dnas.push_back(dt.encodeDNA(dna_seq));
        }
    }
    in.close();
    std::cout << "loaded data: " << pros.size()<< " pro seqs and "<<dnas.size()<<" dna seqs."<<std::endl;
    naiveBaumWelch(pros, dnas, 4, 1);
}

void PairHMM::get_total() {
    double deletionCnts = J_d.cnt + K_d.cnt + F_d.cnt + X_d.cnt + M.cnt;
    double insertionCnts = J_i.cnt + K_i.cnt + F_i.cnt + G_i.cnt + H_i.cnt + X_i.cnt +2*D_d.cnt + E_d.cnt + 3*M.cnt;
    double proCnts = 0;
    double dnaCnts = 0;
    for (int i = 0; i < 4; i ++) dnaCnts += psi_cnt[i];
    for (int i = 0; i < 21; i ++) proCnts += phi_cnt[i];
    double tot = 0;
    for (int i = 0; i < 21; i ++) {
        for (int j = 0; j < 64; j ++) tot += pi_cnt[i][j];
    }
    dnaCnts += 3*tot;
    proCnts += tot;
    std::cout << "total deletions: " << deletionCnts << " " << proCnts << std::endl;
    std::cout << "total insertions: " << insertionCnts << " " << dnaCnts <<std::endl;
}

void PairHMM::pseudocount(int n) {
    for (int i = 0; i < pi_cnt.size(); i ++) {
        for (int j = 0; j < pi_cnt.size(); j ++) {
            pi_cnt[i][j] += n;
        }
    }
    for (int i = 0; i < phi_cnt.size(); i ++) {phi_cnt[i] += n; }
    for (int i = 0; i < psi_cnt.size(); i ++) {psi_cnt[i] += n; }
    std::vector<Transition*> vt {
        &J_d, &J_i, &M, &A, &K_d, &K_i, 
        &F_d, &X_d, &B_d, &D_d, &E_d, &G_d, &H_d,
        &F_i, &X_i, &B_i, &D_i, &E_i, &G_i, &H_i
    };
    for (auto & ptr: vt) {
        ptr->cnt += n;
    }

}