

clear all
close all
tic

M = 10;
K = 10;
L = 2;
N = 2;
tau_c = 200;
nbrOfSetups = 1;
nbrOfRealizations = 800;
w = 2; %Pilot Reuse factor
p = 0.2; %200 mW
%Create the power vector for all UEs (The uplink power is the same
%(p)at each UE)
p_u = p*ones(1,K);
totIter = 10;
Delta = 5e-4; %Threshold


W = ones(K,1);


SE_Optimal_Monte_MMSE_Total = zeros(K,nbrOfSetups,length(N));
SE_sum_Iterative_Monte_MMSE_Level4_Total = zeros(totIter,nbrOfRealizations,nbrOfSetups,length(N));

SE_sum_Iterative_Monte_LMMSE_Level3_Total = zeros(totIter,nbrOfSetups,length(N));
SE_sum_Iterative_Monte_MR_Level3_Total = zeros(totIter,nbrOfSetups,length(N));
SE_sum_Iterative_Theoretical_MR_Level3_Total = zeros(totIter,nbrOfSetups,length(N));

SE_Theoretical_MR_Combining_Level3_Iterative = zeros(K,nbrOfSetups,length(N));
SE_Monte_MMSE_Combining_Level3_Iterative = zeros(K,nbrOfSetups,length(N));
SE_Monte_MR_Combining_Level3_Iterative = zeros(K,nbrOfSetups,length(N));

n = 1;


for i = 1:nbrOfSetups
    
    [channelGain] = RandomAP_generateSetup_Rician_Multi_Antenna(M,K,1,1);
    
        
        F_pre_wo = zeros(N(n),N(n),K);
        
        for k = 1:K
            
            F_pre_wo(:,:,k) = sqrt(p)*1/sqrt(N(n))*eye(N(n),N(n));
            
        end
        
        tau_p = K*N(n)/w;
        
        %------Weichselberger Model
        [HH,R_Vec,H_Vec,Omega_couple,F_precoding_0] = functionChannelGeneration(channelGain,M,K,N(n),L,p_u,nbrOfRealizations);
        
        [Pset] = functionPilotAllocation( R_Vec,M,K,L*N(n),tau_p/N(n),p_u);
        % disp(['PilotAllocation of setup ' num2str(i)]);
        
        [Hhat_MMSE] = functionChannelEstimates_MMSE(H_Vec,R_Vec,F_pre_wo,nbrOfRealizations,M,K,L,N(n),tau_p,Pset);
       
        [Phi,Omega,C_MMSE] = functionMatrixGeneration(R_Vec,F_pre_wo,M,K,L,N(n),tau_p,Pset);
        
        [SE_sum_Iterative_Monte_MMSE,SE_Optimal_Monte_MMSE,SE_Optimal_Monte_MMSE_total,F_Precoding_Optimal_Iterative_Monte_MMSE] = functionOptimalULPrecoding_Fully_Centralized_Level4(Hhat_MMSE,F_pre_wo,C_MMSE,W,M,K,L,N(n),p_u,totIter,nbrOfRealizations,tau_p,tau_c);
        
        SE_Optimal_Monte_MMSE_Total(:,i,n) = SE_Optimal_Monte_MMSE;
        SE_sum_Iterative_Monte_MMSE_Level4_Total(:,:,i,n) = SE_sum_Iterative_Monte_MMSE;
        
        %--Iterative Optimization
        %--MR Monte-Carlo
        
        [SE_sum_Iterative_Monte_MR,SE_Optimal_Monte_MR,F_Precoding_Optimal_Iterative_Monte_MR,GG_Monte_Iterative] = functionOptimalULPrecoding_LSFD_MR_Combining_MonteCarlo(Hhat_MMSE,HH,F_pre_wo,W,M,K,L,N(n),p_u,totIter,nbrOfRealizations,tau_p,tau_c,Delta);
        
        [Gkk,S_k,X_p1_all,X_p2_all_1,X_p2_all_2,X_p3_1,X_p3_2,GG_total] = functionAnalyticalMatrix_Iteration( R_Vec,F_pre_wo,Phi,Omega,M,K,L,N(n),tau_p,Pset);
        [SE_sum_Iterative_Analy,SE_Optimal_Analytical,F_Precoding_Optimal_Iterative_Analytical,GG_Iterative_Analytical] = functionOptimalULPrecoding_LSFD_MR_Combining_Analytical(Gkk,S_k,X_p1_all,X_p2_all_1,X_p2_all_2,X_p3_1,X_p3_2,F_pre_wo,GG_total,W,M,K,N(n),p_u,totIter,tau_p,tau_c,Pset,Delta);
        
        
        % --L-MMSE Monte-Carlo
        [SE_sum_Iterative_Monte_LMMSE,SE_Optimal_Monte_LMMSE,F_Precoding_Optimal_Iterative_Monte_LMMSE] = functionOptimalULPrecoding_LSFD_LMMSE_Combining_MonteCarlo(Hhat_MMSE,HH,F_pre_wo,C_MMSE,W,M,K,L,N(n),p_u,totIter,nbrOfRealizations,tau_p,tau_c,Delta);

        SE_sum_Iterative_Monte_LMMSE_Level3_Total(:,i,n) = SE_sum_Iterative_Monte_LMMSE;
        SE_sum_Iterative_Monte_MR_Level3_Total(:,i,n)  = SE_sum_Iterative_Monte_MR;
        SE_sum_Iterative_Theoretical_MR_Level3_Total(:,i,n) = SE_sum_Iterative_Analy;
        
        SE_Monte_MMSE_Combining_Level3_Iterative(:,i,n) = SE_Optimal_Monte_LMMSE;
        SE_Monte_MR_Combining_Level3_Iterative(:,i,n) = SE_Optimal_Monte_MR;
        SE_Theoretical_MR_Combining_Level3_Iterative(:,i,n) = SE_Optimal_Analytical;
             
        disp([num2str(i) ' setups out of ' num2str(nbrOfSetups)]);
     

end
        
toc  
