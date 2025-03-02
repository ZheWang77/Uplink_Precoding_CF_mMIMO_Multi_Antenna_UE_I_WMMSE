function [SE_sum_Iterative_Monte_LMMSE,SE_Optimal_Monte_LMMSE,F_Precoding_Optimal_Iterative_Monte_LMMSE] = functionOptimalULPrecoding_LSFD_LMMSE_Combining_MonteCarlo(Hhat_MMSE,HH,F_precoding_Initial,C_MMSE,W,M,K,L,N,p,totIter,nbrOfRealizations,tau_p,tau_c,Delta)
warning('off');
%%=============================================================
%This function is used to generate the I-WMMSE precoding matrices and compute the achievable SE for the LSFD processing scheme with L-MMSE combining of the paper:
%
% Z. Wang, J. Zhang, H. Q. Ngo, B. Ai, and M. Debbah, "Uplink Precoding Design for Cell-Free Massive MIMO With Iteratively Weighted MMSE," 
% in IEEE Transactions on Communications, vol. 71, no. 3, pp. 1646-1664, March 2023, doi: 10.1109/TCOMM.2023.3235919.

%
%Download article: https://arxiv.org/abs/2301.02417 or https://ieeexplore.ieee.org/document/10013728
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.
%%=============================================================


%---Initial Value
SE_sum_Iterative_Monte_LMMSE = zeros(totIter,1);
SE_Monte_total = zeros(K,totIter);
F_Precoding_Optimal_Iterative_Monte_LMMSE = zeros(N,N,K,totIter);

% Delta = 1e-4; %Threshold
maxIter = totIter;
iter = 0;

SE_sum_objective_monte = 0.00001;
F_precoding = F_precoding_Initial;


%% ==========Level3  Monte-Carlo  L-MMSE Combining====================%%

%-----Monte-Carlo
while(iter<maxIter)
    
    iter = iter+1;
    
    SE_sum_objective_old_monte = SE_sum_objective_monte;
    
    [V_MMSE_Combining] = functionCompute_MMSE_Combining_Matrix(Hhat_MMSE,F_precoding,C_MMSE,nbrOfRealizations,L,N,K,M);
    [A_Monte,E_Monte,signal_MR_Monte,SE_Monte] = functionCompute_SE_LSFD_Monte_Carlo(V_MMSE_Combining,HH,F_precoding,tau_c,tau_p,nbrOfRealizations,N,L,K,M);
    [GG_Monte] = functionMatrixGeneration_LSFD_Monte_Carlo(V_MMSE_Combining,HH,A_Monte,E_Monte,W,nbrOfRealizations,N,L,K,M);

    %---Monte-Carlo Optimal Precoding
    [F_precoding] = functionOptimalObjective_LSFD(A_Monte,E_Monte,GG_Monte,signal_MR_Monte,W,K,N,p);
    F_Precoding_Optimal_Iterative_Monte_LMMSE(:,:,:,iter) = F_precoding;
    
    %--Monte-Carlo
    SE_sum_objective_monte = sum(SE_Monte.*W);
    SE_sum_Iterative_Monte_LMMSE(iter) = SE_sum_objective_monte;
    SE_Monte_total(:,iter) = SE_Monte;
    
    
    if  SE_sum_objective_monte - SE_sum_objective_old_monte<0

        SE_Monte = SE_Monte_total(:,iter-1);
        SE_sum_Iterative_Monte_LMMSE(iter) = 0; 
        break; 
    end
    
    if  abs(SE_sum_objective_monte-SE_sum_objective_old_monte)/SE_sum_objective_old_monte<=Delta
        
        
        break;  
    end
    
end
SE_Optimal_Monte_LMMSE = SE_Monte;