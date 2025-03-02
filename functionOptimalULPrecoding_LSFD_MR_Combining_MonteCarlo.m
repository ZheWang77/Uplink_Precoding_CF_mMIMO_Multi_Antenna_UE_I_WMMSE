function [SE_sum_Iterative_Monte_MR,SE_Optimal_Monte_MR,F_Precoding_Optimal_Iterative_Monte_MR,GG_Monte_Iterative] = functionOptimalULPrecoding_LSFD_MR_Combining_MonteCarlo(Hhat_MMSE,HH,F_precoding_Initial,W,M,K,L,N,p,totIter,nbrOfRealizations,tau_p,tau_c,Delta)
warning('off');
%%=============================================================
%This function is used to generate the I-WMMSE precoding matrices and compute the achievable SE for the LSFD processing scheme with MR combining of the paper:
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
SE_sum_Iterative_Monte_MR = zeros(totIter,1);
SE_Monte_total = zeros(K,totIter);
F_Precoding_Optimal_Iterative_Monte_MR = zeros(N,N,K,totIter);
GG_Monte_Iterative = zeros(N,N,K,K,totIter);
maxIter = totIter;
iter = 0;

SE_sum_objective_monte = 0.00001;
F_precoding = F_precoding_Initial;

%% ==========Level3  Monte-Carlo  MR Combining====================%%

%-----Monte-Carlo

while(iter<maxIter)
    
    iter = iter+1;
    
    SE_sum_objective_old_monte = SE_sum_objective_monte;

    [A_Monte,E_Monte,signal_MR_Monte,SE_Monte] = functionCompute_SE_LSFD_Monte_Carlo(Hhat_MMSE,HH,F_precoding,tau_c,tau_p,nbrOfRealizations,N,L,K,M);
    [GG_Monte] = functionMatrixGeneration_LSFD_Monte_Carlo(Hhat_MMSE,HH,A_Monte,E_Monte,W,nbrOfRealizations,N,L,K,M);
    
    %---Monte-Carlo
    [F_precoding] = functionOptimalObjective_LSFD(A_Monte,E_Monte,GG_Monte,signal_MR_Monte,W,K,N,p);
    F_Precoding_Optimal_Iterative_Monte_MR(:,:,:,iter) = F_precoding;
    
    %--Monte-Carlo
    SE_sum_objective_monte = sum(SE_Monte.*W);
    SE_sum_Iterative_Monte_MR(iter) = SE_sum_objective_monte;
    GG_Monte_Iterative(:,:,:,:,iter) = GG_Monte;
    SE_Monte_total(:,iter) = SE_Monte;
    
        if  SE_sum_objective_monte - SE_sum_objective_old_monte<0

        SE_Monte = SE_Monte_total(:,iter-1);
        SE_sum_Iterative_Monte_MR(iter) = 0; 
        break; 
    end

    if  abs(SE_sum_objective_monte-SE_sum_objective_old_monte)/SE_sum_objective_old_monte<=Delta
        
        break;  
    end
    

    
end
SE_Optimal_Monte_MR = SE_Monte;
