function [SE_sum_Iterative_Analy,SE_Optimal_Analytical,F_Precoding_Optimal_Iterative_Analytical,GG_Iterative_Analytical] = functionOptimalULPrecoding_LSFD_MR_Combining_Analytical(Gkk,S_k,X_p1_all,X_p2_all_1,X_p2_all_2,X_p3_1,X_p3_2,F_precoding_Initial,GG_total,W,M,K,N,p,totIter,tau_p,tau_c,Pset,Delta)
warning('off');
%%=============================================================
%This function is used to generate the closed-form I-WMMSE precoding matrices and compute the closed-form achievable SE for the LSFD processing scheme with MR combining of the paper:
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
F_precoding = F_precoding_Initial; %--Initial UL Precoding Structure
SE_sum_Iterative_Analy = zeros(totIter,1);
SE_Closed_total = zeros(K,totIter);
F_Precoding_Optimal_Iterative_Analytical = zeros(N,N,K,totIter);
GG_Iterative_Analytical = zeros(N,N,K,K,totIter);

maxIter = totIter;
iter = 0;
SE_sum_objective_analy = 0.00001;


%% ==========Level3  Analytical  MR Combining====================%%
while(iter<maxIter)
    
    
    iter = iter+1;
    SE_sum_objective_old_analy = SE_sum_objective_analy;
    
    %---Analytical
    [A,E,SE_analy] = functionCompute_SE_LSFD_Analytical_For_Iteration(Gkk,S_k,X_p1_all,X_p2_all_1,X_p2_all_2,X_p3_1,X_p3_2,F_precoding,M,K,N,tau_p,tau_c,Pset);
    
    [GG] = functionMatrixGeneration_LSFD_Analytical(A,E,GG_total,W,K,N);
    
    [F_precoding] = functionOptimalObjective_LSFD(A,E,GG,Gkk,W,K,N,p);
    F_Precoding_Optimal_Iterative_Analytical(:,:,:,iter) = F_precoding;
    
    
   %--Analytical
    SE_sum_objective_analy = sum(SE_analy.*W);
    
    SE_sum_Iterative_Analy(iter) = SE_sum_objective_analy;
    GG_Iterative_Analytical(:,:,:,:,iter) = GG;
    SE_Closed_total(:,iter) = SE_analy;
    
    if  abs(SE_sum_objective_analy-SE_sum_objective_old_analy)/SE_sum_objective_old_analy<=Delta
        
        
        break;
        
    end
    
    if  SE_sum_objective_analy - SE_sum_objective_old_analy<0

        SE_analy = SE_Closed_total(:,iter-1);
        SE_sum_Iterative_Analy(iter) = 0; 
        break; 
        
    end
    
    
end
SE_Optimal_Analytical = SE_analy;
