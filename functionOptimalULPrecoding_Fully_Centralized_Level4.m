function [SE_sum_Iterative_Monte_MMSE,SE_Optimal_Monte_MMSE_average,SE_Optimal_Monte_MMSE,F_Precoding_Optimal_Iterative_Monte_MMSE] = functionOptimalULPrecoding_Fully_Centralized_Level4(Hhat_MMSE,F_precoding_Initial,C_MMSE,W,M,K,L,N,p,totIter,nbrOfRealizations,tau_p,tau_c)
warning('off')
%%=============================================================
%This function is used to generate the I-WMMSE precoding matrices and compute the achievable SE for the centralized processing scheme of the paper:
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
SE_sum_Iterative_Monte_MMSE = zeros(totIter,nbrOfRealizations);
SE_Optimal_Monte_MMSE = zeros(K,nbrOfRealizations);
F_Precoding_Optimal_Iterative_Monte_MMSE = zeros(N,N,K,totIter,nbrOfRealizations);

Delta = 5e-4; %Threshold
maxIter = totIter;


[C_MMSE_Total] = functionMatrixGeneration_Total_Error_Matrix(C_MMSE,M,K,L,N);


%-----Monte-Carlo
for i = 1:nbrOfRealizations
    
    SE_sum_objective_monte = 0.0001;
    F_Precoding_Optimal_MMSE = F_precoding_Initial;
    SE_level4_MMSE_Total = zeros(K,totIter);
    
    Hhatallj = reshape(Hhat_MMSE(:,i,:),[M*L K*N]);
    iter = 0;
    
        
        while(iter<maxIter)
            
            iter = iter+1;
            
            SE_sum_objective_old_monte = SE_sum_objective_monte;
            
            [SE_level4_MMSE,E_level4_MMSE,V] = functionCompute_SE_Fully_Centralized_Level4(Hhatallj,F_Precoding_Optimal_MMSE,C_MMSE,tau_c,tau_p,N,L,K,M);
            [GG_level4_MR] = functionMatrixGeneration_Fully_Centralized_Level4(V,Hhatallj,E_level4_MMSE,W,C_MMSE_Total,N,K);
            
            
            %---Monte-Carlo Optimal Precoding
            [F_Precoding_Optimal_MMSE] = functionOptimalObjective_FullyCentralized_Level4(V,Hhatallj,E_level4_MMSE,GG_level4_MR,W,K,N,p);
            F_Precoding_Optimal_Iterative_Monte_MMSE(:,:,:,iter,i) = F_Precoding_Optimal_MMSE;
            
            %--Monte-Carlo
            SE_sum_objective_monte = sum(SE_level4_MMSE.*W);
            SE_sum_Iterative_Monte_MMSE(iter,i) = SE_sum_objective_monte;
            SE_level4_MMSE_Total(:,iter) = SE_level4_MMSE;
            
            if  SE_sum_objective_monte-SE_sum_objective_old_monte<0
                
                SE_sum_Iterative_Monte_MMSE(iter,i) = 0;
                SE_Optimal_Monte_MMSE(:,i) = SE_level4_MMSE_Total(:,iter-1);
                
                1;
                break;
            end
            
            if  abs(SE_sum_objective_monte-SE_sum_objective_old_monte)/SE_sum_objective_old_monte<=Delta
                
                SE_Optimal_Monte_MMSE(:,i) = SE_level4_MMSE;

                break;
            end
            
            SE_Optimal_Monte_MMSE(:,i) = SE_level4_MMSE;

        end
   

end
SE_Optimal_Monte_MMSE_average = mean(SE_Optimal_Monte_MMSE,2);
