function [GG] = functionMatrixGeneration_Fully_Centralized_Level4(H_Combining,Hhatallj,E_level4,W,C_MMSE_Total,N,K)
%%=============================================================
%This function is used to generate matrices applied in optimizing UL Precoding for the centralized processing of the paper:
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
GG_hat = zeros(N,N,K,K);
GG_error = zeros(N,N,K,K);

for k = 1:K
    for l = 1:K
        
        GG_hat(:,:,k,l) = W(l)*Hhatallj(:,(k-1)*N+1:k*N)'*H_Combining(:,(l-1)*N+1:l*N)/E_level4(:,:,l)*H_Combining(:,(l-1)*N+1:l*N)'*Hhatallj(:,(k-1)*N+1:k*N);
        
        for n = 1:N
            for i = 1:N
                
                GG_error(i,n,k,l) = W(l)*trace(H_Combining(:,(l-1)*N+1:l*N)/E_level4(:,:,l)*H_Combining(:,(l-1)*N+1:l*N)'*C_MMSE_Total(:,:,n,i,k));
                
            end
        end
    end
end

GG = GG_hat + GG_error;