function [C_MMSE_Total] = functionMatrixGeneration_Total_Error_Matrix(C_MMSE,M,K,L,N)
%%=============================================================
%This function is used to generate the useful matrix for the centralized processing scheme of the paper:
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


%Prepare to store the result
C_MMSE_Total = zeros(M*L,M*L,N,N,K); % h_error_kn*h_error_ki'


%Go through all APs          
for m = 1:M
    %Go through all UEs
    for k = 1:K
        for n = 1:N
            for i = 1:N
                
                C_MMSE_Total((m-1)*L+1:m*L,(m-1)*L+1:m*L,n,i,k) = C_MMSE((n-1)*L+1:n*L,(i-1)*L+1:i*L,m,k);
                
            end
        end
    end
end

         
                        

            
