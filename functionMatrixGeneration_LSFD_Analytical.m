function [GG] = functionMatrixGeneration_LSFD_Analytical(A,E,GG_total,W,K,N)
%%=============================================================
%This function is used to compute the required closed-form matrices for the LSFD scheme with MR combining of the paper:
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
GG = zeros(N,N,K,K);


for k = 1:K
    for l = 1:K
        for n = 1:N
            for i = 1:N
                
                GG(i,n,k,l) = W(l)*trace(A(:,:,l)/E(:,:,l)*A(:,:,l)'*GG_total(:,:,i,n,k,l));
                
            end
        end
    end
end
