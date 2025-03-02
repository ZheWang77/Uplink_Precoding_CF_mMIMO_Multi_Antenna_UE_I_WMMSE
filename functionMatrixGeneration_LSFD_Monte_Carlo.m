function [GG] = functionMatrixGeneration_LSFD_Monte_Carlo(H_Combining,H,A,E,W,nbrOfRealizations,N,L,K,M)
%%=============================================================
%This function is used to generate matrices applied in optimizing UL Precoding for the LSFD processing of the paper:
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
%-------Compute terms in optimal F structure

%Go through all channel realizations
for i = 1:nbrOfRealizations

    %-----------------Levels 1-3
    gp_MR = zeros(M*N,N,K,K);
    
    %Go through all APs
    for m = 1:M
        
        %Extract channel realizations from all UEs to AP m
        Hallj = reshape(H(1+(m-1)*L:m*L,i,:),[L K*N]);
        %Extract channel estimate realizations from all UEs to AP m
        Hhatallj = reshape(H_Combining(1+(m-1)*L:m*L,i,:),[L K*N]);
        
        
        %Compute MR/L-MMSE combining
        V = Hhatallj;
        
        for l = 1:K
            
            v = V(:,(l-1)*N+1:l*N); %Extract combining matrix

            for k = 1:K
                
                gp_MR((m-1)*N+1:m*N,:,k,l) = v'*Hallj(:,(k-1)*N+1:k*N);
                
            end
        end
    end

    
    for k = 1:K
        for l = 1:K
            
            
            GG(:,:,k,l) = GG(:,:,k,l) + W(l)*(gp_MR(:,:,k,l)'*A(:,:,l)/E(:,:,l)*A(:,:,l)'*gp_MR(:,:,k,l))/nbrOfRealizations;
            
        end
    end
end