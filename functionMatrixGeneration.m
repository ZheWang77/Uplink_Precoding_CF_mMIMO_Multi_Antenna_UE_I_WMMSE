function [Phi,Omega,C_MMSE] = functionMatrixGeneration(R_Vec,F_Precoding_Pilot,M,K,L,N,tau_p,Pset)
%%=============================================================
%This function is used to generate the matrix used in the next
%subsequent calculation of the paper:
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
Phi = zeros(L*N,L*N,M,K);
Omega = zeros(L*N,L*N,M,K);
C_MMSE = zeros(L*N,L*N,M,K);
%Generate the precoding matrix
F_Pre = zeros(L*N,L*N,K);

for k = 1:K
    
    F_Pre(:,:,k) = kron(F_Precoding_Pilot(:,:,k).',eye(L));  %kron(F^T,I)
    
end


%Go through all APs          
for m = 1:M
    
    %Go through all UEs
    for k = 1:K
        
        %Compute the UEs indexes that use the same pilot as UE k
        inds = Pset(:,k);
        PsiInv = zeros(L*N,L*N);
        
        %Go through all UEs that use the same pilot as UE k 
        for z = 1:length(inds)   
            
            PsiInv = PsiInv + tau_p*F_Pre(:,:,inds(z))*R_Vec(:,:,m,inds(z))*F_Pre(:,:,inds(z))';


        end
            PsiInv = PsiInv + eye(L*N);

            
            for z = 1:length(inds)
                
                Phi(:,:,m,inds(z)) = PsiInv;
            
            end
            
            Omega(:,:,m,k) = tau_p*R_Vec(:,:,m,k)*F_Pre(:,:,k)'/PsiInv*F_Pre(:,:,k)*R_Vec(:,:,m,k);
            
    end
end

%Generate estimation error correlation matrices
for k = 1:K
    
    C_MMSE(:,:,:,k) = R_Vec(:,:,:,k) - Omega(:,:,:,k);
   
end
                        
                        

            
