function [Hhat_MMSE] = functionChannelEstimates_MMSE(H_Vec,R_Vec,F_precoding,nbrOfRealizations,M,K,L,N,tau_p,Pset)
%%=============================================================
%The file is used to generate the MMSE channel estimates of the paper:
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


%Prepare to store MMSE channel estimates
Hhat_MMSE_vec = zeros(L*N,nbrOfRealizations,M,K);
Hhat_MMSE = zeros(M*L,nbrOfRealizations,K*N);
    
%Store identity matrix of size N x N
eyeLN = eye(L*N);

%Generate realizations of normalized noise 
Np = sqrt(0.5)*(randn(N*L,nbrOfRealizations,M,K) + 1i*randn(N*L,nbrOfRealizations,M,K));

%Generate the precoding matrix
F_Pre = zeros(L*N,L*N,K);

for k = 1:K
    
    F_Pre(:,:,k) = kron(F_precoding(:,:,k).',eye(L));  %kron(F^T,I)
    
end
    
for m = 1:M
    for k = 1:K
        
        yp = zeros(L*N,nbrOfRealizations);
%         yMean = zeros(L*N,nbrOfRealizations);
        PsiInv = zeros(L*N,L*N);
        inds = Pset(:,k); 

        
        for z = 1:length(inds)
            
            yp = yp + tau_p*F_Pre(:,:,inds(z))*H_Vec(:,:,m,inds(z));
%             yMean = yMean + sqrt(p(inds(z)))*tau_p*HMean((m-1)*N+1:m*N,:,inds(z));
            PsiInv = PsiInv + tau_p*F_Pre(:,:,inds(z))*R_Vec(:,:,m,inds(z))*F_Pre(:,:,inds(z))';
            
        end
        
        yp = yp + sqrt(tau_p)*Np(:,:,m,k);
        PsiInv = PsiInv + eyeLN;
 
      
        for z = 1:length(inds)
            
            RPsi = R_Vec(:,:,m,inds(z))*F_Pre(:,:,inds(z))'/PsiInv;
            Hhat_MMSE_vec(:,:,m,inds(z)) = RPsi*yp;
            
        end
    end
end


for m = 1:M
    for k = 1:K
        for n = 1:N
            
            Hhat_MMSE((m-1)*L+1:m*L,:,(k-1)*N+n) = Hhat_MMSE_vec((n-1)*L+1:n*L,:,m,k);
            
        end
    end
end
        




