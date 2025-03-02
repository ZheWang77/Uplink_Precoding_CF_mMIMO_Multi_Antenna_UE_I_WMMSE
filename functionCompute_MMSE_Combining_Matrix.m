function [V_MMSE_Combining,C_MMSE_MMSE_Combining] = functionCompute_MMSE_Combining_Matrix(Hhat,F_precoding_optimal,C_MMSE,nbrOfRealizations,L,N,K,M)
%%=============================================================
%This function is used to generate the L-MMSE combining vector used in the next subsequent calculation of the paper:
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

%Store identity matrices of different sizes
eyeL = eye(L);


%Prepare to store the result
% Phi = zeros(L*N,L*N,M,K);
% Omega = zeros(L*N,L*N,M,K);
% F_Pre = zeros(L*N,L*N,K);
FF = zeros(N,N,K); % F*F'
F_total = zeros(K*N,K*N); % F*F'
FF_total = zeros(K*N,K*N); % F*F'
% C_MMSE = zeros(L*N,L*N,M,K);
C_MMSE_MMSE_Combining = zeros(L,L,M,K); %Matrix for MMSE Combining
%%

for k = 1:K
    
    FF(:,:,k) = F_precoding_optimal(:,:,k)*F_precoding_optimal(:,:,k)'; %F*F'
    FF_total((k-1)*N+1:k*N,(k-1)*N+1:k*N) = F_precoding_optimal(:,:,k)*F_precoding_optimal(:,:,k)'; %F*F'
    F_total((k-1)*N+1:k*N,(k-1)*N+1:k*N) = F_precoding_optimal(:,:,k);
    
    for m = 1:M
        for j = 1:L
            for q = 1:L
                
                for mm = 1:N
                    for p = 1:N
                        
                        C_MMSE_MMSE_Combining(j,q,m,k) = C_MMSE_MMSE_Combining(j,q,m,k) + C_MMSE((p-1)*L+j,(mm-1)*L+q,m,k)*FF(p,mm,k);
                        
                    end
                end
                
            end
        end
    end
end


V_MMSE_Combining = zeros(M*L,nbrOfRealizations,K*N);


%Compute sum of all estimation error correlation matrices at every BS
C_tot = sum(C_MMSE_MMSE_Combining,4);




%% Go through all channel realizations
for n = 1:nbrOfRealizations

    %Go through all APs
    for m = 1:M
        
        
        %Extract channel estimate realizations from all UEs to AP l
        Hhatallj = reshape(Hhat(1+(m-1)*L:m*L,n,:),[L K*N]);

        %Compute MMSE combining
        V_MMSE_Combining((m-1)*L+1:m*L,n,:) = ((Hhatallj*FF_total*Hhatallj') + C_tot(:,:,m)+eyeL)\(Hhatallj*F_total);
        
    end
end