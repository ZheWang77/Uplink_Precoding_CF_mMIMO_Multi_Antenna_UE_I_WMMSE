function [SE_level4,E_level4,V] = functionCompute_SE_Fully_Centralized_Level4(Hhatallj,F_precoding,C_MMSE,tau_c,tau_p,N,L,K,M)
warning('off')
%%=============================================================
%This function is used to compute the achievable SE for the centralized processing scheme and generate MMSE combining of the paper:
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


%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);
SE_level4 = zeros(K,1);
E_level4 = zeros(N,N,K);
FF = zeros(N,N,K); % F*F'
FF_total = zeros(K*N,K*N); % F*F'
F_total = zeros(K*N,K*N); % F*F'

for k = 1:K
    
    FF(:,:,k) = F_precoding(:,:,k)*F_precoding(:,:,k)'; %F*F'
    FF_total((k-1)*N+1:k*N,(k-1)*N+1:k*N) = F_precoding(:,:,k)*F_precoding(:,:,k)'; %F*F'
    F_total((k-1)*N+1:k*N,(k-1)*N+1:k*N)  = F_precoding(:,:,k);
    
end

%Store identity matrices of different sizes
eyeML = eye(M*L);

%% Compute sum of all estimation error correlation matrices at every BS
C_MMSE_MMSE_Combining = zeros(L,L,M,K); %Matrix for MMSE Combining

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

C_tot = sum(C_MMSE_MMSE_Combining,4);

C_tot_blk = zeros(M*L,M*L);
for m = 1:M
    
    C_tot_blk(1+(m-1)*L:m*L,1+(m-1)*L:m*L) = C_tot(:,:,m);
    
end

V = ((Hhatallj*FF_total*Hhatallj')+C_tot_blk+eyeML)\(Hhatallj*F_total);

%% Compute SE and MSE Matrix

for k = 1:K
    
    v = V(:,(k-1)*N+1:k*N); %Extract combining matrix
    
    %Compute numerator and denominator of instantaneous SINR at Level 4
    numerator = (v'*Hhatallj(:,(k-1)*N+1:k*N)*F_total((k-1)*N+1:k*N,(k-1)*N+1:k*N));
    denominator = v'*(Hhatallj*FF_total*Hhatallj' + C_tot_blk + eyeML)*v - numerator*numerator';
    
    SE_level4(k) = prelogFactor*real(log2(det(eye(N) + numerator'/denominator*numerator)));
    
    
    E_level4(:,:,k) = eye(N) - v'*Hhatallj(:,(k-1)*N+1:k*N)*F_total((k-1)*N+1:k*N,(k-1)*N+1:k*N) - F_total((k-1)*N+1:k*N,(k-1)*N+1:k*N)'*Hhatallj(:,(k-1)*N+1:k*N)'*v ...
        + v'*(Hhatallj*FF_total*Hhatallj' + C_tot_blk + eyeML)*v;
        
end
