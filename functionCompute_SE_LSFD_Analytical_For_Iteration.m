function [A,E,SE] = functionCompute_SE_LSFD_Analytical_For_Iteration(Gkk,S_k,X_p1_all,X_p2_all_1,X_p2_all_2,X_p3_1,X_p3_2,F_precoding,M,K,N,tau_p,tau_c,Pset)
warning('off');
%%=============================================================
%This function is used to compute the closed-form achievable SE in each iteration for the LSFD processing scheme with MR combining of the paper:
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


%Compute the pre-log factor
%assuming only uplink transmission
prelogFactor = (tau_c-tau_p)/(tau_c);

%Prepare to store the result
FF = zeros(N,N,K); % F*F'
Gama_kl = zeros(M*N,M*N,K,K); % Gkl
A = zeros(M*N,N,K); % LSFD coefficient matrix
E = zeros(N,N,K); % MSE matrix
SE = zeros(K,1); 
X_p1 = zeros(N,N,K,K,M); %Term Gama(1)
X_p2 = zeros(N,N,K,K,M); %Term Gama(2)


%%
%----Generate Matrices Applied Below

for k = 1:K
    
    FF(:,:,k) = F_precoding(:,:,k)*F_precoding(:,:,k)'; %F*F'
    
end
%%
%----Compute G A E in closed-form

%---Compute G_kk G_kl
for k = 1:K
    
    for m = 1:M
        
        for n = 1:N
            for nn = 1:N

                for l = 1:K
                    
                    for i = 1:N
                        for ii = 1:N   % i'
                            
                            X_p1(n,nn,k,l,m) = X_p1(n,nn,k,l,m) + FF(ii,i,l)*X_p1_all(ii,i,nn,n,l,k,m); %Term Gama(1)
                            
                            if any(l == Pset(:,k))
                                for p1 = 1:N
                                    for p2 = 1:N
                                        
                                        X_p2(n,nn,k,l,m) = X_p2(n,nn,k,l,m) + FF(ii,i,l)*X_p2_all_2(p1,p2,ii,i,nn,n,l,k,m);
                                        
                                    end
                                end
                                
                                X_p2(n,nn,k,l,m) = X_p2(n,nn,k,l,m) + FF(ii,i,l)*X_p2_all_1(ii,i,nn,n,l,k,m);
                                
                            end
                        end
                    end
                end
            end
            
            
        end
    end
end

%----------------------------------
for k = 1:K
    for m = 1:M
        for mm = 1:M
            for l = 1:K
                
                if any(l == Pset(:,k))
                    if any(mm == m)
                        
                        Gama_kl((m-1)*N+1:m*N,(m-1)*N+1:m*N,k,l) = X_p2(:,:,k,l,m);
                    else
                        Gama_kl((m-1)*N+1:m*N,(mm-1)*N+1:mm*N,k,l) = X_p3_1(:,:,k,l,m)*FF(:,:,l)*X_p3_2(:,:,k,l,mm);
                    end
                else
                    
                    if any(mm == m)
                        
                        Gama_kl((m-1)*N+1:m*N,(m-1)*N+1:m*N,k,l) = X_p1(:,:,k,l,m);
                        
                    end
                    
                end
            end
            S_k((m-1)*N+1:m*N,(m-1)*N+1:m*N,k) = Gkk((m-1)*N+1:m*N,:,k);
        end
    end
end


%------------Calculate LSFD coefficient matrix
Gama_kl_total = sum(Gama_kl,4);
for k = 1:K
    
    A(:,:,k) = ((Gama_kl_total(:,:,k)) + S_k(:,:,k))\(Gkk(:,:,k)*F_precoding(:,:,k));
    E(:,:,k) = eye(N) - (Gkk(:,:,k)*F_precoding(:,:,k))'*A(:,:,k);
    
end


%%

%-----------Compute the SE
for k = 1:K
    
    numerator = A(:,:,k)'*Gkk(:,:,k)*F_precoding(:,:,k);
    denominator = A(:,:,k)'*(Gama_kl_total(:,:,k))*A(:,:,k) - numerator*numerator' + A(:,:,k)'*S_k(:,:,k)*A(:,:,k);
    SE(k) = prelogFactor*real(log2(det(eye(N) + numerator'/denominator*numerator)));
    
end

