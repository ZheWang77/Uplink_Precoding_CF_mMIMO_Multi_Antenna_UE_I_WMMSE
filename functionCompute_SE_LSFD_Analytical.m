function [A,E,Gkk,SE] = functionCompute_SE_LSFD_Analytical( R_Vec,F_precoding,F_precoding_Pilot,Phi,Omega,M,K,L,N,tau_p,tau_c,Pset)

%---This function is used to computes the theoretical uplink SE for
%MMSE estimator of level3.
%Each AP is equipped with L antennas and each UE is equipped with N antennas.
%This is version 1.0 (Last edited: 2021-02-25)


%INPUT:
%R_AP                 = Matrix with dimension N x N x M x K  where(:,:,m,k) is
%                       the spatial correlation matrix between AP m and UE k
%                       in setup n, normalized by the noise power
%                       normalized by the noise power
%Lk                   = Matrix with dimension MN x MN x K
%
%Phi                  = Matrix with dimension MN x MN x K
%
%Omega                = Matrix with dimension MN x MN x K
%
%A                    = Diagonal matrix with dimension M x M x K where (:,:,k)
%                       is the LSFD coefficients of UE k (when MMSE
%                       estimator is used.)
%M                   = Number of APs
%K                   = Number of UEs
%p                   = 1xK vector, uplink power at each UE
%tau_p               = Pilot length
%tau_c               = Length of the coherence block
%Pset                = Pilot allocation set
%
%
%OUTPUT:
%
%SE_CC              = Vector with dimension K x 1 where (k) is the SE of UE k


%Compute the pre-log factor
%assuming only uplink transmission
prelogFactor = (tau_c-tau_p)/(tau_c);

%Prepare to store the result
% Phi = zeros(L*N,L*N,M,K);
% Omega = zeros(L*N,L*N,M,K);
F_Pre = zeros(L*N,L*N,K);
FF = zeros(N,N,K); % F*F'
Gkk = zeros(M*N,N,K);     % Gkk
Gama_kl = zeros(M*N,M*N,K,K); % Gkl
A = zeros(M*N,N,K); % LSFD coefficient matrix
S_k = zeros(M*N,M*N,K); % Sk
E = zeros(N,N,K); % MSE matrix
SE = zeros(K,1); 


X_p1 = zeros(N,N,K,K,M); %Term Gama(1)
X_p2 = zeros(N,N,K,K,M); %Term Gama(2)
X_p3_1 = zeros(N,N,K,K,M); %Term m unequal m'
X_p3_2 = zeros(N,N,K,K,M);%Term m unequal m'

%%
%----Generate Matrices Applied Below

for k = 1:K
    
    FF(:,:,k) = F_precoding(:,:,k)*F_precoding(:,:,k)'; %F*F'
    F_Pre(:,:,k) = kron(F_precoding_Pilot(:,:,k).',eye(L));  %kron(F^T,I)
    
end

% for m = 1:M
%     
%     for k = 1:K
%         
%         %Compute the UEs indexes that use the same pilot as UE k
%         inds = Pset(:,k);
%         PsiInv = zeros(L*N,L*N);
%         
%         %Go through all UEs that use the same pilot as UE k
%         for z = 1:length(inds)
%             
%             PsiInv = PsiInv + tau_p*F_Pre(:,:,inds(z))*R_Vec(:,:,m,inds(z))*F_Pre(:,:,inds(z))';
%             
%         end
%         
%         PsiInv = PsiInv + eye(L*N);
%         
%         for z = 1:length(inds)
%             
%             Phi(:,:,m,inds(z)) = PsiInv;
%             
%         end
%         
%         Omega(:,:,m,k) = tau_p*R_Vec(:,:,m,k)*F_Pre(:,:,k)'/PsiInv*F_Pre(:,:,k)*R_Vec(:,:,m,k);
%         
%     end
% end
%%
%----Compute G A E in closed-form

%---Compute G_kk G_kl
for k = 1:K
    
    for m = 1:M
        
        for n = 1:N
            for nn = 1:N
                
                Gkk((m-1)*N+n,nn,k) = trace(Omega((nn-1)*L+1:nn*L,(n-1)*L+1:n*L,m,k)); %Term Gkk
                
                for l = 1:K
                    
                    for i = 1:N
                        for ii = 1:N   % i'
                            
                            X_p1(n,nn,k,l,m) = X_p1(n,nn,k,l,m) + FF(ii,i,l)*trace(R_Vec((ii-1)*L+1:ii*L,(i-1)*L+1:i*L,m,l)*Omega((nn-1)*L+1:nn*L,(n-1)*L+1:n*L,m,k)); %Term Gama(1)
                            
                        end
                    end
                    
                    
                    if any(l == Pset(:,k))
                        
                        S_mk = R_Vec(:,:,m,k)*F_Pre(:,:,k)'/Phi(:,:,m,k);
                        P_mkl_1 = tau_p*S_mk*(Phi(:,:,m,k) - tau_p*F_Pre(:,:,l)*R_Vec(:,:,m,l)*F_Pre(:,:,l)')*S_mk';
                        P_mkl_2 = S_mk*F_Pre(:,:,l)*R_Vec(:,:,m,l)*F_Pre(:,:,l)'*S_mk';
                        
                        P_mkl_2_sqrt = P_mkl_2^(1/2);
                        R_ml_sqrt = R_Vec(:,:,m,l)^(1/2);
                        
                        for i = 1:N
                            for ii = 1:N
                                
                                for p1 = 1:N
                                    for p2 = 1:N
                                        
                                        X_p2(n,nn,k,l,m) = X_p2(n,nn,k,l,m) + FF(ii,i,l)*tau_p^2*trace(P_mkl_2_sqrt((p1-1)*L+1:p1*L,(n-1)*L+1:n*L)*R_ml_sqrt((ii-1)*L+1:ii*L,(p2-1)*L+1:p2*L)*R_ml_sqrt((p2-1)*L+1:p2*L,(i-1)*L+1:i*L)*P_mkl_2_sqrt((nn-1)*L+1:nn*L,(p1-1)*L+1:p1*L))...
                                            + FF(ii,i,l)*tau_p^2*trace(P_mkl_2_sqrt((p1-1)*L+1:p1*L,(n-1)*L+1:n*L)*R_ml_sqrt((ii-1)*L+1:ii*L,(p1-1)*L+1:p1*L))*trace(P_mkl_2_sqrt((nn-1)*L+1:nn*L,(p2-1)*L+1:p2*L)*R_ml_sqrt((p2-1)*L+1:p2*L,(i-1)*L+1:i*L));
                                        
                                    end
                                end
                                
                                X_p2(n,nn,k,l,m) = X_p2(n,nn,k,l,m) + FF(ii,i,l)*trace(R_Vec((ii-1)*L+1:ii*L,(i-1)*L+1:i*L,m,l)*P_mkl_1((nn-1)*L+1:nn*L,(n-1)*L+1:n*L));
                            end
                        end
                    
                    
                    RR1 =  R_Vec(:,:,m,l)*F_Pre(:,:,l)'/Phi(:,:,m,k)*F_Pre(:,:,k)*R_Vec(:,:,m,k);
                    RR2 =  R_Vec(:,:,m,k)*F_Pre(:,:,k)'/Phi(:,:,m,k)*F_Pre(:,:,l)*R_Vec(:,:,m,l);
                    
                    
                    X_p3_1(n,nn,k,l,m) = tau_p*trace(RR1((nn-1)*L+1:nn*L,(n-1)*L+1:n*L));
                    X_p3_2(n,nn,k,l,m) = tau_p*trace(RR2((nn-1)*L+1:nn*L,(n-1)*L+1:n*L));
                    
                    clear RR1 RR2 S_mk F_mkl_1 F_mkl_2 F_mkl_2_sqrt R_ml_sqrt
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

