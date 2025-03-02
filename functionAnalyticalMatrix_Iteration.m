function [Gkk,S_k,X_p1_all,X_p2_all_1,X_p2_all_2,X_p3_1,X_p3_2,GG_total] = functionAnalyticalMatrix_Iteration( R_Vec,F_precoding_Pilot,Phi,Omega,M,K,L,N,tau_p,Pset)
%%=============================================================
%This function is used to compute closed-form matrices applied in optimizing UL Precoding for the LSFD processing with MR combining of the paper:
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

Gkk = zeros(M*N,N,K);     % Gkk
S_k = zeros(M*N,M*N,K); % Sk
F_Pre = zeros(L*N,L*N,K);
X_p1_all = zeros(N,N,N,N,K,K,M); % Analytical Matrices in Gama(1) ii i nn n l k m
X_p2_all_1 = zeros(N,N,N,N,K,K,M); % Analytical Matrices in Gama(2) ii i nn n l k m
X_p2_all_2 = zeros(N,N,N,N,N,N,K,K,M); % Analytical Matrices in Gama(2) q1 q2 ii i nn n l k m

X_p3_1 = zeros(N,N,K,K,M); %Term m unequal m'
X_p3_2 = zeros(N,N,K,K,M);%Term m unequal m'
GG_total = zeros(M*N,M*N,N,N,K,K);

for k = 1:K
    
    F_Pre(:,:,k) = kron(F_precoding_Pilot(:,:,k).',eye(L));  %kron(F^T,I)
    
end


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
                            
                            X_p1_all(ii,i,nn,n,l,k,m) = trace(R_Vec((ii-1)*L+1:ii*L,(i-1)*L+1:i*L,m,l)*Omega((nn-1)*L+1:nn*L,(n-1)*L+1:n*L,m,k));
                            
                            
                            
                            if any(l == Pset(:,k))
                                
                                S_mk = R_Vec(:,:,m,k)*F_Pre(:,:,k)'/Phi(:,:,m,k);
                                P_mkl_1 = tau_p*S_mk*(Phi(:,:,m,k) - tau_p*F_Pre(:,:,l)*R_Vec(:,:,m,l)*F_Pre(:,:,l)')*S_mk';
                                P_mkl_2 = S_mk*F_Pre(:,:,l)*R_Vec(:,:,m,l)*F_Pre(:,:,l)'*S_mk';
                                
                                P_mkl_2_sqrt = P_mkl_2^(1/2);
                                R_ml_sqrt = R_Vec(:,:,m,l)^(1/2);
                                
                                
                                S_ml = R_Vec(:,:,m,l)*F_Pre(:,:,l)'/Phi(:,:,m,l);
                                P_mlk_1 = tau_p*S_ml*(Phi(:,:,m,l) - tau_p*F_Pre(:,:,k)*R_Vec(:,:,m,k)*F_Pre(:,:,k)')*S_ml';
                                P_mlk_2 = S_ml*F_Pre(:,:,k)*R_Vec(:,:,m,k)*F_Pre(:,:,k)'*S_ml';
                                P_mlk_2_sqrt = P_mlk_2^(1/2);
                                R_mk_sqrt = R_Vec(:,:,m,k)^(1/2);
                                
                                
                                %                         for i = 1:N
                                %                             for ii = 1:N
                                
                                X_p2_all_1(ii,i,nn,n,l,k,m) = trace(R_Vec((ii-1)*L+1:ii*L,(i-1)*L+1:i*L,m,l)*P_mkl_1((nn-1)*L+1:nn*L,(n-1)*L+1:n*L));
                                
                                for p1 = 1:N
                                    for p2 = 1:N
                                        
                                        % Analytical Matrices in Gama(2) q1 q2 ii i nn n l k m
                                        X_p2_all_2(p1,p2,ii,i,nn,n,l,k,m) = tau_p^2*trace(P_mkl_2_sqrt((p1-1)*L+1:p1*L,(n-1)*L+1:n*L)*R_ml_sqrt((ii-1)*L+1:ii*L,(p2-1)*L+1:p2*L)*R_ml_sqrt((p2-1)*L+1:p2*L,(i-1)*L+1:i*L)*P_mkl_2_sqrt((nn-1)*L+1:nn*L,(p1-1)*L+1:p1*L))...
                                            + tau_p^2*trace(P_mkl_2_sqrt((p1-1)*L+1:p1*L,(n-1)*L+1:n*L)*R_ml_sqrt((ii-1)*L+1:ii*L,(p1-1)*L+1:p1*L))*trace(P_mkl_2_sqrt((nn-1)*L+1:nn*L,(p2-1)*L+1:p2*L)*R_ml_sqrt((p2-1)*L+1:p2*L,(i-1)*L+1:i*L));
                                        
                                        GG_total((m-1)*N+i,(m-1)*N+ii,nn,n,k,l) = GG_total((m-1)*N+i,(m-1)*N+ii,nn,n,k,l) + tau_p^2*trace(P_mlk_2_sqrt((p1-1)*L+1:p1*L,(i-1)*L+1:i*L)*R_mk_sqrt((n-1)*L+1:n*L,(p2-1)*L+1:p2*L)*R_mk_sqrt((p2-1)*L+1:p2*L,(nn-1)*L+1:nn*L)*P_mlk_2_sqrt((ii-1)*L+1:ii*L,(p1-1)*L+1:p1*L))...
                                            + tau_p^2*trace(P_mlk_2_sqrt((p1-1)*L+1:p1*L,(i-1)*L+1:i*L)*R_mk_sqrt((n-1)*L+1:n*L,(p1-1)*L+1:p1*L))*trace(P_mlk_2_sqrt((ii-1)*L+1:ii*L,(p2-1)*L+1:p2*L)*R_mk_sqrt((p2-1)*L+1:p2*L,(nn-1)*L+1:nn*L));
                                    end
                                end
                                
                                GG_total((m-1)*N+i,(m-1)*N+ii,nn,n,k,l) = GG_total((m-1)*N+i,(m-1)*N+ii,nn,n,k,l) + trace(R_Vec((n-1)*L+1:n*L,(nn-1)*L+1:nn*L,m,k)*P_mlk_1((ii-1)*L+1:ii*L,(i-1)*L+1:i*L));
                                
                                for mm = 1:M
                                    if mm~= m
                                        
                                        RR11 =  tau_p*R_Vec(:,:,m,k)*F_Pre(:,:,k)'/Phi(:,:,m,k)*F_Pre(:,:,l)*R_Vec(:,:,m,l);
                                        RR22 =  tau_p*R_Vec(:,:,mm,l)*F_Pre(:,:,l)'/Phi(:,:,mm,k)*F_Pre(:,:,k)*R_Vec(:,:,mm,k);
                                        
                                        GG_total((m-1)*N+i,(mm-1)*N+ii,nn,n,k,l) = trace(RR11((n-1)*L+1:n*L,(i-1)*L+1:i*L))*trace(RR22((ii-1)*L+1:ii*L,(nn-1)*L+1:nn*L));
                                        
                                    end
                                end
                                
                                
                                
                                RR1 =  R_Vec(:,:,m,l)*F_Pre(:,:,l)'/Phi(:,:,m,k)*F_Pre(:,:,k)*R_Vec(:,:,m,k);
                                RR2 =  R_Vec(:,:,m,k)*F_Pre(:,:,k)'/Phi(:,:,m,k)*F_Pre(:,:,l)*R_Vec(:,:,m,l);
                                
                                X_p3_1(n,nn,k,l,m) = tau_p*trace(RR1((nn-1)*L+1:nn*L,(n-1)*L+1:n*L));
                                X_p3_2(n,nn,k,l,m) = tau_p*trace(RR2((nn-1)*L+1:nn*L,(n-1)*L+1:n*L));
                                
                                
                                
                                clear RR1 RR2 S_mk F_mkl_1 F_mkl_2 F_mkl_2_sqrt R_ml_sqrt
                                
                            else
                                GG_total((m-1)*N+i,(m-1)*N+ii,nn,n,k,l) = trace(R_Vec((n-1)*L+1:n*L,(nn-1)*L+1:nn*L,m,k)*Omega((ii-1)*L+1:ii*L,(i-1)*L+1:i*L,m,l));
                                
                            end
                        end
                    end
                end
            end
        end
        S_k((m-1)*N+1:m*N,(m-1)*N+1:m*N,k) = Gkk((m-1)*N+1:m*N,:,k);
    end
end
