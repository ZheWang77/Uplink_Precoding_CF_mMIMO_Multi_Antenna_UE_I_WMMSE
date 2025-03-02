function [Gp_MR_all,GG_total_Monte,signal_MR,scaling_MR] = functionMonte_Carlo_Matrix_Iteration(H_Combining,H,nbrOfRealizations,N,L,K,M)
%%=============================================================
%This function is used to generate the useful closed-form matrices in each iteration for the LSFD processing scheme of the paper:
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

Gp_MR_all = zeros(N,N,N,N,K,K,M,M); % Monte-Carlo Matrices in Gama ii i nn n l k m
GG_total_Monte = zeros(M*N,M*N,N,N,K,K);
signal_MR = zeros(M*N,N,K);
scaling_MR = zeros(M*N,M*N,K);


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
        
        for k = 1:K
            
            v = V(:,(k-1)*N+1:k*N); %Extract combining matrix
            signal_MR((m-1)*N+1:m*N,:,k) = signal_MR((m-1)*N+1:m*N,:,k) + (v'*Hallj(:,(k-1)*N+1:k*N))/nbrOfRealizations; % (v_mk)'h_mk
            scaling_MR((m-1)*N+1:m*N,(m-1)*N+1:m*N,k) = scaling_MR((m-1)*N+1:m*N,(m-1)*N+1:m*N,k) + v'*v/nbrOfRealizations; % ||v_mk||^2
            
            
            for l = 1:K
                
                v_l = V(:,(l-1)*N+1:l*N);
                h_k = Hallj(:,(k-1)*N+1:k*N);
                h = Hallj(:,(l-1)*N+1:l*N);
                
                gp_MR((m-1)*N+1:m*N,:,k,l) = v_l'*h_k;
                
                for mm = 1:M
                    
                    Hallj_mm = reshape(H(1+(mm-1)*L:mm*L,i,:),[L K*N]);
                    Hhatallj_mm = reshape(H_Combining(1+(mm-1)*L:mm*L,i,:),[L K*N]);
                    
                    h_mm = Hallj_mm(:,(l-1)*N+1:l*N);
                    v_mm = Hhatallj_mm(:,(l-1)*N+1:l*N);
                    
                    for n = 1:N
                        for nn = 1:N
                            for ii = 1:N %i
                                for iii = 1:N   % i'
                                    
                                    Gp_MR_all(iii,ii,nn,n,l,k,m,mm) = Gp_MR_all(iii,ii,nn,n,l,k,m,mm) + v(:,n)'*h(:,iii)*h_mm(:,ii)'*v_mm(:,nn)/nbrOfRealizations;
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    for k = 1:K
        for l = 1:K
            for n = 1:N
                for nn = 1:N
                    
                    GG_total_Monte(:,:,n,nn,k,l) = GG_total_Monte(:,:,n,nn,k,l) + gp_MR(:,nn,k,l)*gp_MR(:,n,k,l)'/nbrOfRealizations;
                    
                end
            end
        end
    end
end











