function [A,E,SE] = functionCompute_SE_LSFD_Monte_Carlo_for_Iteration(Gp_MR_all,signal_MR,scaling_MR,F_precoding,tau_c,tau_p,N,K,M)
%%=============================================================
%This function is used to compute the achievable SE in each iteration for the LSFD processing scheme of the paper:
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

%Prepare to store simulation results
SE = zeros(K,1);
FF = zeros(N,N,K); % F*F'
Gp_MR_total = zeros(M*N,M*N,K,K);
A = zeros(M*N,N,K);
E = zeros(N,N,K);


%%
%----Generate Matrices Applied Below

for k = 1:K
    
    FF(:,:,k) = F_precoding(:,:,k)*F_precoding(:,:,k)'; %F*F'
    
end

for k = 1:K
    for l = 1:K
        for m = 1:M
            for mm = 1:M
                
                for n = 1:N
                    for nn = 1:N
                        for ii = 1:N %i
                            for iii = 1:N   % i'
                                
                                Gp_MR_total((m-1)*N+n,(mm-1)*N+nn,k,l) = Gp_MR_total((m-1)*N+n,(mm-1)*N+nn,k,l) + FF(iii,ii,l)*Gp_MR_all(iii,ii,nn,n,l,k,m,mm);
                                
                            end
                        end
                    end
                end
            end
        end
    end
end
                                
Gp_MR = sum(Gp_MR_total,4);


%------------Calculation of LSFD coefficients
for k = 1:K
    
    b = signal_MR(:,:,k)*F_precoding(:,:,k);
    A(:,:,k) = ((Gp_MR(:,:,k)) + scaling_MR(:,:,k))\b;
    E(:,:,k) = eye(N) - b'*A(:,:,k);
    
    numerator = A(:,:,k)'*b;
    denominator = A(:,:,k)'*(Gp_MR(:,:,k))*A(:,:,k) - numerator*numerator' + A(:,:,k)'*scaling_MR(:,:,k)*A(:,:,k);
    
    SE(k) = prelogFactor*real(log2(det(eye(N) + numerator'/denominator*numerator)));
    
end
