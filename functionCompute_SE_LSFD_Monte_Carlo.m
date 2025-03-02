function [A,E,signal_MR,SE] = functionCompute_SE_LSFD_Monte_Carlo(H_Combining,H,F_precoding,tau_c,tau_p,nbrOfRealizations,N,L,K,M)
warning('off');
%%=============================================================
%This function is used to compute the achievable SE for the LSFD processing scheme of the paper:
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
signal_MR = zeros(M*N,N,K);
scaling_MR = zeros(M*N,M*N,K);
Gp_MR_total = zeros(M*N,M*N,K,K);
A = zeros(M*N,N,K);
% GG = zeros(N,N,K,K);
E = zeros(N,N,K);

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
                
                gp_MR((m-1)*N+1:m*N,:,k,l) = v'*Hallj(:,(l-1)*N+1:l*N)*F_precoding(:,:,l);
                
            end
        end
    end
    

    
    for k = 1:K
        for l = 1:K
            
            Gp_MR_total(:,:,k,l) = Gp_MR_total(:,:,k,l) + (gp_MR(:,:,k,l)*gp_MR(:,:,k,l)')/nbrOfRealizations;
            
        end
    end
end


Gp_MR = sum(Gp_MR_total,4);
%------------Calculation of LSFD coefficients
for k = 1:K
    
    b = signal_MR(:,:,k)*F_precoding(:,:,k);
    A(:,:,k) = ((Gp_MR(:,:,k)) + scaling_MR(:,:,k))\b;
    E(:,:,k) = eye(N) - b'*A(:,:,k);
    
end

%Compute the SE
for k = 1:K
    
    b = signal_MR(:,:,k)*F_precoding(:,:,k);
    numerator = A(:,:,k)'*b;
    denominator = A(:,:,k)'*(Gp_MR(:,:,k))*A(:,:,k) - numerator*numerator' + A(:,:,k)'*scaling_MR(:,:,k)*A(:,:,k);
    
    SE(k) = prelogFactor*real(log2(det(eye(N) + numerator'/denominator*numerator)));
    
end



