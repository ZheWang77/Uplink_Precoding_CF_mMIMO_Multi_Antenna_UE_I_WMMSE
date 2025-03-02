function [HH,R_Vec,H_Vec,Omega,F_precoding] = functionChannelGeneration(channelGain,M,K,N,L,p_u,nbrOfRealizations)
%%=============================================================
%The file is used to generate the Weichselberger model based channel of the paper:
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

%Prepare to store the results   
U_APR = zeros(M*L,M*L,K);
U_UET = zeros(K*N,K*N,M);

%---UL Precoding
F_precoding = zeros(N,N,K);

Omega = zeros(M*L,K*N);
WW = sqrt(0.5)*(randn(M*L,K*N,nbrOfRealizations)+1i*randn(M*L,K*N,nbrOfRealizations));
R_Vec = zeros(L*N,L*N,M,K);
H_Vec = zeros(L*N,nbrOfRealizations,M,K);
HH = zeros(M*L,nbrOfRealizations,K*N);


for k = 1:K
    for m = 1:M
        
        %----Unitary Matrix
        AA = rand(L,L) + 1i*rand(L,L);
        BB = rand(N,N) + 1i*rand(N,N);
        [U_APR((m-1)*L+1:m*L,(m-1)*L+1:m*L,k),~,~] = svd(AA*AA');
        [U_UET((k-1)*N+1:k*N,(k-1)*N+1:k*N,m),~,~] = svd(BB*BB'); 
        

        
        %----Coupling Matrix
        Main_Direction = 0.8 + 0.15*rand(1); % One direction capturing dominant channel power
        C1 = sort(rand(L,1),'descend'); %--Dominant Path
        C1 = C1./sum(C1(:))*L*N*Main_Direction;

        C2 = rand(L,N-1);
        C2 = C2./sum(C2(:))*L*N*(1-Main_Direction);
        Omega0 = [C1,C2];
        
        Omega((m-1)*L+1:m*L,(k-1)*N+1:k*N) = channelGain(m,k)*Omega0; %--
           
    end
end


%---Channel Generation
H = permute((sqrt(Omega).*WW),[1,3,2]);
for k = 1:K
    for i = (k-1)*N+1:k*N
        
        HH(:,:,i) = U_APR(:,:,k)*H(:,:,i);
        
    end
end

for m = 1:M
    for t = (m-1)*L+1:m*L
        
        H11 = reshape(HH(t,:,:),nbrOfRealizations,K*N);
        HH(t,:,:) = reshape(H11*U_UET(:,:,m)',nbrOfRealizations,K*N);
        
    end
end



%---Full Correlation Matrix & Channel Vectorization
HH_v = permute(HH,[1,3,2]);
for m = 1:M
    for k = 1:K
        
        Omeg = reshape(Omega((m-1)*L+1:m*L,(k-1)*N+1:k*N),L*N,1);
        R_Vec(:,:,m,k) = kron(conj(U_UET((k-1)*N+1:k*N,(k-1)*N+1:k*N,m)),U_APR((m-1)*L+1:m*L,(m-1)*L+1:m*L,k))*diag(Omeg)*kron(conj(U_UET((k-1)*N+1:k*N,(k-1)*N+1:k*N,m)),U_APR((m-1)*L+1:m*L,(m-1)*L+1:m*L,k))';
        H_Vec(:,:,m,k) = reshape(HH_v((m-1)*L+1:m*L,(k-1)*N+1:k*N,:),L*N,nbrOfRealizations);
        
    end
end


UU_t = sum(U_UET,3);

for k = 1:K
    
    F_precoding(:,:,k) = sqrt(p_u(k))*1/norm(UU_t((k-1)*N+1:k*N,(k-1)*N+1:k*N),'fro')*UU_t((k-1)*N+1:k*N,(k-1)*N+1:k*N);
    
end

