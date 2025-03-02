function [H_Vec,R_Vec] = functionChannelVectorization(H,R_AP,R_UE,M,K,N,L,nbrOfRealizations)

H_Vec = zeros(L*N,nbrOfRealizations,M,K);
R_Vec = zeros(L*N,L*N,M,K);

HH = permute(H,[1,3,2]);

for m = 1:M
    for k = 1:K
        
       H_Vec(:,:,m,k) = reshape(HH((m-1)*L+1:m*L,(k-1)*N+1:k*N,:),L*N,nbrOfRealizations);
       R_Vec(:,:,m,k) = kron(R_UE(:,:,m,k).',R_AP(:,:,m,k));
%        R_Vec(:,:,m,k) = kron(R_UE(:,:,m,k),R_AP(:,:,m,k));
    end
end