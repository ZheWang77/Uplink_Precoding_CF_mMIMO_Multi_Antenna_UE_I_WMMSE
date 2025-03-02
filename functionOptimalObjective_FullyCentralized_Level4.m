function [F_Precoding_Optimal] = functionOptimalObjective_FullyCentralized_Level4(H_Combining,Hhatallj,E,GG,W,K,N,p)
%%=============================================================
%This function is used to generate I-WMMSE precoding matrices for the
%centralized processing scheme of the paper
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


F_Precoding_Optimal = zeros(N,N,K);

GG_total = sum(GG,4);

for k = 1:K
    

    
    lambda_up = 0;
    lambda_down = 20;
    error = 0.00001;
    
    F_power_up = norm((GG_total(:,:,k) + lambda_up*eye(N))\(W(k)*Hhatallj(:,(k-1)*N+1:k*N)'*H_Combining(:,(k-1)*N+1:k*N)/E(:,:,k)),'fro')^2 - p(k);
   
    
    if F_power_up < 0
        
        1;
        
        lambda = 0;
        F_Precoding_Optimal(:,:,k) = (GG_total(:,:,k) + lambda*eye(N))\(W(k)*Hhatallj(:,(k-1)*N+1:k*N)'*H_Combining(:,(k-1)*N+1:k*N)/E(:,:,k));
         
    else
 
        
    F_power_down = norm((GG_total(:,:,k) + lambda_down*eye(N))\(W(k)*Hhatallj(:,(k-1)*N+1:k*N)'*H_Combining(:,(k-1)*N+1:k*N)/E(:,:,k)),'fro')^2 - p(k);
    
    while(F_power_up * F_power_down > 0)
        
        lambda_down = lambda_down +10;
        
        F_power_down = norm((GG_total(:,:,k) + lambda_down*eye(N))\(W(k)*Hhatallj(:,(k-1)*N+1:k*N)'*H_Combining(:,(k-1)*N+1:k*N)/E(:,:,k)),'fro')^2 - p(k);
        
        if(F_power_up * F_power_down < 0)
            break;
        end
    end
    
    
    %---The bisection method to search the Lagrange multiplier
    
    while(F_power_up * F_power_down < 0)
        
        lambda = (lambda_up + lambda_down)/2;
        
        F_power_now = norm((GG_total(:,:,k) + lambda*eye(N))\(W(k)*Hhatallj(:,(k-1)*N+1:k*N)'*H_Combining(:,(k-1)*N+1:k*N)/E(:,:,k)),'fro')^2 - p(k);
        
        if( F_power_now*F_power_up < 0 )
            
            lambda_down = lambda;
        else
            lambda_up = lambda;
        end
        
        if( abs(lambda_down-lambda_up) < error )
            break;
        end
    end
    
    lambda_optimal = (lambda_down + lambda_up)/2;
    
    F_Precoding_Optimal(:,:,k) = (GG_total(:,:,k) + lambda_optimal*eye(N))\(W(k)*Hhatallj(:,(k-1)*N+1:k*N)'*H_Combining(:,(k-1)*N+1:k*N)/E(:,:,k));
    
    end
end


