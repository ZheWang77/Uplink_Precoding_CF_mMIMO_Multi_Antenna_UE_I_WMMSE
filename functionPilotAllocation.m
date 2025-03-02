function [Pset] = functionPilotAllocation( R_AP,M,K,N,tau_p,p)
%%=============================================================
%This function is used to generate the pilot allocation strategy of the
%paper:
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

%Pilot set initialize
Pset=1:tau_p;


for z=1:(K/tau_p)-1
    Pset=[Pset;((tau_p*z)+1)*ones(1,tau_p)];
    ind=[];
    for s=1:tau_p
        %Check fot the coherent interference levels
        [coherentx,~] = functionMMSE_interferenceLevels( R_AP,M,tau_p,N,tau_p,p,Pset);
        %Select the UE index that creates least interference
        if s ~=1
            coherentx(ind)=nan;
        end
        [~,ind(s)]=min(coherentx);
        x=1:tau_p;
        x(ind)=[];
        Pset(z+1,x)=(z*tau_p)+s+1;
        
    end
    
end
%Order the pilot allocation set
for i=1:K
    [~,c]=find(Pset==i);
    temp=Pset(:,c);
    temp(temp==i)=[];
    PsetOrdered(:,i)=[i;temp];
end

%The output file
Pset=PsetOrdered;

end

