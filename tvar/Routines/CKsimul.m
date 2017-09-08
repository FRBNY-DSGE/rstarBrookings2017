function [xdraw,xdraw0] = CKsimul(y,CC,RR,AA,QQ,ResKF);

[T,N]= size(y);
ns = size(CC,2); 

ResKS = FIS(y,CC,RR,AA,QQ,ResKF);

xsmooth = ResKS.AmT(:,2:end);
Vsmooth = ResKS.PmT(:,:,2:end);

xsmooth0 = ResKS.AmT(:,1);
Vsmooth0 = ResKS.PmT(:,:,1);


%FF=xsmooth(:,2:end);V=Vsmooth(:,:,2:end);
%init_F1=xsmooth(:,1);init_V1=Vsmooth(:,:,1);
    

xdraw = zeros(ns,T);

xdraw(:,T) = cholred(Vsmooth(:,:,end))'*randn(ns,1) + xsmooth(:,end);


for t = T-1:-1:1
    F_b = xsmooth(:,t)  + Vsmooth(:,:,t)*AA'/(AA*Vsmooth(:,:,t)*AA'+QQ)*(xdraw(:,t+1)-AA*xsmooth(:,t));
    V_b = Vsmooth(:,:,t) - Vsmooth(:,:,t)*AA'/(AA*Vsmooth(:,:,t)*AA'+QQ)*AA*Vsmooth(:,:,t);
    xdraw(:,t) = cholred(V_b)'*randn(ns,1) + F_b; 
end

F_b = xsmooth0  + Vsmooth0*AA'/(AA*Vsmooth0*AA'+QQ)*(xdraw(:,1)-AA*xsmooth0);
V_b = Vsmooth0 - Vsmooth0*AA'/(AA*Vsmooth0*AA'+QQ)*AA*Vsmooth0;
xdraw0 = cholred(V_b)'*randn(ns,1) + F_b;


    
function C = cholred(S)


[C,temp] = chol(S/2+S'/2+1e-6*eye(size(S,1)));
