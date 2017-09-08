
% df0tr = 100;
% %           pi^*      r^*        ts^*      liq        safe
% SC0tr =    diag([   2     1/sqrt(2)      1     1/sqrt(4)     1/sqrt(4)   ]).^2/400;

n = length(SC0tr);

nSim = 1e+6; %size(QQ,3);

for jm = 1:nSim
    if mod(jm,1000)==0
        disp(['Now running ',num2str(jm),' out of ',num2str(nSim)])
    end;
    
    [V E]=eig(diag(SC0tr*(df0tr+n+1)));
    %Sinv=V*diag(1./abs(diag(E)))*V';
    eta = randn(df0tr,n)*diag(1./abs(sqrt(diag(E))))*V';
    QQprior(:,:,jm) = (eta'*eta)\eye(n);
end;
