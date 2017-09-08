function Sigma = CovarianceDraw(z,df0,mS0)

 n = size(z,2);
 T = size(z,1);
 Sc0 = mS0*(df0+n+1); %From mode to scale
 
 S=(z'*z)+Sc0;
 
[V E]=eig(S);
%Sinv=V*diag(1./abs(diag(E)))*V';
eta = randn(T+df0,n)*diag(1./abs(sqrt(diag(E))))*V';
Sigma=(eta'*eta)\eye(n);