function kc = KC(kf);

S  = kf.S;
P  = kf.P;
Sf = kf.Sf;
Pf = kf.Pf;
S0 = kf.S0;
P0 = kf.P0;
A = kf.A;
C = kf.C;
R = kf.R;
Q = kf.Q;

[T,ns] = size(S);
drS = zeros(T,ns)*NaN;
mu = S(T,:)';
sigma = P(:,:,T)';
draw = mu + cholred(sigma)'*randn(size(mu));
drS(T,:) = draw';

for t=T-1:-1:1
    iPf = pinv(Pf(:,:,t+1));
    mu =  S(t,:)'+P(:,:,t)*A'*iPf*(drS(t+1,:)-Sf(t+1,:))';
    sigma  =  P(:,:,t) - P(:,:,t)*A'*iPf*A*P(:,:,t);
    draw = mu + cholred(sigma)'*randn(size(mu));
    drS(t,:) = draw';
end;

iPf = inv(Pf(:,:,1));
mu =  S0+P0*A'*iPf*(drS(1,:)'-Sf(1,:)');
sigma  =  P0 - P0*A'*iPf*A*P0;
draw = mu + cholred(sigma)'*randn(size(mu));
drS0 = draw;
kc.S  = drS;
kc.S0 = drS0;



function C = cholred(S)
[v,d] = eig(S);
dd = diag(d);
dd(dd<1e-12)=0;
C = real(diag(sqrt(dd))*v');

%[C,temp] = chol(S/2+S'/2+1e-2*eye(size(S,1)));
