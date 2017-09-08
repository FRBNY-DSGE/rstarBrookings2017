function kf = KF(y,C,R,A,Q,S0,P0);

T = size(y,1);
ns = size(C,2);
Sprev = S0;
Pprev  = P0;

LogLik = 0;

 S = zeros(T,ns)*NaN;
 P = zeros(ns,ns,T)*NaN;
Sf = zeros(T,ns)*NaN;
Pf = zeros(ns,ns,T)*NaN;;
    
    
for t = 1:T
    Sft = A*Sprev;
    Pft = A*Pprev*A'+Q;
    yt = y(t,:)';
    Miss = isnan(yt);
    yt = yt(Miss==0);
    Ct = C(Miss==0,:);
    Rt = R(Miss==0,Miss==0);
    yf = Ct*Sft;
    iV  = inv(Ct*Pft*Ct'+Rt);
    Gain = Pft*Ct'*iV;
    St  = Sft + Gain*(yt-yf);
    Pt  = Pft - Gain*Ct*Pft;
    LogLik = LogLik + .5*log(det(iV))-.5*(yt-yf)'*iV*(yt-yf)-.5*(2*pi);
    S(t,:) = St';
    P(:,:,t) = Pt;
    Sf(t,:) = Sft';
    Pf(:,:,t) = Pft;
    Sprev = St;
    Pprev = Pt;
end;

kf.LogLik = LogLik;
kf.S=S;
kf.P=P;
kf.S0=S0;
kf.P0=P0;
kf.Sf=Sf;
kf.Pf=Pf;
kf.A=A;
kf.Q=Q;
kf.C =C;
kf.R = R;