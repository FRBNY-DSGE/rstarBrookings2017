%Pi + Yields + AAA + BAA with inflation affecting the term spread
clear;
clc;
close all;
set(0,'defaultAxesFontName', 'Times');
% set(0,'DefaultTextInterpreter', 'latex')
set(0, 'DefaultAxesFontSize',15)
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1.5)
setappdata(0, 'defaultAxesXTickFontSize', 1)
setappdata(0, 'defaultAxesYTickFontSize', 1)

addpath Routines
rng('shuffle')


RunEstimation = 1;
OutputName = 'OutputModel3TsInfl';
FigSubFolder = 'FiguresModel3TsInfl';
if ~exist(FigSubFolder,'dir')
    mkdir(FigSubFolder);
end

if RunEstimation
    Ndraws  =  100000;
    p = 4; %Number of lags in the VAR for the cycle;
    
    if ispc
        [DATA,TEXT] = xlsread('DataCompleteLatest.xls');
        Mnem = TEXT(1,2:end);
        Time = datenum(TEXT(2:end,1),'mm/dd/yyyy');
        Y = DATA;
    else
        [DATA,TEXT] = xlsread('DataCompleteLatest.xls');
        Mnem = TEXT(2:end);
        Time = DATA(:,1) + datenum('12-31-1899');
        Y = DATA(:,2:end);
    end;
    
    T0pre = min(find(year(Time)==1954));
    T1pre = max(find(year(Time)==1959));
    disp(['Avg. and std in the presample: 1954-1959'])
    disp(['' ,Mnem])
    disp(mean(Y(T0pre:T1pre,:)))
    disp(std(Y(T0pre:T1pre,:)))
    
    
    FirstY = 1960;
    LastY  = 2016;
    
    T0 = min(find(year(Time)==FirstY));
    T1 = max(find(year(Time)==LastY));
    
    
    Select = [1:7];
    Mnem = Mnem(Select);
    Y = Y(T0:T1,Select);
    
    y=Y;
    Time = Time(T0:T1);
    [T,n] = size(y);
    T70 = max(find(year(Time)==1970));
    Tzlb = max(find(year(Time)==2008));
    y(Tzlb:end,ismember(Mnem,'BILL'))=NaN;
    y(1:T70,2)=NaN;
    %       1       2        3        4       5
    %      pi       r       ts      liq     safe
    Ctr =[
        1       0       0        0       0       ;%1: Inflatiom
        1       0       0        0       0       ;%2: Inflation Expectations
        1       1       0       -1      -1       ;%3: Short term rate (3m TBILL)
        1       1       0       -1      -1       ;%4: Expected Short Rate
        1       1       1       -1      -1       ;%5: Long Rate (20y TBOND)
        1       1       1        0      -1       ;%6: AAA
        1       1       1        0       0       ;%7: BAA
        ];
    
    
    
    r = size(Ctr,2);
    
    b0 = zeros(n*p,n); b0(1:n,1:n) = eye(n)*0;
    
    lexp0 = 10;
    lambda = 0;
    % Trend = hpfilter(Y*iCtr,1600);
    df0tr = 100;
    %              pi^*      r^*        ts^*      liq        safe
    SC0tr =    ([   2     1/sqrt(2)      1     1/sqrt(4)     1/sqrt(4)   ]).^2/400;%var(diff(Trend));
    S0tr =      [   2        1.5         1        .25         .75         ]';
    P0tr = diag([   1         1          1         1           1         ]);
    
    Psi =                  [ 2    1     1    .5      1        1       1       ];
    %                       pi   Epi   bill  Eill   bond     AAA     BAA
    S0cyc = zeros(n*p,1);
    P0cyc =diag(kron(ones(1,p),Psi));
    
    
    Ccyc = zeros(n,n*p); Ccyc(1:n,1:n) = eye(n);
    C = [Ctr Ccyc];
    
    Atr  = eye(r);
    Acyc = zeros(n*p);
    Acyc(n+1:end,1:end-n) = eye(n*(p-1));
    
    A = zeros(r+n*p);
    A(1:r,1:r) = Atr;
    A(r+1:end,r+1:end) = Acyc;
    
    R = eye(n)*1e-12;
    
    Q0cyc = zeros(n*p); Q0cyc(1:n,1:n) = diag(Psi);
    Q0tr  = diag(SC0tr);
    
    Qcyc = Q0cyc;
    Qtr  = Q0tr;
    Q = zeros(r+n*p);
    Q(1:r,1:r) = Qtr;
    Q(r+1:end,r+1:end) = Qcyc;
    
    
    S0 = [S0tr;S0cyc];
    P0 = zeros(r+n*p);P0(1:r,1:r)=P0tr; P0(r+1:end,r+1:end)=P0cyc;
    
    tic
    
    P_acc = ones(1,Ndraws)*NaN;
    States = ones(T,r+n*p,Ndraws)*NaN;
    Trends = ones(T,n,Ndraws)*NaN;
    LogLik = ones(1,Ndraws)*NaN;
    SS0 = ones(r,Ndraws)*NaN;
    Theta = ones(1,Ndraws);
    AA = ones(r+n*p,r+n*p,Ndraws)*NaN;
    QQ = ones(r+n*p,r+n*p,Ndraws)*NaN;
    CC = ones(n,r+n*p,Ndraws)*NaN;
    RR  = ones(n,n,Ndraws)*NaN;
    
    
    
    for jm = 1:Ndraws
        
        kf = KF(y,C,R,A,Q,S0,P0);
        loglik = kf.LogLik;
        
        
        lambda_new = lambda + randn/sqrt(lexp0);
        
        if lambda_new >=0
            C_new = C; C_new(5:7,1) = 1+lambda_new;
            
            kf_new = KF(y,C_new,R,A,Q,S0,P0);
            loglik_new = kf_new.LogLik;
            log_rat =   (loglik_new - lexp0*lambda_new) ...
                - (loglik     - lexp0*lambda);
            
            p_acc = min(exp(log_rat),1);
            if rand<=p_acc
                lambda = lambda_new;
                C = C_new;
                loglik = loglik_new;
                kf = kf_new;
            end;
            
        else
            p_acc =0;
        end;
        
        
        kc = KC(kf);
        Ycyc = kc.S(:,r+1:r+n);
        for jp=1:p
            Ycyc = [kc.S0(r+(jp-1)*n+1:r+n*jp)'; Ycyc];
        end;
        [beta,sigma] = BVAR(Ycyc,p,b0,Psi,.2,1);
        A(r+1:r+n,r+1:end) = beta';
        Q(r+1:r+n,r+1:r+n) = sigma;
        Ytr = [kc.S0(1:r)'; kc.S(:,1:r)];
        
        %     for jr =1:size(Ytr,2)
        %         dSCtr(jr,:) = CovarianceDraw(diff(Ytr(:,jr)),df0tr,SC0tr(jr));
        %     end;
        %     SCtr = diag(dSCtr);
        
        SCtr = CovarianceDraw(diff(Ytr),df0tr,diag(SC0tr));
        
        Q(1:r,1:r) = SCtr;
        vecP0full = pinv(eye((r+n*p)^2) - kron(A,A))*Q(:);
        P0full = reshape(vecP0full,r+n*p,r+n*p);
        P0(r+1:end,r+1:end) = P0full(r+1:end,r+1:end);
        
        States(:,:,jm) = kc.S;
        Trends(:,:,jm) = kc.S(:,1:r)*C(:,1:r)';
        TrendsReal(:,:,jm) = kc.S(:,2:r)*C(:,2:r)';
        LogLik(jm) = loglik;
        Lambda(:,jm) = lambda;
        AA(:,:,jm) = A;
        QQ(:,:,jm) = Q;
        CC(:,:,jm) = C;
        RR(:,:,jm) = R;
        P_acc(jm)= p_acc;
        
        if mod(jm,50)==0
            disp([num2str(jm),'th draw of ',num2str(Ndraws),'; Elappsed time: ',num2str(toc),' seconds'])
            
            if jm <=1000
                disp(['Acceptance rate in so far: ',num2str(mean(P_acc(1:jm)))])
            elseif jm <10000
                disp(['Acceptance rate of the last 1k draws: ',num2str(mean(P_acc(jm-1000+1:jm)))])
            else
                disp(['Acceptance rate of the last 10k draws: ',num2str(mean(P_acc(jm-1000+1:jm)))])
            end;
            
            
        end;
        
    end;
    
    
    skip = 1;
    Discard = ceil(Ndraws/2);
    
    States = States(:,:,Discard+1:end);
    Trends = Trends(:,:,Discard+1:end);
    TrendsReal = TrendsReal(:,:,Discard+1:end);
    AA = AA(:,:,Discard+1:end);
    QQ = QQ(:,:,Discard+1:end);
    CC = CC(:,:,Discard+1:end);
    RR = RR(:,:,Discard+1:end);
    LogLik = LogLik(:,:,Discard+1:end);
    SS0 = SS0(:,:,Discard+1:end);
    
    CommonTrends = States(:,1:r,:);
    Cycles       = States(:,r+1:r+n,:);
    
    
    save(OutputName,'CommonTrends','Trends','TrendsReal','Cycles','AA','QQ','CC','RR','LogLik','SS0','P_acc','Ndraws','Discard','SC0tr','S0tr','P0tr','df0tr','Psi','Time','Y','y','Mnem')
    
    
else
    
    load(OutputName)
    
end;




Quant = [.025 .16 .5 .84  .975];
sCommonTrends=sort(CommonTrends,3);
sCycles=sort(Cycles,3);
sTrends=sort(Trends,3);
sTrendsReal=sort(TrendsReal,3);


M = size(sCycles,3);
qCommonTrends = sCommonTrends(:,:,ceil(Quant*M));
qCycles = sCycles(:,:,ceil(Quant*M));
qTrends = sTrends(:,:,ceil(Quant*M));
qTrendsReal = sTrendsReal(:,:,ceil(Quant*M));
 
 

 
Pi_bar    =  squeeze(CommonTrends(:,1,:));
M_bar     =  squeeze(CommonTrends(:,2,:));
Ts_bar    =  squeeze(CommonTrends(:,3,:));
Liq_bar   =  squeeze(CommonTrends(:,4,:));
Safe_bar  =  squeeze(CommonTrends(:,5,:));
Cy_bar    =  Liq_bar + Safe_bar;
R_bar     =  M_bar - Cy_bar;

 
sPi_bar    = sort(Pi_bar,2);
sM_bar     = sort(M_bar,2);
sTs_bar    = sort(Ts_bar,2);
sLiq_bar   = sort(Liq_bar,2);
sSafe_bar  = sort(Safe_bar,2);
sCy_bar    = sort(Cy_bar,2);
sR_bar     = sort(R_bar,2);
 
 
qPi_bar    = sPi_bar(:,ceil(Quant*M));
qM_bar     = sM_bar(:,ceil(Quant*M));
qTs_bar    = sTs_bar(:,ceil(Quant*M));
qLiq_bar   = sLiq_bar(:,ceil(Quant*M));  
qSafe_bar  = sSafe_bar(:,ceil(Quant*M));  
qCy_bar    = sCy_bar(:,ceil(Quant*M));
qR_bar     = sR_bar(:,ceil(Quant*M));
 

save OutMod3TsInflforCharts Time qR_bar qCy_bar qM_bar qPi_bar qTs_bar qSafe_bar qLiq_bar y

 
 
Ytr = [Y(:,2) Y(:,3)-Y(:,2)+Y(:,7)-Y(:,5)  Y(:,5)-Y(:,3) Y(:,6)-Y(:,5) Y(:,7)-Y(:,6) ];
 
 
f = figure;
subplot(1,3,1)
PlotStatesShaded(Time,qR_bar)
axis([-inf inf -.5 4.5])
title('R=M-Cy')
subplot(1,3,2)
PlotStatesShaded(Time,qM_bar)
axis([-inf inf -.5 4.5])
title('M')
subplot(1,3,3)
PlotStatesShaded(Time,qCy_bar)
title('Cy')
axis([-inf inf -.5 4.5])
filename=strcat(['./',FigSubFolder,'/Rdecomp.pdf']);
printpdf(f,filename);
 
 
close all
 
f = figure;
PlotStatesShaded(Time,qPi_bar)
filename=strcat(['./',FigSubFolder,'/PIbar.pdf']);
printpdf(f,filename);
title('\pi^*')
 
 
f = figure;
PlotStatesShaded(Time,qPi_bar)
hold on
plot(Time,y(:,2),'b-','linewidth',2.5)
hold on
plot(Time,y(:,1),'b:','linewidth',1)
axis([-inf inf -5 12.5])
hold off
filename=strcat(['./',FigSubFolder,'/PIbar_obs.pdf']);
printpdf(f,filename);
title('\pi^* and \pi')
 
 
 
 
f = figure;
PlotStatesShaded(Time,qM_bar)
filename=strcat(['./',FigSubFolder,'/Mbar.pdf']);
printpdf(f,filename);
title('M^*')
 
 
 
 
f = figure;
PlotStatesShaded(Time,qM_bar)
hold on
plot(Time,y(:,3)-y(:,2)+y(:,7)-y(:,5),'b:','linewidth',1)
hold off
filename=strcat(['./',FigSubFolder,'/Mbar_obs.pdf']);
axis([-inf inf 0 10])
printpdf(f,filename);
title('m^*, r-\pi^e + (baa-r^L)')
 
 
 
f = figure;
PlotStatesShaded(Time,qR_bar)
filename=strcat(['./',FigSubFolder,'/Rbar.pdf']);
printpdf(f,filename);
title('r^*')
 
 
f = figure;
PlotStatesShaded(Time,qR_bar)
hold on
plot(Time,y(:,4)-y(:,2),'b*-','linewidth',2.5)
hold on
plot(Time,y(:,3)-y(:,2),'b:','linewidth',1)
hold off
filename=strcat(['./',FigSubFolder,'/Rbar_obs.pdf']);
printpdf(f,filename);
title('m^*-cy^*, r-\pi^e   and r^e-\pi^e')
 
 
f = figure;
PlotStatesShaded(Time,qTs_bar)
filename=strcat(['./',FigSubFolder,'/TSbar.pdf']);
printpdf(f,filename);
title('Ts^*')
 
 
f = figure;
PlotStatesShaded(Time,qTs_bar)
hold on
plot(Time,y(:,5)-y(:,3),'b:','linewidth',1)
hold off
filename=strcat(['./',FigSubFolder,'/TSbar_obs.pdf']);
printpdf(f,filename);
title('ts^{*}, r-r^L ')
 
 
 
f = figure;
PlotStatesShaded(Time,qCy_bar)
filename=strcat(['./',FigSubFolder,'/CYbar.pdf']);
printpdf(f,filename);
title('CY^*')
 
 
f = figure;
PlotStatesShaded(Time,qCy_bar)
hold on
plot(Time,y(:,7)-y(:,5),'b:','linewidth',1)
hold off
filename=strcat(['./',FigSubFolder,'/CYbar_obs.pdf']);
printpdf(f,filename);
title('CY^*, r^{BAA}-r^L ')
 
 
f = figure;
PlotStatesShaded(Time,qLiq_bar)
filename=strcat(['./',FigSubFolder,'/LIQbar.pdf']);
printpdf(f,filename);
title('CY^*')
 
f = figure;
PlotStatesShaded(Time,qLiq_bar)
hold on
plot(Time,y(:,6)-y(:,5),'b:','linewidth',1)
ylim([0 2])
hold off
filename=strcat(['./',FigSubFolder,'/LIQbar_obs.pdf']);
printpdf(f,filename);
title('Liq^*, r^{AAA}-r^L ')
 
 
 
 
 
f = figure;
PlotStatesShaded(Time,qSafe_bar)
filename=strcat(['./',FigSubFolder,'/SAFEbar.pdf']);
printpdf(f,filename);
title('Safe^*')
 
f = figure;
PlotStatesShaded(Time,qSafe_bar)
hold on
plot(Time,y(:,7)-y(:,6),'b:','linewidth',1)
hold off
filename=strcat(['./',FigSubFolder,'/SAFEbar_obs.pdf']);
printpdf(f,filename);
title('Safe^*, r^{BAA}-r^{AAA} ')
 
 
 
 
 
 
 
 
 
 
 
f = figure;
PlotStatesShaded(Time,-qLiq_bar)
axis([-inf inf -3 1])
filename=strcat(['./',FigSubFolder,'/LIQscaled.pdf']);
printpdf(f,filename);
 
f = figure;
PlotStatesShaded(Time,-qSafe_bar)
axis([-inf inf -3 1])
filename=strcat(['./',FigSubFolder,'/SAFEscaled.pdf']);
printpdf(f,filename);
 
 
f = figure;
PlotStatesShaded(Time,-qCy_bar)
axis([-inf inf -3 1])
filename=strcat(['./',FigSubFolder,'/CYscaled.pdf']);
printpdf(f,filename);
 
f = figure;
PlotStatesShaded(Time,qM_bar)
axis([-inf inf 1 5])
filename=strcat(['./',FigSubFolder,'/Mscaled.pdf']);
printpdf(f,filename);
 
f = figure;
PlotStatesShaded(Time,qR_bar)
axis([-inf inf  -.5 3.5])
filename=strcat(['./',FigSubFolder,'/Rscaled.pdf']);
printpdf(f,filename);
 
 
 
f = figure;
PlotStatesShaded(Time,qTs_bar)
axis([-inf inf  -.5 3.5])
filename=strcat(['./',FigSubFolder,'/TSscaled.pdf']);
printpdf(f,filename);
 
 
 
 
 
for j = 1:size(qCommonTrends,2)
    f=figure(100+j)
    PlotStatesShaded(Time,squeeze(qCommonTrends(:,j,:)))
    hold on
    plot(Time,Ytr(:,j))
    hold off
    title(num2str(j))
    filename=strcat(['./',FigSubFolder,'/fig',num2str(j),'.pdf']);
    printpdf(f,filename);
end;

for j = 1:size(qTrends,2)
    f=figure(110+j)
    subplot(2,1,1)
    PlotStatesShaded(Time,squeeze(qTrends(:,j,:)))
    hold on
    plot(Time,Y(:,j),'-*')
    hold off
    title(Mnem{j})
    subplot(2,1,2)
    PlotStatesShaded(Time,squeeze(qCycles(:,j,:)))
    filename=strcat(['./',FigSubFolder,'/fig',num2str(10+j),'.pdf']);
    printpdf(f,filename);
end;


for j = 1:size(qTrends,2)
    f=figure(120+j)
    subplot(2,1,1)
    PlotStatesShaded(Time,squeeze(qTrendsReal(:,j,:)))
    hold on
    plot(Time,Y(:,j)-Y(:,2)*(j<8),'-*')
    hold off
    title(Mnem{j})
    subplot(2,1,2)
    PlotStatesShaded(Time,squeeze(qCycles(:,j,:)))
    filename=strcat(['./',FigSubFolder,'/fig',num2str(20+j),'.pdf']);
    printpdf(f,filename);
end;



 


SimulInvWishart
 
f = figure;
histogram([sqrt(squeeze(QQ(1,1,:)))]*20,30,'Normalization','pdf')
hold on
[dens,grid] = ksdensity([sqrt(squeeze(QQprior(1,1,:)))]*20);
plot(grid,dens);
hold off
axis([0 4 0 inf])
filename=strcat(['./',FigSubFolder,'/HistSpi.pdf']);
print('-dpdf', filename);
 
f = figure;
histogram([sqrt(squeeze(QQ(2,2,:)))]*20,30,'Normalization','pdf')
hold on
[dens,grid] = ksdensity([sqrt(squeeze(QQprior(2,2,:)))]*20);
plot(grid,dens);
hold off
axis([0 4 0 inf])
filename=strcat(['./',FigSubFolder,'/HistSm.pdf']);
print('-dpdf', filename);
 
f = figure;
histogram([sqrt(squeeze(QQ(3,3,:)))]*20,30,'Normalization','pdf')
hold on
[dens,grid] = ksdensity([sqrt(squeeze(QQprior(3,3,:)))]*20);
plot(grid,dens);
hold off
axis([0 4 0 inf])
filename=strcat(['./',FigSubFolder,'/HistSts.pdf']);
print('-dpdf', filename);
 
 
f = figure;
histogram([sqrt(squeeze(QQ(4,4,:)))]*20,30,'Normalization','pdf')
hold on
[dens,grid] = ksdensity([sqrt(squeeze(QQprior(4,4,:)))]*20);
plot(grid,dens);
hold off
axis([0 4 0 inf])
filename=strcat(['./',FigSubFolder,'/HistSliq.pdf']);
print('-dpdf', filename);
 
f = figure;
histogram([sqrt(squeeze(QQ(5,5,:)))]*20,30,'Normalization','pdf')
hold on
[dens,grid] = ksdensity([sqrt(squeeze(QQprior(5,5,:)))]*20);
plot(grid,dens);
hold off
axis([0 4 0 inf])
filename=strcat(['./',FigSubFolder,'/HistSsafe.pdf']);
print('-dpdf', filename);
 
 



tmax = min(find(year(Time)==1998));

disp('Change in R')
dR = R_bar(end,:)-R_bar(tmax,:);
disp(quantile(dR,Quant))

disp('Change in M')
dM = M_bar(end,:)-M_bar(tmax,:);
disp(quantile(dM,Quant))

disp('Change in Cy')
dCy = Cy_bar(end,:)-Cy_bar(tmax,:);
disp(quantile(dCy,Quant))

disp('Change in Liq')
dLiq = Liq_bar(end,:)-Liq_bar(tmax,:);
disp(quantile(dLiq,Quant))

disp('Change in Safe')
dSafe = Safe_bar(end,:)-Safe_bar(tmax,:);
disp(quantile(dSafe,Quant))

