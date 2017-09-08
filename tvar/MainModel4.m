%Pi + Yields + AAA + BAA + Consumption
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
OutputName = 'OutputModel4';
FigSubFolder = 'FiguresModel4';
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
    
    
    Select = [1:7 9];
    Mnem = Mnem(Select);
    Y = Y(T0:T1,Select);
    
    y=Y;
    Time = Time(T0:T1);
    [T,n] = size(y);
    T70 = find(year(Time)==1970, 1, 'last' );
    Tzlb = max(find(year(Time)==2008));
    y(Tzlb:end,ismember(Mnem,'BILL'))=NaN;
    y(1:T70,2)=NaN;
    %       1       2        3        4       5     6      7
    %      pi      g       ts      liq     safe    beta    g2
    Ctr =[
           1       0       0        0       0       0      0;%1: Inflatiom
           1       0       0        0       0       0      0;%2: Inflation Expectations
           1       1       0       -1      -1       1      0;%3: Short term rate (3m TBILL)
           1       1       0       -1      -1       1      0;%4: Expected Short Rate
           1       1       1       -1      -1       1      0;%5: Long Rate (20y TBOND)
           1       1       1        0      -1       1      0;%6: AAA
           1       1       1        0       0       1      0;%7: BAA
           0       1       0        0       0       0      1;%8: dC
        ];
    
    
    
    r = size(Ctr,2);
    
    b0 = zeros(n*p,n); b0(1:n,1:n) = eye(n)*0;
    
    lambda = 1;
    % Trend = hpfilter(Y*iCtr,1600);
    df0tr = 100;   
    %              pi^*      G^*           ts^*        liq        safe         beta            g2

    SC0tr =    ([   2     sqrt(1)          1      1/sqrt(4)    1/sqrt(4)  sqrt(1/2*1/4)       sqrt(1) ]).^2/400;%var(diff(Trend));
    S0tr =      [   2         1.5          1          .25        .75           0                  0]';
    P0tr = diag([   1         1            1           1           1           1                  1]);
    
    Psi =                  [ 2    1     1    .5      1        1       1     4     ];
    %                       pi   Epi   bill  Eill   bond     AAA     BAA    dC
    
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
    
    States = ones(T,r+n*p,Ndraws)*NaN;
    Trends = ones(T,n,Ndraws)*NaN;
    LogLik = ones(1,Ndraws)*NaN;
    SS0 = ones(r,Ndraws)*NaN;
    AA = ones(r+n*p,r+n*p,Ndraws)*NaN;
    QQ = ones(r+n*p,r+n*p,Ndraws)*NaN;
    CC = ones(n,r+n*p,Ndraws)*NaN;
    RR  = ones(n,n,Ndraws)*NaN;
    
    
    
    for jm = 1:Ndraws
        
        kf = KF(y,C,R,A,Q,S0,P0);
        loglik = kf.LogLik;
        
        
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
        SS0(:,jm) = S0(1:r);
        AA(:,:,jm) = A;
        QQ(:,:,jm) = Q;
        CC(:,:,jm) = C;
        RR(:,:,jm) = R;
        
        
        if mod(jm,50)==0
            disp([num2str(jm),'th draw of ',num2str(Ndraws),'; Elappsed time: ',num2str(toc),' seconds'])
            
            
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
    
    
    save(OutputName,'CommonTrends','Trends','TrendsReal','Cycles','AA','QQ','CC','RR','LogLik','SS0','Ndraws','Discard','SC0tr','S0tr','P0tr','df0tr','Psi','Time','Y','y','Mnem')
    
    
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


 
 
Pi_bar     =  squeeze(CommonTrends(:,1,:));
G_bar      =  squeeze(CommonTrends(:,2,:));
Ts_bar     =  squeeze(CommonTrends(:,3,:));
Liq_bar    =  squeeze(CommonTrends(:,4,:));
Safe_bar   =  squeeze(CommonTrends(:,5,:));
B_bar      =  squeeze(CommonTrends(:,6,:));
G2_bar     =  squeeze(CommonTrends(:,7,:));
Cy_bar     =  Liq_bar + Safe_bar;
M_bar      =  G_bar + B_bar;
R_bar      =  M_bar - Cy_bar;
DC_bar     =  G_bar + G2_bar;
 

sPi_bar    = sort(Pi_bar,2);
sM_bar     = sort(M_bar,2);
sG_bar     = sort(G_bar,2);
sG2_bar    = sort(G2_bar,2);
sB_bar     = sort(B_bar,2);
sTs_bar    = sort(Ts_bar,2);
sLiq_bar   = sort(Liq_bar,2);
sSafe_bar  = sort(Safe_bar,2);
sCy_bar    = sort(Cy_bar,2);
sR_bar     = sort(R_bar,2);
sDC_bar    = sort(DC_bar,2);
 
 
qPi_bar    = sPi_bar(:,ceil(Quant*M));
qM_bar     = sM_bar(:,ceil(Quant*M));
qG_bar     = sG_bar(:,ceil(Quant*M));
qG2_bar    = sG2_bar(:,ceil(Quant*M));
qB_bar     = sB_bar(:,ceil(Quant*M));
qTs_bar    = sTs_bar(:,ceil(Quant*M));
qLiq_bar   = sLiq_bar(:,ceil(Quant*M));  
qSafe_bar  = sSafe_bar(:,ceil(Quant*M));  
qCy_bar    = sCy_bar(:,ceil(Quant*M));
qR_bar     = sR_bar(:,ceil(Quant*M));
qDC_bar    = sDC_bar(:,ceil(Quant*M));



save OutMod4forCharts Time qR_bar qCy_bar qM_bar qPi_bar qB_bar qTs_bar qSafe_bar qLiq_bar qDC_bar qG2_bar qG_bar y
% load OutMod4forCharts
% Cy = -qCy_bar(:,3);
% G  =  qG_bar(:,3);
% B  =  qB_bar(:,3);
% xlswrite('rstar_dec_trendvar',[Time-1 Cy G B])
% xlswrite('rstar_dec_trendvar',cellstr(datestr(Time-1,'mm/dd/yyyy')))
 
Ytr = [Y(:,2) Y(:,3)-Y(:,2)+Y(:,7)-Y(:,5)  Y(:,5)-Y(:,3) Y(:,6)-Y(:,5) Y(:,7)-Y(:,6) Y(:,8)];

 
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
PlotStatesShaded(Time,qB_bar)
filename=strcat(['./',FigSubFolder,'/Bbar.pdf']);
printpdf(f,filename);
title('B^*')


f = figure;
PlotStatesShaded(Time,qG_bar)
filename=strcat(['./',FigSubFolder,'/Gbar.pdf']);
printpdf(f,filename);
title('G^*')


f = figure;
PlotStatesShaded(Time,qG2_bar)
filename=strcat(['./',FigSubFolder,'/G2bar.pdf']);
printpdf(f,filename);
title('G2^*')

f = figure;
PlotStatesShaded(Time,qDC_bar)
filename=strcat(['./',FigSubFolder,'/DCbar.pdf']);
printpdf(f,filename);
title('DC^*')

 
f = figure;
PlotStatesShaded(Time,qM_bar)
filename=strcat(['./',FigSubFolder,'/Mbar.pdf']);
printpdf(f,filename);
title('M^*')
 
f = figure;
PlotStatesShaded(Time,qR_bar)
filename=strcat(['./',FigSubFolder,'/Rbar.pdf']);
printpdf(f,filename);
title('r^*')
 
 

f = figure;
PlotStatesShaded(Time,qDC_bar)
hold on
plot(Time,y(:,end))
hold off
filename=strcat(['./',FigSubFolder,'/DCbar_obs.pdf']);
printpdf(f,filename);
title('G^*')



f = figure;
PlotStatesShaded(Time,qDC_bar)
hold on
temp = filter(ones(4,1)/4,1,y(:,end)); temp(1:4) = NaN;
plot(Time,temp)
hold off
filename=strcat(['./',FigSubFolder,'/DCbar_obsYoY.pdf']);
printpdf(f,filename);
title('G^*')




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
printpdf(f,filename);
 
f = figure;
histogram([sqrt(squeeze(QQ(2,2,:)))]*20,30,'Normalization','pdf')
hold on
[dens,grid] = ksdensity([sqrt(squeeze(QQprior(2,2,:)))]*20);
plot(grid,dens);
hold off
axis([0 4 0 inf])
filename=strcat(['./',FigSubFolder,'/HistSg.pdf']);
printpdf(f,filename);
 
f = figure;
histogram([sqrt(squeeze(QQ(3,3,:)))]*20,30,'Normalization','pdf')
hold on
[dens,grid] = ksdensity([sqrt(squeeze(QQprior(3,3,:)))]*20);
plot(grid,dens);
hold off
axis([0 4 0 inf])
filename=strcat(['./',FigSubFolder,'/HistSts.pdf']);
printpdf(f,filename); 
 
f = figure;
histogram([sqrt(squeeze(QQ(4,4,:)))]*20,30,'Normalization','pdf')
hold on
[dens,grid] = ksdensity([sqrt(squeeze(QQprior(4,4,:)))]*20);
plot(grid,dens);
hold off
axis([0 4 0 inf])
filename=strcat(['./',FigSubFolder,'/HistSliq.pdf']);
printpdf(f,filename);
 
f = figure;
histogram([sqrt(squeeze(QQ(5,5,:)))]*20,30,'Normalization','pdf')
hold on
[dens,grid] = ksdensity([sqrt(squeeze(QQprior(5,5,:)))]*20);
plot(grid,dens);
hold off
axis([0 4 0 inf])
filename=strcat(['./',FigSubFolder,'/HistSsafe.pdf']);
printpdf(f,filename);


 
f = figure;
histogram([sqrt(squeeze(QQ(6,6,:)))]*20,30,'Normalization','pdf')
hold on
[dens,grid] = ksdensity([sqrt(squeeze(QQprior(6,6,:)))]*20);
plot(grid,dens);
hold off
axis([0 4 0 inf])
filename=strcat(['./',FigSubFolder,'/HistSbeta.pdf']);
printpdf(f,filename);


f = figure;
histogram([sqrt(squeeze(QQ(7,7,:)))]*20,30,'Normalization','pdf')
hold on
[dens,grid] = ksdensity([sqrt(squeeze(QQprior(7,7,:)))]*20);
plot(grid,dens);
hold off
axis([0 4 0 inf])
filename=strcat(['./',FigSubFolder,'/HistSg2.pdf']);
    printpdf(f,filename);




tmax = min(find(year(Time)==1998));


f = figure;
bands(Time,qR_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis([-inf inf 0 3.5])
datetick('x', 17,'keeplimits', 'keepticks')
hold off;
filename=strcat(['./',FigSubFolder,'/Figure_Rbar.pdf']);
printpdf(f,filename);



f = figure;
bands(Time,qR_bar,-qCy_bar-(-qCy_bar(tmax,3)-qR_bar(tmax,3)))
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis([-inf inf 0 3.5])
datetick('x', 17,'keeplimits', 'keepticks')
hold off;
filename=strcat(['./',FigSubFolder,'/Figure_Rbar-Cybar.pdf']);
    printpdf(f,filename);





f = figure;
bands(Time,qR_bar,qM_bar-(qM_bar(tmax,3)-qR_bar(tmax,3)))
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis([-inf inf 0 3.5])
datetick('x', 17,'keeplimits', 'keepticks')
hold off;
filename=strcat(['./',FigSubFolder,'/Figure_Rbar-Mbar.pdf']);
printpdf(f,filename);


f = figure;
bands(Time,qR_bar,qG_bar-(qG_bar(tmax,3)-qR_bar(tmax,3)))
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis([-inf inf 0 3.5])
datetick('x', 17,'keeplimits', 'keepticks')
hold off;
filename=strcat(['./',FigSubFolder,'/Figure_Rbar-Gbar.pdf']);
printpdf(f,filename);




disp('Change in R')
dR = R_bar(end,:)-R_bar(tmax,:);
disp(quantile(dR,Quant))

disp('Change in M')
dM = M_bar(end,:)-M_bar(tmax,:);
disp(quantile(dM,Quant))

disp('Change in G')
dG = G_bar(end,:)-G_bar(tmax,:);
disp(quantile(dG,Quant))

disp('Change in B')
dB = B_bar(end,:)-B_bar(tmax,:);
disp(quantile(dB,Quant))

disp('Change in Cy')
dCy = Cy_bar(end,:)-Cy_bar(tmax,:);
disp(quantile(dCy,Quant))

disp('Change in Liq')
dLiq = Liq_bar(end,:)-Liq_bar(tmax,:);
disp(quantile(dLiq,Quant))

disp('Change in Safe')
dSafe = Safe_bar(end,:)-Safe_bar(tmax,:);
disp(quantile(dSafe,Quant))

disp('Change in Growth')
dC = DC_bar(end,:)-DC_bar(tmax,:);
disp(quantile(dC,Quant))

disp('Change in Hand-to-Mouth')
dG2 = G2_bar(end,:)-G2_bar(tmax,:);
disp(quantile(dG2,Quant))
 


