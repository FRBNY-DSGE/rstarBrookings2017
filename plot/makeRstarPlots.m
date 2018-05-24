%%% MAKERSTARPLOTS.M: Plots for r* paper
%
% This script creates figures for the paper "Safety, Liquidity, and the Natural
% Rate of Interest" by Marco Del Negro, Domenico Giannone, Marc Giannoni, and
% Andrea Tambalotti.
%
% The following plots are produced:
%  - DSGE {20,30}-year forward r* vs TVAR rbar (Figure 9)
%  - DSGE {5,10}-year forward r vs r* (Figure 10)
%  - DSGE 5-year forward r* vs Laubach Williams r* (Figure 11)
%  - DSGE 20-year forward r* and its drivers (Figure 12)
%  - DSGE short-term r and r* (Figure 13)
%  - DSGE short-term r* and its drivers (Figure 14)
%
% Created 2017-03-02  Erica Moszkowski and Brandyn Bok
%                     Federal Reserve Bank of New York


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
addpath('helperFunctions');

% Set up output directory for figures
figurespath = 'Figures/';

if ~exist(figurespath)
  mkdir(figurespath);
end

% File suffix common to all files
tvarTablespath = fullfile('..', 'tvar', 'output_data');
dsgeTablespath = fullfile('..', 'dsge', 'output_data', 'm1010', 'ss20', ...
    'forecast', 'tables');
otherTablespath = 'Tables/';
suffix = '_cond=none_para=full_vint=161223.csv';

% Do we want to label shockdecs with trend (adjustLevel = 1) or starting at 0?
adjustLevel = 1;

% Use Bloomberg data?
useBloombergData = 0;

% This numbers the columns of the DSGE output tables (NOT including the dates).
% spec1010_18_2016Q3_1223.jl is set up to print the bands from narrowest on the
% outside to widest on the inside, with the mean in the right-hand
% column (so, in the CSV, the 68% bands are in columns 1 and 4, the 95% bands
% are in columns 2 and 3, and the mean is in column 5). The plotting code
% requires that the mean be in the middle, so we reorder the columns here.
columns = [2,1,5,4,3];

% Start/end dates
shockdecStart = datenum('1960-03-31');
shockdecEnd = datenum('2016-09-30');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ IN DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp  = readtable(fullfile(dsgeTablespath,['hist_RealNaturalRate',suffix]));
time = datenum(tmp.date);

%%%%%%% Smoothed histories
%% Forward real rate smoothed history
% Short term r
tmp = csvread(fullfile(dsgeTablespath,['hist_ExAnteRealRate',suffix]),...
    1, 1);
DSGE_rST = tmp(:, columns);

% 5yr fwd r
tmp = csvread(fullfile(dsgeTablespath,['hist_Forward5YearRealRate',suffix]),...
    1, 1);
DSGE_r5fwd = tmp(:, columns);

% 10yr fwd r
tmp = csvread(fullfile(dsgeTablespath,['hist_Forward10YearRealRate',suffix]),...
    1, 1);
DSGE_r10fwd = tmp(:, columns);

%% Forward real natural rate
% ST r*
tmp = csvread(fullfile(dsgeTablespath,['hist_RealNaturalRate',suffix]),...
    1, 1);
DSGE_rstarST = tmp(:, columns);

% 5yr fwd rstar
tmp = csvread(fullfile(dsgeTablespath,['hist_Forward5YearRealNaturalRate',suffix]),...
    1, 1);
DSGE_rstar5fwd = tmp(:, columns);

% 10yr fwd rstar
tmp = csvread(fullfile(dsgeTablespath,['hist_Forward10YearRealNaturalRate',suffix]),...
    1, 1);
DSGE_rstar10fwd = tmp(:, columns);

% 20yr fwd rstar
tmp = csvread(fullfile(dsgeTablespath,['hist_Forward20YearRealNaturalRate',suffix]),...
    1, 1);
DSGE_rstar20fwd = tmp(:, columns);

% 30yr fwd rstar
tmp = csvread(fullfile(dsgeTablespath,['hist_Forward30YearRealNaturalRate',suffix]),...
    1, 1);
DSGE_rstar30fwd = tmp(:, columns);

%% Trendy VAR
% Import trendy VAR data.
load(fullfile(tvarTablespath, 'trendyVARrbar'));
rstarTrendyVAR = qR_bar(1:end-1,:);

% Import trendy VAR shockdec
TVAR_shockdec = readtable(fullfile(tvarTablespath, 'rstar_dec_trendvar_less1998Q4value.csv'));
TVAR_cybar    = TVAR_shockdec.x_Cy_bar;
TVAR_gbar     = TVAR_shockdec.G_bar;
TVAR_bbar     = TVAR_shockdec.B_bar;
TVAR_rbar     = TVAR_shockdec.r_bar;
TVAR_time     = datenum(TVAR_shockdec.date);

%% Laubach Williams
%% Import Laubach Williams data
rstarLW = csvread([otherTablespath, 'LaubachWilliamsrstar.csv'],1,1);

%% Refcorp

[DATA,TEXT] = xlsread(fullfile(otherTablespath, 'Refcorp.xlsx'),'Sheet1');

TimeRefcorp = datenum(DATA(:,1))+datenum('12-30-1899','mm-dd-yyyy');
Refcorp5Y  = DATA(:,9);
Refcorp10Y = DATA(:,11);
Refcorp20Y = DATA(:,12);


%% Set figure specs
% axes
set(0, 'defaultAxesFontName', 'Times');
set(0, 'defaultAxesFontSize', 20);
set(0, 'defaultLineLineWidth', 1);
setappdata(0, 'defaultAxesXTickFontSize', 1);
setappdata(0, 'defaultAxesYTickFontSize', 1);

% colors
red = [.667 .29 .224]; dark_red = [.502 .149 .082]; light_red = [.831 .482 .416];
blue = [.18 .259 .447]; dark_blue = [.086 .161 .333]; light_blue = [.31 .384 .557];
green = [.161 .486 .275]; dark_green = [.063 .365 .165]; light_green = [.306 .612 .408];
gold = [.667 .518 .224]; yellow = [.906 .929 .271]; light_yellow = [.992, .992, .588];
purple = [.494 .204 .624]; orange = [.937 .506 .275]; sea_green = [.18; .62; .486];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT
%% Intro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Figure 1 - DSGE 30-year forward r* vs TVAR rbar
fig1 = figure(1);
l1 = PlotStatesShadedv3_oneband(time,rstarTrendyVAR(:,2:4));
hold on;

l2 = PlotStatesShadedv3_oneband(time,DSGE_rstar30fwd(:,2:4),[0 0 1], 0.2);
hold off;
leg = legend([l1,l2], 'VAR', 'DSGE', 'location', 'SouthWest','interpreter','latex');
set(leg,'interpreter', 'latex')
ylim([0,3.5]);
set(gca,'YTick',[0:0.5:3.5])
printpdf(fig1,[figurespath, 'Figure1'], 'square', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Figure 2
%%%%%%%% Charts for the model with Inflation, Tbill, Tbond %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(tvarTablespath, 'OutMod1forCharts'));

f = figure('Name','Figure 2(a)','NumberTitle','off');
PlotStatesShaded(Time,qR_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
hold on
plot(Time,y(:,4)-y(:,2),'b*-','linewidth',2.5)
hold on
plot(Time,y(:,3)-y(:,2),'b:','linewidth',1)
axis tight; box on;
plot(Time,Time*0,'k','LineWidth',.25)
ylim([-1 6])
hold off
filename= fullfile(figurespath, 'Figure2a');
printpdf(f,filename, 'square', 1);

f = figure('Name','Figure 2(b)','NumberTitle','off');
PlotStatesShaded(Time,qPi_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
hold on
plot(Time,y(:,2),'b-','linewidth',2.5)
hold on
plot(Time,y(:,1),'b:','linewidth',1)
axis tight; box on;
plot(Time,Time*0,'k','LineWidth',.25)
ylim([-1 10])
hold off
filename = fullfile(figurespath, 'Figure2b');
printpdf(f,filename, 'square', 1);

% request for marco
clear Time
load(fullfile(tvarTablespath, 'OutMod1forCharts'));

f = figure('Name','Figure 2(a) alt','NumberTitle','off');
bands(Time,qR_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis([-inf inf 0 3.5])
hold off;
filename= fullfile(figurespath,'Figure2a_alt');
printpdf(f, filename, 'square', 1);

f = figure('Name','Figure 2(b) alt','NumberTitle','off');
bands(Time,qPi_bar);
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis([-inf inf 0 7.5])
hold off;
filename= fullfile(figurespath,'Figure2b_alt');
printpdf(f, filename, 'square', 1);

%% Figure 3
clear Time
%%%%%%%% Charts for the model with Inflation, Tbill, Tbond, BAA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(tvarTablespath, 'OutMod2forCharts')); % Data from the model with Pi, PiExp, IntShort IntLong and BAA
%
% The data that enter the model are stored in y
% (in the order Inflation, Inlfation Expectations, Short rate, Expected Short Rate, Long Rate, BAA)
%
% The quantiles (.025, .16 .5 .84 .975) of the posterior of the trends are stored in
% qR_bar (trend in the real rate), qPi_bar (trend inflation),  qTs_bar(trend of the Term Spread),
% qCy_bar (trend of the Convenience yield) and qM_bar (trend of the discount factor)

%%Charts of Trends and Observables, Convenience Yield Model

f = figure('Name','Figure 3(a)','NumberTitle','off');
PlotStatesShaded(Time,qR_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis tight
hold on
plot(Time,y(:,4)-y(:,2),'b*-','linewidth',2.5)
hold on
plot(Time,y(:,3)-y(:,2),'b:','linewidth',1)
plot(Time,Time*0,'k','LineWidth',.25)
ylim([-1 6])
hold off
filename= fullfile(figurespath, 'Figure3a');
printpdf(f,filename, 'square', 1);

f = figure('Name','Figure 3(b)','NumberTitle','off');
PlotStatesShaded(Time,qCy_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis tight
hold on
plot(Time,y(:,6)-y(:,5),'b:','linewidth',1)
plot(Time,Time*0,'k','LineWidth',.25)
ylim([-1 6])
hold off
filename= fullfile(figurespath,'Figure3b');
printpdf(f,filename, 'square', 1);

f = figure('Name','Figure 3(c)','NumberTitle','off');
PlotStatesShaded(Time,qM_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis tight
hold on
plot(Time,y(:,3)-y(:,2)+y(:,6)-y(:,5),'b:','linewidth',1)
plot(Time,Time*0,'k','LineWidth',.25)
ylim([1 8])
hold off
filename= fullfile(figurespath,'Figure3c');
printpdf(f,filename, 'square', 1);



%% Figure 4
clear Time
%%%%%%%% Charts for the model with Inflation, Tbill, Tbond, AAA, BAA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(tvarTablespath, 'OutMod3forCharts')) % Data from the model with Pi, PiExp, IntShort IntLong, AAA and BAA
%
% The data that enter the model are stored in y
% (in the order Inflation, Inlfation Expectations, Short rate, Expected Short Rate, Long Rate, BAA, AAA)
%
% The quantiles (.025, .16 .5 .84 .975) of the posterior of the trends are stored in
% qR_bar (trend in the real rate), qPi_bar (trend inflation),  qTs_bar(trend of the Term Spread),
% qCy_bar (trend of the Convenience yield) and qM_bar (trend of the discount factor)
% qLiq_bar (trend iquidity), qSafe_bar (trend safety)

tmax = min(find(year(Time)==1998));

f = figure('Name','Figure 4(a)','NumberTitle','off');
bands(Time,qR_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis([-inf inf 0 3.5])
hold off;
filename= fullfile(figurespath,'Figure4a');
printpdf(f, filename, 'square', 1);

f = figure('Name','Figure 4(b)','NumberTitle','off');
bands(Time,qR_bar,-qCy_bar-(-qCy_bar(tmax,3)-qR_bar(tmax,3)))
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis([-inf inf 0 3.5])
hold off;
filename= fullfile(figurespath,'Figure4b');
printpdf(f, filename, 'square', 1);

f = figure('Name','Figure 4(c)','NumberTitle','off');
bands(Time,qR_bar,qM_bar-(qM_bar(tmax,3)-qR_bar(tmax,3)))
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis([-inf inf 0 3.5])
hold off;
filename= fullfile(figurespath,'Figure4c');
printpdf(f, filename, 'square', 1);


%% Figure 5
%%Charts on Trends in Compensation for Safety and Liquidity, and
%%Observables

f = figure('Name','Figure 5(a)','NumberTitle','off');
PlotStatesShaded(Time,qCy_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis([-inf inf 0 3.5])
hold on
plot(Time,y(:,7)-y(:,5),'b:','linewidth',1)
hold off
filename=fullfile(figurespath,'Figure5a');
printpdf(f,filename, 'square', 1);

f = figure('Name','Figure 5(b)','NumberTitle','off');
PlotStatesShaded(Time,qSafe_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis([-inf inf 0 3.5])
hold on
plot(Time,y(:,7)-y(:,6),'b:','linewidth',1)
hold off
filename=fullfile(figurespath,'Figure5b');
printpdf(f,filename, 'square', 1);

f = figure('Name','Figure 5(c)','NumberTitle','off');
PlotStatesShaded(Time,qLiq_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis([-inf inf 0 3.5])
hold on
plot(Time,y(:,6)-y(:,5),'b:','linewidth',1)
hold off
filename= fullfile(figurespath,'Figure5c');
printpdf(f,filename, 'square', 1);


%% Figure 6
clear Time
%%%%%%%% Charts for the model with Consumption%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(tvarTablespath, 'OutMod4forCharts'));

f = figure('Name','Figure 6','NumberTitle','off');
PlotStatesShaded(Time,qDC_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis tight
hold on
temp = filter(ones(4,1)/4,1,y(:,end)); temp(1:4)=NaN;
plot(Time,temp,'b:','linewidth',1)
ylim([-2 6])
plot(Time,Time*0,'k','LineWidth',.25)
hold off
filename= fullfile(figurespath,'Figure6');
printpdf(f,filename, 'square', 0);


%% Figure 7

load(fullfile(tvarTablespath,'OutMod3forCharts')); % Data from the model with Pi, PiExp, IntShort IntLong, AAA and BAA

%%% GZ Spreads

[DATA,TEXT] = xlsread(fullfile(otherTablespath,'gz_spreads.xlsx'),'Sheet1');
TimeGZ = datenum(DATA(:,1))+datenum('12-30-1899','mm-dd-yyyy');

GZspread = DATA(:,2);
excessBP = DATA(:,3);
GZdefault= DATA(:,4);

f = figure('Name','Figure 7(a)','NumberTitle','off');
bands(Time,qCy_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis tight
hold on
plot(Time,y(:,7)-y(:,5),'b:','linewidth',1)
ylim([0.0 5])
set(gca,'YTick',[0:1.0:5])
hold on
yyaxis right
plot(TimeGZ,GZspread)
ylim([-1.2 6.5])
set(gca,'YTick',[-4:1:8])
set(gca,'XLim',[715877 736605]);
hold off
% %%% Switch left and right y-axes (color, tickmarks, and ticklabels).
%     yyaxis left
%     ax_left = gca; y_color1 = ax_left.YColor; y_tick1 = ax_left.YTick; y_tick_label1 = ax_left.YTickLabel; y_lim1 = ax_left.YLim;
%     yyaxis right
%     ax_right = gca; y_color2 = ax_right.YColor; y_tick2 = ax_right.YTick; y_tick_label2 = ax_right.YTickLabel; y_lim2 = ax_right.YLim;
%     yyaxis left
%     set(gca,'YColor',y_color2,'YTick',y_tick2,'YTickLabel',y_tick_label2,'YLim',y_lim2);
%     yyaxis right
%     set(gca,'YColor',y_color1,'YTick',y_tick1,'YTickLabel',y_tick_label1,'YLim',y_lim1);
% %%%
filename = fullfile(figurespath, 'Figure7a');
printpdf(f,filename, 'square', 1);

%%% Industrial spreads

if useBloombergData
  [DATA,TEXT] = xlsread(fullfile(otherTablespath, 'spreads.xlsx'),'monthly');
  IndustrialA_spread = DATA(:,2);
  TimeIndustrial_monthly = datenum(DATA(:,1))+datenum('12-30-1899','mm-dd-yyyy');

  [DATA,TEXT] = xlsread(fullfile(otherTablespath, 'spreads.xlsx'),'quarterly');
  IndustrialBBB_spread = DATA(:,3);
  TimeIndustrial_quarterly = datenum(DATA(:,1))+datenum('12-30-1899','mm-dd-yyyy');
end

f = figure('Name','Figure 7(b)','NumberTitle','off');
bands(Time,qCy_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis tight
hold on
plot(Time,y(:,7)-y(:,5),'b:','linewidth',1)
% datetick('x', 17,'keeplimits', 'keepticks')
ylim([0 5])
set(gca,'YTick',[-1:1.0:5])
hold on
yyaxis right
if useBloombergData
  plot(TimeIndustrial_quarterly,IndustrialBBB_spread)
end
ylim([0 5])
set(gca,'YTick',[0:1:5])
set(gca,'XLim',[715877 736605]);
hold off
filename = fullfile(figurespath, 'Figure7b');
printpdf(f,filename, 'square', 1);

f = figure('Name','Figure 7(c)','NumberTitle','off');
bands(Time,qLiq_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis tight
hold on
plot(Time,y(:,6)-y(:,5),'b:','linewidth',1)
% datetick('x', 17,'keeplimits', 'keepticks')
ylim([0 2])
set(gca,'YTick',[0:0.5:2])
hold on
yyaxis right
if useBloombergData
  plot(TimeIndustrial_monthly,IndustrialA_spread)
end
ylim([-0.3 3.2])
set(gca,'YTick',[-3:0.5:7])
set(gca,'XLim',[715877 736605]);
hold off
filename = fullfile(figurespath, 'Figure7c');
printpdf(f,filename, 'square', 1);


%%% Refcorp

f = figure('Name','Figure 7(d)','NumberTitle','off');
bands(Time,qLiq_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis tight
hold on
plot(Time,y(:,6)-y(:,5),'b:','linewidth',1)
% datetick('x', 17,'keeplimits', 'keepticks')
ylim([0 2])
set(gca,'YTick',[0:0.5:2])
hold on
yyaxis right
plot(TimeRefcorp,Refcorp5Y, 'Color',[203 137 203]/255,'LineStyle','-','LineWidth',0.25)
hold on
plot(TimeRefcorp,Refcorp10Y,'Color',[087 038 087]/255,'LineStyle','-','LineWidth',0.25);
hold on
plot(TimeRefcorp,Refcorp20Y,'Color',[163 073 163]/255,'LineStyle','-','LineWidth',0.25)
ylim([-.5 1])
set(gca,'YTick',[-.5:0.5:1.1])
hold off
filename = fullfile(figurespath, 'Figure7d');
printpdf(f,filename, 'square', 1);

%% Figure 8
load(fullfile(tvarTablespath,'OutModDDforCharts'));

f = figure('Name','Figure 8','NumberTitle','off');
PlotStatesShaded(Time,nan(length(Time),5))
hold on
plot(Time,y(:,end),'b-','linewidth',1)
axis([-inf inf 2 14])
xlim([datenum('01-Jan-1960') datenum('01-Oct-2016')])
hold off
filename= fullfile(figurespath, 'Figure8'); %FiguresModelDD/DD
printpdf(f,filename, 'square', 0);

%% Figure 9
clear Time
%%%%%%%% Charts for the model with Inflation, Tbill, Tbond, AAA, BAA with loose prior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(tvarTablespath, 'OutMod3LooseforCharts')); % Data from the model with Pi, PiExp, IntShort IntLong, AAA and BAA
%
% The data that enter the model are stored in y
% (in the order Inflation, Inlfation Expectations, Short rate, Expected Short Rate, Long Rate, BAA, AAA)
%
% The quantiles (.025, .16 .5 .84 .975) of the posterior of the trends are stored in
% qR_bar (trend in the real rate), qPi_bar (trend inflation),  qTs_bar(trend of the Term Spread),
% qCy_bar (trend of the Convenience yield) and qM_bar (trend of the discount factor)
% qLiq_bar (trend iquidity), qSafe_bar (trend safety)

%%Charts on Trends in Compensation for Safety and Liquidity, and
%%Observables

%%Charts of Trends and Observables, Convenience Yield Model

f = figure('Name','Figure 9(a)','NumberTitle','off');
PlotStatesShaded(Time,qR_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis tight
hold on
plot(Time,y(:,4)-y(:,2),'b*-','linewidth',2.5)
hold on
plot(Time,y(:,3)-y(:,2),'b:','linewidth',1)
ylim([-1 6])
plot(Time,Time*0,'k','LineWidth',.25)
hold off
filename=fullfile(figurespath, 'Figure9a');
printpdf(f,filename, 'square', 1);

f = figure('Name','Figure 9(b)','NumberTitle','off');
PlotStatesShaded(Time,qCy_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis tight
hold on
plot(Time,y(:,7)-y(:,5),'b:','linewidth',1)
ylim([-1 6])
plot(Time,Time*0,'k','LineWidth',.25)
hold off
filename=fullfile(figurespath, 'Figure9b');
printpdf(f,filename, 'square', 1);

f = figure('Name','Figure 9(c)','NumberTitle','off');
PlotStatesShaded(Time,qM_bar)
set(gca,'XTick', Time(1:40:end),'XMinorTick','on')
box on; axis tight
hold on
plot(Time,y(:,3)-y(:,2)+y(:,7)-y(:,5),'b:','linewidth',1)
ylim([1 8])
plot(Time,Time*0,'k','LineWidth',.25)
hold off
filename= fullfile(figurespath, 'Figure9c');
printpdf(f,filename, 'square', 1);

%% Figure 10: Drivers of TVAR r*
load(fullfile(tvarTablespath, 'OutMod4forCharts'))
shockdec = [TVAR_cybar TVAR_gbar TVAR_bbar];
shockcats     = {[1]; [2]; [3]};
shockdecColors = {red; blue; 'LightGray'};
% shockdecColors = {red; blue; light_yellow};
shockcatNames = {};
trend = 2.5329;

qR_bar_bands = [qR_bar(:, 1:2) qR_bar(:, 4:end)] - trend;
fig10 = figure('Name','Figure 10','NumberTitle','off');
% hold on;
% bands(Time,qR_bar_bands);
hold on;
[barPos, tmp, tmp] = plotStackedShockdec(TVAR_time, shockdec, TVAR_rbar, shockcats, ...
shockdecColors, shockcatNames, 'fig', fig10);
legend(barPos, {'Convenience Yield'; 'Growth ($\overline{g}$)'; 'Other ($\overline{\beta}$)'}, ...
    'interpreter', 'latex', 'location', 'southoutside', 'orientation', 'horizontal');
legend boxoff
hold on
% label y-axis centered at level of r*
if adjustLevel

  % DSGE r* trend value
  accuracy = 0.5;
  ax = gca;
  ymax = 3.1;
  ymin = 1.0;
  ax1 = addTrendToShockdec(ax, trend, accuracy, 'ymax', ymax, 'ymin', ymin, ...
  'forceUpper', 1);
end

% qR_bar_bands = [qR_bar(:, 1:2) qR_bar(:, 4:end)] - trend;
% bands(Time,qR_bar_bands);

hold off
printpdf(fig10, [figurespath, 'Figure10'], 'square', 0, 'fontsize', 12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DSGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Figure 11
% Panel A: DSGE 20yr fwd r and rstar with trendy VAR
fig11a = figure('Name','Figure 11(a)','NumberTitle','off');
l1 = PlotStatesShadedv3_oneband(time,rstarTrendyVAR(:,2:4)); box on;
hold on;
l2 = PlotStatesShadedv3_oneband(time,DSGE_rstar20fwd(:,2:4),[0 0 1], 0.2);
hold off;
leg = legend([l1,l2], 'VAR', 'DSGE', 'location', 'SouthWest','interpreter','latex');
set(leg,'interpreter', 'latex')
ylim([0,3.5]);
printpdf(fig11a,fullfile(figurespath, 'Figure11a'), 'square', 1);

% Panel B: DSGE 30yr fwd r and rstar with TVAR
fig11b = figure('Name','Figure 11(b)','NumberTitle','off');
l1 = PlotStatesShadedv3_oneband(time,rstarTrendyVAR(:,2:4));
hold on;
l2 = PlotStatesShadedv3_oneband(time,DSGE_rstar30fwd(:,2:4),[0 0 1], 0.2);
hold off;
leg = legend([l1,l2], 'VAR', 'DSGE', 'location', 'SouthWest','interpreter','latex');
set(leg,'interpreter', 'latex')
ylim([0,3.5]);
printpdf(fig11b,[figurespath, 'Figure11b'], 'square', 1);

%% Figure 12
% Panel A: DSGE 5yr fwd r and rstar
fig12a = figure('Name','Figure 12(a)','NumberTitle','off');
l1 = PlotStatesShadedv3_oneband(time,DSGE_r5fwd(:,2:4),[1 0 0]);
hold on;
l2 = PlotStatesShadedv3_oneband(time,DSGE_rstar5fwd(:,2:4),[0 0 1], 0.2);
box on;
hold off;
leg = legend([l1,l2], 'Actual', 'Natural', 'location', 'SouthWest','interpreter','latex');
set(leg,'interpreter', 'latex')
ylim([-1,5]);
printpdf(fig12a,[figurespath, 'Figure12a'], 'square', 1);

% Panel B: DSGE 10yr fwd r and rstar
fig12b = figure('Name','Figure 12(b)','NumberTitle','off');
l1 = PlotStatesShadedv3_oneband(time,DSGE_r10fwd(:,2:4),[1 0 0]);
hold on;
l2 = PlotStatesShadedv3_oneband(time,DSGE_rstar10fwd(:,2:4),[0 0 1], 0.2);
box on;
hold off;
leg = legend([l1,l2], 'Actual', 'Natural', 'location', 'SouthWest','interpreter','latex');
set(leg,'interpreter', 'latex')
ylim([-1,5]);
printpdf(fig12b,[figurespath, 'Figure12b'], 'square', 1);


%% Figure 13: DSGE 5-year fwd r* vs LW r*
fig13 = figure('Name','Figure 13','NumberTitle','off');
hold on;
l3 = PlotStatesShadedv3_oneband(time,DSGE_rstar5fwd(:,2:4),[0 0 1], 0.2);
hold on;
l5 = plot(time,rstarLW,'Color',green,'LineWidth',2.0);
ylim([-2 7])
set(gca, 'YTick', -2:1:7)
hold off;
leg = legend([l3,l5], 'DSGE', 'Laubach-Williams');
set(leg,'interpreter', 'latex')
ylim([-1,7]);
printpdf(fig13,[figurespath, 'Figure13'], 'square', 0);

%% Figure 14: 30-year fwd r* and component attributable to convenience yield
%shockcats     = {[b_liqtil_sh b_liqp_sh b_safetil_sh b_safep_sh]; ...
%                 [mu_sh z_sh zp_sh g_sh sigw_sh]};
%shockdecColors = {[191, 68, 68]./256; light_yellow};
%shockcatNames  = {'Convenience Yield', 'Other'};

shockcats   = {{'b_liqtil_sh', 'b_liqp_sh', 'b_safetil_sh', 'b_safep_sh'}; ...
                 {'sigma_omega_sh'}; {'z_sh', 'zp_sh'}; {'mu_sh', 'g_sh'}};
shockdecColors = {red; green; blue; 'LightGray'};
shockcatNames = {'Convenience Yield'; 'Risk'; 'Productivity'; 'Other'};

filename = [dsgeTablespath, '/shockdec_Forward30YearRealNaturalRate', suffix];
trendFilename = [dsgeTablespath, '/trend_Forward30YearRealNaturalRate', suffix];

fig14 = figure('Name','Figure 14','NumberTitle','off');
[barPos, tmp, tmp] = prepareStackedShockdec(filename,  trendFilename, shockdecStart,...
    shockdecEnd, shockcats, shockcatNames, shockdecColors, 'fig', fig14);
legend(barPos, shockcatNames, 'interpreter', 'latex', ...
    'location', 'southoutside', 'orientation', 'horizontal');
legend boxoff
if adjustLevel
  % DSGE r* trend value
  trend = 1.92886;
  accuracy = 0.5;
  ax = gca;
  ymax = 3.0;
  ymin = 1.0;
  ax1 = addTrendToShockdec(ax, trend, accuracy, 'ymax', ymax, 'ymin', ymin);
end
printpdf(fig14, [figurespath, 'Figure14'], 'square', 0, 'fontsize', 12);


%% Figure 15: DSGE ST r and r*
fig15 = figure('Name','Figure 15','NumberTitle','off');
l1 = PlotStatesShadedv3_oneband(time,DSGE_rST(:,2:4),[1 0 0]);
hold on;
l2 = PlotStatesShadedv3_oneband(time,DSGE_rstarST(:,2:4),[0 0 1], 0.2);
ylim([-9,12]);
set(gca, 'YTick', [-9:3:12]);
box on;
hold off;
leg = legend([l1,l2], 'r', 'r*');
set(leg,'interpreter', 'latex');
printpdf(fig15,[figurespath, 'Figure15'], 'square', 0);


%% Figure 16: Shock decomposition of r*
shockcats     = {{'b_liqtil_sh', 'b_liqp_sh', 'b_safetil_sh', 'b_safep_sh'}; ...
                 {'sigma_omega_sh'}; {'z_sh', 'zp_sh'}; {'mu_sh', 'g_sh'}};
shockdecColors = {red; green; blue; 'LightGray'};
shockcatNames = {'Convenience Yield'; 'Risk'; 'Productivity'; 'Other'};
filename = [dsgeTablespath, '/shockdec_RealNaturalRate', suffix];
trendFilename = [dsgeTablespath, '/trend_RealNaturalRate', suffix];
fig16 = figure('Name','Figure 16','NumberTitle','off');
[barPos, tmp, tmp] = prepareStackedShockdec(filename, trendFilename, shockdecStart,...
    shockdecEnd, shockcats, shockcatNames, shockdecColors, 'fig', fig16);
ylim([-10,12]);
set(gca, 'YTick', [-10:2:12]);

legend(barPos, {'Convenience Yield'; 'Risk'; 'Productivity'; 'Other'}, 'interpreter', 'latex', ...
    'location', 'southoutside', 'orientation', 'horizontal');
legend boxoff
printpdf(fig16, [figurespath, 'Figure16'], 'square', 0, 'fontsize', 12);
