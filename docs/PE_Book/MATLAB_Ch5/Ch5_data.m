% ----------------------------  Ch4_data.g --------------------------------
%
% provides some of the statistics of Chapter 5
%
% author: Burkhard Heer
%
% this version: June 19, 2018
%
% -------------------------------------------------------------------------



clear; clc;



us_tax=xlsread('Fig_5_1_data_bh.xlsx',3,'E42:DC42');
% Deletion of the rows with no entry
TF = (us_tax>0);
TF1 = (TF<1);
us_tax(TF1) = [];

can_tax=xlsread('Fig_5_1_data_bh.xlsx',3,'E11:DC11');
% Deletion of the rows with no entry
TF = (can_tax>0);
TF1 = (TF<1);
can_tax(TF1) = [];

fra_tax=xlsread('Fig_5_1_data_bh.xlsx',3,'E17:DC17');
% Deletion of the rows with no entry
TF = (fra_tax>0);
TF1 = (TF<1);
fra_tax(TF1) = [];

ger_tax=xlsread('Fig_5_1_data_bh.xlsx',3,'E18:DC18');
% Deletion of the rows with no entry
TF = (ger_tax>0);
TF1 = (TF<1);
ger_tax(TF1) = [];

ita_tax=xlsread('Fig_5_1_data_bh.xlsx',3,'E24:DC24');
% Deletion of the rows with no entry
TF = (ita_tax>0);
TF1 = (TF<1);
ita_tax(TF1) = [];

jap_tax=xlsread('Fig_5_1_data_bh.xlsx',3,'E25:DC25');
% Deletion of the rows with no entry
TF = (jap_tax>0);
TF1 = (TF<1);
jap_tax(TF1) = [];
jap_tax(52) = NaN;

spa_tax=xlsread('Fig_5_1_data_bh.xlsx',3,'E37:DC37');
% Deletion of the rows with no entry
TF = (spa_tax>0);
TF1 = (TF<1);
spa_tax(TF1) = [];

uk_tax=xlsread('Fig_5_1_data_bh.xlsx',3,'E41:DC41');
% Deletion of the rows with no entry
TF = (uk_tax>0);
TF1 = (TF<1);
uk_tax(TF1) = [];

nint=52;
period_tax=linspace(1965,1965+nint-1,nint);


figure
plot(period_tax,can_tax,period_tax,fra_tax,period_tax,ger_tax,period_tax,ita_tax,period_tax,jap_tax,period_tax,spa_tax,period_tax,uk_tax,period_tax,us_tax);
legend('Canada','France','Germany','Italy','Japan','Spain','UK','US');
xlabel('Year');
ylabel('Tax revenue');
pause;

data=xlsread('Ch5_data_matlab','A4:I251');

period_tax=linspace(1947.0,2008.75,248);

tauk=data(:,8);
logtauk=log(tauk);
tauktrend=hpfilter(logtauk,1600);
taukcycle=logtauk-tauktrend;

y=taukcycle(42:248);
x=taukcycle(41:247);
x1=ones(207,1);
x2=horzcat(x1,x);
[btauk,binttauk,rtauk,rinttauk,statstauk] = regress(y,x2);
disp('regression coefficients AR(2) tau_k');
btauk
pause;

taul=data(:,9);
logtaul=log(taul);
taultrend=hpfilter(logtaul,1600);
taulcycle=logtaul-taultrend;

y=taukcycle(43:248);
x=taukcycle(42:247);
xp=taukcycle(41:246);
x1=ones(206,1);
x2=horzcat(x1,x,xp);
[btaul,binttaul,rtaul,rinttaul,statstaul] = regress(y,x2);
disp('regression coefficients AR(2) tau_l');
btaul
pause;


figure
plot(period_tax,tauk,period_tax,taul); hold on
xlabel('Period')
legend('Capital income tax \tau_k','Labor income tax \tau_l')



gdp=data(:,2);
loggdp=log(gdp);
gdptrend=hpfilter(loggdp,1600);
gdpcycle=loggdp-gdptrend;

hours=data(:,3);
loghours=log(hours);
hourstrend=hpfilter(loghours,1600);
hourscycle=loghours-hourstrend;


cp=data(:,5);
logcp=log(cp);
cptrend=hpfilter(logcp,1600);
cpcycle=logcp-cptrend;


gc=data(:,6);
loggc=log(gc);
gctrend=hpfilter(loggc,1600);
gccycle=loggc-gctrend;


disp('statistics 1948-2014:');
disp('Correlations: GDP, G, L, Cp, tauk, taul: ');
x=horzcat(gdpcycle,gccycle,hourscycle,cpcycle,taukcycle,taulcycle);
corrcoef(x)
pause;
disp('statistics 1948-2014:');
disp('standard deviations');
std(x)
pause;


tauk=data(41:248,8);
logtauk=log(tauk);
tauktrend=hpfilter(logtauk,1600);
taukcycle=logtauk-tauktrend;

taul=data(41:248,9);
logtaul=log(taul);
taultrend=hpfilter(logtaul,1600);
taulcycle=logtaul-taultrend;

gdp=data(41:248,2);
loggdp=log(gdp);
gdptrend=hpfilter(loggdp,1600);
gdpcycle=loggdp-gdptrend;

hours=data(41:248,3);
loghours=log(hours);
hourstrend=hpfilter(loghours,1600);
hourscycle=loghours-hourstrend;


cp=data(41:248,5);
logcp=log(cp);
cptrend=hpfilter(logcp,1600);
cpcycle=logcp-cptrend;


gc=data(41:248,6);
loggc=log(gc);
gctrend=hpfilter(loggc,1600);
gccycle=loggc-gctrend;


disp('statistics 1956-2014:');
disp('Correlations: GDP, G, L, Cp, tauk, taul: ');
x=horzcat(gdpcycle,gccycle,hourscycle,cpcycle,taukcycle,taulcycle);
corrcoef(x)
pause;
disp('statistics 1956-2014:');
disp('standard deviations');
std(x)

