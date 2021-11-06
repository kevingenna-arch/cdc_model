% ----------------------------  Ch7_data.m --------------------------------
%
% provides some of the statistics of Chapter 7
%
% author: Burkhard Heer
%
% this version: June 20, 2018
%
% -------------------------------------------------------------------------




clear; clc;


deficitperiod=linspace(2000,2016,17);
deficitperiod=deficitperiod';
ausdeficit=xlsread('gov_deficit_OECD.xls','B1:B16');
x0=NaN;
ausdeficit=vertcat(ausdeficit,x0);
candeficit=xlsread('gov_deficit_OECD.xls','B51:B67');
fradeficit=xlsread('gov_deficit_OECD.xls','B119:B135');
gerdeficit=xlsread('gov_deficit_OECD.xls','B136:B152');
gredeficit=xlsread('gov_deficit_OECD.xls','B153:B169');
itadeficit=xlsread('gov_deficit_OECD.xls','B221:B237');
japdeficit=xlsread('gov_deficit_OECD.xls','B238:B248');
x1=NaN * ones(5,1);
japdeficit=vertcat(x1,japdeficit,x0);
spadeficit=xlsread('gov_deficit_OECD.xls','B381:B397');
ukdeficit=xlsread('gov_deficit_OECD.xls','B439:B455');
usdeficit=xlsread('gov_deficit_OECD.xls','B456:B471');
usdeficit=vertcat(usdeficit,x0);

disp('mean Italian deficit 2001-2016): ');
mean(itadeficit(2:17))
pause;

figure
plot(deficitperiod,ausdeficit,deficitperiod,candeficit,deficitperiod,fradeficit,deficitperiod,gerdeficit,deficitperiod,gredeficit);
xlabel('Year');
ylabel('Government budget');
legend('Australia','Canada','France','Germany','Greece');
pause;

plot(deficitperiod,itadeficit,deficitperiod,japdeficit,deficitperiod,spadeficit,deficitperiod,ukdeficit,deficitperiod,usdeficit);
legend('Italy','Japan','Spain','UK','US');
xlabel('Year');
ylabel('Government budget');
pause;


debt_gdp=xlsread('Debt_GDP_US_GFDEGDQ188S.xls','B12:B217');
period=linspace(1966,2017.25,206);
figure
plot(period,debt_gdp);
xlabel('Year');
ylabel('US Debt-GDP ratio');

pause;


periodgdp=linspace(1961,2016,56);
italy_gdp_growth=xlsread('world_bank_econ_growth_data.xls','F119:BI119');
italy_inflation=xlsread('world_bank_inflation_rate.xls','F119:BI119');
figure
plot(periodgdp,italy_gdp_growth);
xlabel('Year');
ylabel('Real GDP growth rate Italy');
pause;

disp('mean GDP growth Italy, 2001-2016: ');
mean(italy_gdp_growth(40:56))
disp('mean inflation Italy, 2001-2016: ');
mean(italy_inflation(40:56))

pause;

periodusgdp=linspace(1948,2017.5,279);
us_gdp_growth=xlsread('growth_real_gdp_US.xls','C16:C294');
figure 
plot(periodusgdp,us_gdp_growth);
xlabel('Year');
ylabel('Real GDP growth rate US');
pause;

periodusinflation=linspace(1960,2016,57);
us_inflation=xlsread('inflation_US.xls','B12:B68');
plot(periodusinflation,us_inflation);
xlabel('Year');
ylabel('Inflation rate US');
pause;


disp('mean inflation rate US, 2000-2016: ');
mean(us_inflation(41:57))
disp('mean growth rate US, 2000-2016: ');
mean(us_gdp_growth(229:276))
pause;

aus_ndebt=xlsread('IMF_netdebt1.xls','B4:AC4');
%aus_gdebt=aus_gdebt';
can_ndebt=xlsread('IMF_netdebt1.xls','B5:AC5');
%can_gdebt=can_gdebt';
fra_ndebt=xlsread('IMF_netdebt1.xls','B6:AC6');
ger_ndebt=xlsread('IMF_netdebt1.xls','B7:AC7');
x0=NaN*ones(1,6);
ger_ndebt=horzcat(x0,ger_ndebt);
ita_ndebt=xlsread('IMF_netdebt1.xls','B8:AC8');
%ita_gdebt=ita_gdebt';
jap_ndebt=xlsread('IMF_netdebt1.xls','B9:AC9');
%jap_gdebt=jap_gdebt';
spa_ndebt=xlsread('IMF_netdebt1.xls','B10:AC10');
%spa_gdebt=spa_gdebt';
gbr_ndebt=xlsread('IMF_netdebt1.xls','B11:AC11');
%gbr_gdebt=gbr_gdebt';
us_ndebt=xlsread('IMF_netdebt1.xls','B12:AC12');
x0=NaN*ones(1,12);
us_ndebt=horzcat(x0,us_ndebt);
%us_gdebt=us_gdebt';


perioddebt=linspace(1989,2006,28);


figure
plot(perioddebt,aus_ndebt,perioddebt,can_ndebt,perioddebt,fra_ndebt, perioddebt,ger_ndebt, perioddebt,ita_ndebt);
xlabel('Year');
ylabel('Net debt-GDP ratio');
legend('Australia','Canada','France','Germany','Italy');
pause;

figure
plot(perioddebt,jap_ndebt,perioddebt,spa_ndebt ,perioddebt,gbr_ndebt,perioddebt,us_ndebt); 
xlabel('Year');
ylabel('Net debt-GDP ratio');
legend('Japan','Spain','UK','US');
pause;



aus_gdebt=xlsread('IMF_grossdebt1.xls','B4:AC4');
%aus_gdebt=aus_gdebt';
can_gdebt=xlsread('IMF_grossdebt1.xls','B5:AC5');
%can_gdebt=can_gdebt';
chn_gdebt=xlsread('IMF_grossdebt1.xls','B6:AC6');
x0=NaN*ones(1,6);
chn_gdebt=horzcat(x0,chn_gdebt);
%chn_gdebt=chn_gdebt';
fra_gdebt=xlsread('IMF_grossdebt1.xls','B7:AC7');
%fra_gdebt=fra_gdebt';
ger_gdebt=xlsread('IMF_grossdebt1.xls','B8:AC8');
x0=NaN*ones(1,2);
ger_gdebt=horzcat(x0,ger_gdebt);
%ger_gdebt=ger_gdebt';
gre_gdebt=xlsread('IMF_grossdebt1.xls','B9:AC9');
%gre_gdebt=gre_gdebt';
ita_gdebt=xlsread('IMF_grossdebt1.xls','B10:AC10');
%ita_gdebt=ita_gdebt';
jap_gdebt=xlsread('IMF_grossdebt1.xls','B11:AC11');
%jap_gdebt=jap_gdebt';
spa_gdebt=xlsread('IMF_grossdebt1.xls','B12:AC12');
%spa_gdebt=spa_gdebt';
gbr_gdebt=xlsread('IMF_grossdebt1.xls','B13:AC13');
%gbr_gdebt=gbr_gdebt';
us_gdebt=xlsread('IMF_grossdebt1.xls','B14:AC14');
x0=NaN*ones(1,12);
us_gdebt=horzcat(x0,us_gdebt);
%us_gdebt=us_gdebt';

perioddebt=linspace(1989,2006,28);


figure
plot(perioddebt,aus_gdebt,perioddebt, can_gdebt,perioddebt, chn_gdebt ,perioddebt,fra_gdebt, perioddebt,ger_gdebt, perioddebt,gre_gdebt);
xlabel('Year');
ylabel('Gross debt-GDP ratio');
legend('Australia','Canada','China','France','Germany','Greece');
pause;

figure
plot(perioddebt,ita_gdebt,perioddebt,jap_gdebt,perioddebt,spa_gdebt ,perioddebt,gbr_gdebt,perioddebt,us_gdebt); xlabel('Year');
ylabel('Gross debt-GDP ratio');
legend('Italy','Japan','Spain','UK','US');
pause;

% loading data for interest rate 2001-2017
fra_i=xlsread('10_year_government_bond_yields.xls','B511:B712');
ger_i=xlsread('10_year_government_bond_yields.xls','C511:C712');
gre_i=xlsread('10_year_government_bond_yields.xls','D511:D712');
ire_i=xlsread('10_year_government_bond_yields.xls','E511:e712');
ita_i=xlsread('10_year_government_bond_yields.xls','F511:F712');
por_i=xlsread('10_year_government_bond_yields.xls','G511:G712');
spa_i=xlsread('10_year_government_bond_yields.xls','H511:H712');
nint=size(fra_i);
nint=nint(1);

period_i=linspace(2001,2017.75,nint);


figure
plot(period_i,fra_i,period_i,gre_i,period_i,ire_i,period_i,ita_i,period_i,por_i,period_i,spa_i);
legend('France','Greece','Ireland','Italy','Portugal','Spain');
xlabel('Period');
ylabel('10-year government bond yield');
pause;
