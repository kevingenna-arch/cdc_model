% Original Guass code by Burkhard Heer, March 13, 2015;
% Translated into Matlab code, Sijmen Duineveld, February 23, 2018

% If  Econometrics toolbox unavialable: set "addpath tools" and replace 'hpfilter' with 'hpfltr' (both have cyclical component as second output argument)

clear all;
close all;
clc;

dbstop if error;

%%
data = dlmread('FRED_data_2_24okt2015.txt');
data1 = dlmread('FRED_data_2_28okt2015.txt');

period  = 1947 + 0.25*[0:273];
period1 = 1959.25 + 0.25*[0:225];

price=data1(:,11);
qinflation=(price(2:226)-price(1:225))./price(1:225);	% quarterly inflation
inflation=(1+qinflation).^4-1;	% annual inflation
riskfree=data1(2:226,9);		% 3-Month Treasury Bill: Secondary Market Rate, Percent, Quarterly, Not Seasonally Adjusted	
riskfree=riskfree/100;			
realriskfree=riskfree-inflation;
fprintf('mean real risk-free rate: %.4f \n',mean(realriskfree));


returnequity1=data1(2:226,10); % Standard & Poors Total Return, yield, Percent, Quarterly, Not seasonally Adjusted, own calculation
returnequity1=returnequity1/100;
%returnequity1=(1+returnequity1)^4-1;	% annualized
realreturnequity=returnequity1-qinflation;
fprintf('mean real equity return: %.4f \n\n',mean(realreturnequity));


%%
gdp2    = data1(2:226,2);		%Real Gross Domestic Product, Billions of Chained 2009 Dollars, Quarterly, Seasonally Adjusted Annual Rate	
hours2  = data1(2:226,3);	%Nonfarm Business Sector: Hours of All Persons, Index 2009=100, Quarterly, Seasonally Adjusted

loggdp2         = log(gdp2);
[~,gdpcycle2]       = hpfilter(loggdp2,1600);


[~,hourscycle2]    = hpfilter(log(hours2),1600);
[~,realrfcycle1]    = hpfilter(realriskfree,1600);
[~,realrecycle1]    = hpfilter(realreturnequity,1600);

fprintf('standard deviations, subperiod 1959.3-2015.2 \n');
fprintf('gdp: %.4f \n',std(gdpcycle2));
fprintf('hours: %.4f \n', std(hourscycle2));
fprintf('real risk-free rate: %.4f \n', std(realrfcycle1));
fprintf('real equity return: %.4f \n\n', std(realrecycle1));


fprintf('correlations, subperiod 1959.3-2015.2, GDP~Hours~rf~re: \n');
RR = corrcoef([gdpcycle2,hourscycle2,realrfcycle1,realrecycle1]);
%RR
fprintf('%+.4f %+.4f %+.4f %+.4f \n', RR);
fprintf('\n');


%%
gdp         = data(:,2);		%Real Gross Domestic Product, Billions of Chained 2009 Dollars, Quarterly, Seasonally Adjusted Annual Rate	
gdp1        = data(37:274,2);
hours       = data(:,3);	%Nonfarm Business Sector: Hours of All Persons, Index 2009=100, Quarterly, Seasonally Adjusted
hours1      = data(37:274,3);
cp          = data(:,5);		%	Real Personal Consumption Expenditures, Billions of Chained 2009 Dollars, 
					%  Quarterly, Seasonally Adjusted Annual Rate	
cp1         = data(37:274,5);
gc          = data(:,6);		% Government Consumption Expenditures, Billions of Chained 2009 Dollars, 
					% Quarterly, Seasonally Adjusted Annual Rate, own calculation		
gc1         = data(37:274,6);

wage        = data(:,7);		% Nonfarm Business Sector: Compensation Per Hour, Index 2009=100, Quarterly, Seasonally Adjusted		
wage1       = data(37:274,7);

investment  = data(:,8);	% Private Nonresidential Fixed Investment, Billions of Dollars, 
						% Quarterly, Seasonally Adjusted Annual Rate
investment1 = data(37:274,8);

riskfree    = data(:,9);		% 3-Month Treasury Bill: Secondary Market Rate, Percent, Quarterly, Not Seasonally Adjusted	
riskfree1   = data(37:274,9);

returnequiy     = data(:,10); % Standard & Poors Total Return, yield, Percent, Quarterly, Not seasonally Adjusted, own calculation
returnequity1   = data(37:274,10);

%%
loggdp                  = log(gdp);
[gdptrend,gdpcycle]     = hpfilter(loggdp,1600);

[hourstrend,hourscycle] = hpfilter(log(hours),1600);
[cptrend,cpcycle]       = hpfilter(log(cp),1600);
[gctrend,gccycle]       = hpfilter(log(gc),1600);


fprintf('Statistics 1948-2014 \n');
fprintf('Correlations: GDP, G, L, Cp: \n');
fprintf('%+.4f %+.4f %+.4f %+.4f \n',corrcoef([gdpcycle,gccycle,hourscycle,cpcycle]));
fprintf('\n');
fprintf('standard deviations %.4f %.4f %.4f %.4f \n',std(gdpcycle),std(gccycle),std(hourscycle),std(cpcycle));
fprintf('\n');


%%
loggdp1                     = log(gdp1);
[gdptrend1,gdpcycle1]       = hpfilter(loggdp1,1600);
[hourstrend1,hourscycle1]   = hpfilter(log(hours1),1600);
[cptrend1,cpcycle1]         = hpfilter(log(cp1),1600);
[gctrend1,gccycle1]         = hpfilter(log(gc1),1600);
[wagetrend1,wagecycle1]     = hpfilter(log(wage1),1600);
[itrend1,icycle1]           = hpfilter(log(investment1),1600);
[rftrend1,rfcycle1]         = hpfilter(riskfree1,1600);
[retrend1,recycle1]         = hpfilter(returnequity1,1600);

fprintf('statistics 1953-2014: \n');
fprintf('Correlations: GDP, G, L, Cp, w, I, rf, re: \n');

st = 1;%starting period
en = length(gdpcycle1);%last period
fprintf('%+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f \n',corrcoef([gdpcycle1(st:en),gccycle1(st:en),hourscycle1(st:en),cpcycle1(st:en),wagecycle1(st:en),icycle1(st:en),rfcycle1(st:en),recycle1(st:en)]));
fprintf('\n');
fprintf('covariances: \n');
fprintf('%+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f %+.4f \n',cov([gdpcycle1(st:en),gccycle1(st:en),hourscycle1(st:en),cpcycle1(st:en),wagecycle1(st:en),icycle1(st:en),rfcycle1(st:en),recycle1(st:en)]));
clear st en;

%% Figure 2.11
LnWd = 1;

figure
plot(period,loggdp,'--','LineWidth',LnWd);
hold all;
plot(period,gdptrend,'-','LineWidth',LnWd);
xlabel('Period');
ylabel('GDP');
legend({'Original series','Trend'})
xlim([period(1) period(end)]);
title('GDP');


%% Figure 2.13
figure
plot(period,gdpcycle,'-','LineWidth',LnWd);
hold all;
plot(period,hourscycle,'--','LineWidth',LnWd)
xlabel('Period');
ylabel('Cyclical Component');
legend({'GDP','Hours'});
xlim([period(1) period(end)]);
title('GDP and Hours');


%% Cyclical component Government Consumption
figure
plot(period,gdpcycle,'-','LineWidth',LnWd);
hold all;
plot(period,gccycle,'--','LineWidth',LnWd)
xlabel('Period');
ylabel('Cyclical Component');
legend({'GDP','Gov. Cons.'});
xlim([period(1) period(end)]);
title('Gov. Cons.');


%% Original and trend in Government Consumption
figure
plot(period,log(gc),'--','LineWidth',LnWd);
hold all;
plot(period,gctrend,'-','LineWidth',LnWd);
xlabel('Period');
ylabel('Gov. Cons.');
legend({'Original series','Trend'})
xlim([period(1) period(end)]);
title('Gov. Cons.');

fprintf('\n');
for ii = 0:12    
    zz = corrcoef(gdpcycle1(1:226-ii),gccycle1(1+ii:226));
    fprintf('Correlation GDP(t- %d) & G(t): %.4f \n',ii,zz(1,2));
end


%% Figure 2.12
figure
nperiods = 52;
period_tmp = 1947 + 0.25*[0:nperiods-1];
plot(period_tmp,loggdp(1:nperiods),'--','LineWidth',LnWd);
hold all;
plot(period_tmp,gdptrend(1:nperiods),'-','LineWidth',LnWd);
xlabel('Period');
ylabel('GDP');
legend({'Original series','Trend'})
xlim([period_tmp(1) period_tmp(end)]);
title('GDP');














