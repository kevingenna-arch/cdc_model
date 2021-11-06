% ----------------------------  Ch4_data.g --------------------------------
%
% provides some of the statistics of Chapter 4
%
% author: Burkhard Heer
%
% this version: June 6, 2018
%
% -------------------------------------------------------------------------



clear; clc;


% loading data for interest rate 2001-2017
can_GY=xlsread('Fig_4_1_data.xlsx','B5:AN5');
fra_GY=xlsread('Fig_4_1_data.xlsx','B6:AN6');
ger_GY=xlsread('Fig_4_1_data.xlsx','B7:AN7');
uk_GY=xlsread('Fig_4_1_data.xlsx','B8:AN8');
us_GY=xlsread('Fig_4_1_data.xlsx','B9:AN9');
can_GY=can_GY';
fra_GY=fra_GY';
ger_GY=ger_GY';
uk_GY=uk_GY';
us_GY=us_GY';


nint=size(can_GY,1);
periods=linspace(1980,1980+nint-1,nint);

figure
plot(periods,can_GY,periods,fra_GY,periods,ger_GY,periods,uk_GY,periods,us_GY); hold on
xlabel('Year')
ylabel('Government expenditures')
legend('Canada','France','Germany','UK','US')

%pause;

dataA=importdata('FRED_data_2_24okt2015.txt');  % data 1947.1-2015.2, excluding inflation
dataB=importdata('FRED_data_2_28okt2015.txt');  % data 1959.1-2015.2, including PCE index for prices

t1=1947+273*0.25;
period=linspace(1947,t1,274);
t2=1959.25+225*0.25;
periods1=linspace(1959.25,t2,226);

gdp=dataA(:,2); % Real Gross Domestic Product, Billions of Chained 2009 Dollars, Quarterly, Seasonally Adjusted Annual Rate	
loggdp=log(gdp);
gdptrend=hpfilter(loggdp,1600);
gdpcycle=loggdp-gdptrend;

gc=dataA(:,6);  % Government Consumption Expenditures, Billions of Chained 2009 Dollars, 
% Quarterly, Seasonally Adjusted Annual Rate, own calculation		
loggc=log(gc);
gctrend=hpfilter(loggc,1600);
gccycle=loggc-gctrend;	


figure
plot(period,gc,period,gctrend); hold on
xlabel('Year')
ylabel('Government consumption')
legend('Original series','Trend');

%pause;
gdp1=dataA(37:274,2);
loggdp1=log(gdp1);
gdptrend1=hpfilter(loggdp1,1600);
gdpcycle1=loggdp1-gdptrend1;

gc1=dataA(37:274,6);
loggc1=log(gc1);
gctrend1=hpfilter(loggc1,1600);
gccycle1=loggc1-gctrend1;

figure
plot(period,gdpcycle,period,gccycle); hold on
xlabel('Year')
ylabel('Cyclical component')
legend('GDP','Government consumption');

disp('Correlation coefficients');

corgy=zeros(8,1);

for i=0:7
    temp=corrcoef(gdpcycle1(1:226-i),gccycle1(1+i:226));
    corgy(i+1)=temp(1,2);
end

temp1=1:1:8;
corrmatrix=[temp1' corgy]   



