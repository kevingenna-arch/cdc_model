function problem_P5_4d();
% 
%	Problem 5.4d  on the Laffer Curve
%   from my book on "Public Economics"
%   Author: Burkhard Heer
%   Last Change: Jan 3, 2020


clc;
close all;


% Parameterization 
%
par.alpha=0.38;
par.rbar=1.04;      % annual real interest rate on bonds
par.nbar=0.25;		% steady-state labor supply
par.kappa=3.46;     % utility parameter: labor supply
par.psi=1.02;          % annual growth factor
par.varphi=1;		% Frisch labor supply elasticity
par.sigma=2.0;		% 1/IES
par.taul=0.28;      % labor income tax
par.taulbar = 0.28; 
par.tauk=0.36;      % capital income tax
par.tauc=0.05;      % consumption tax
par.by = 0.63;      % debt-GDP ratio
par.gy = 0.18;      % G/Y
par.delta=0.07;     % depreciation
psi0 = par.psi;
rbar = par.rbar;
alpha = par.alpha;
tauk = par.tauk;
tauc = par.tauc;
taul = par.taul;
delta = par.delta;
varphi = par.varphi;
nbar = par.nbar;
gybar = par.gy;
bybar = par.by;
sigma = par.sigma;

% calibration of the steady state 
beta0=psi0^(sigma)/rbar;
par.beta = beta0;
%ky=(rbar-1)/(alpha*(1-tauk))+delta/alpha;
ky=(rbar-1)/(alpha*(1-tauk))+delta/(alpha*(1-tauk));
ky=1/ky;
yn=ky^(alpha/(1-alpha));
ybar=yn*nbar;
kbar=ky*ybar;
gbar=gybar*ybar;
par.gbar = gbar;
bbar=bybar*ybar;
par.bbar = bbar;

xbar=(psi0-1+delta)*kbar;
temp = (1+tauc)/(1-taul)*(1+1/varphi)/(1-alpha);
% calibration of kappa
cybar=1-xbar/ybar-gybar;
kappa=nbar^(1+1/varphi)*(sigma*temp*cybar-sigma+1);
kappa=1/kappa
par.kappa = kappa;

disp('parameterization');
par
pause;


cy=1/(sigma*kappa*nbar^(1+1/varphi))+1-1/sigma;
cy=cy/temp;
cbar=cy*ybar;
wbar=(1-alpha)*ybar/nbar;

disp('steady state calibration: ');
beta0
ybar
kbar 
xbar
cbar
gbar
bbar
1-xbar/ybar-cbar/ybar
pause;
disp('resource constraint of economy: ');
ybar
cbar+gbar+xbar
pause;

disp('steady state eqs = 0? ');
x=[cbar, ybar, kbar, xbar, nbar];
xinitial = x;
steadystate(x,par)
x0=x;
pause;

% Laffer curve tau^l
ntau=89;
taul0=0;
taul1=linspace(0.1,0.9,ntau);
results=zeros(ntau,3);

for i=1:1:ntau
	taul=taul1(i);
    par.taul = taul;
	results(i,1)=taul*100;
%	rk=(rbar-1)/(1-tauk);
    rk=(rbar-1+delta)/(1-tauk);
    x1 = fsolve(@(x)steadystate(x,par),x0);
	w=(1-alpha)*x1(2)/x1(5);
	results(i,2)=taul*w*x1(5);
	results(i,3)=taul*w*x1(5)+tauk*rk*x1(3)+tauc*x1(1);
	x0=x1;
end


figure
plot(results(:,1),results(:,2:3))
title('Tax Revenue'); 
xlabel('Labor income tax \tau^l');
legend('Labor income tax','All taxes');
pause;

% Laffer curve for capital income tax
taul = par.taulbar;
par.taul = taul;
x=xinitial;
steadystate(x,par)
x0=x;

ntau=87;
tauk1=linspace(0.1,0.9,ntau);
results=zeros(ntau,3);

for i=1:1:ntau
	tauk = tauk1(i);
    par.tauk = tauk;
	results(i,1)= tauk*100;
	
    x1 = fsolve(@(x)steadystate(x,par),x0);
	
%	rk=(rbar-1)/(1-tauk);	
    rk=(rbar-1+delta)/(1-tauk);
	results(i,2)=tauk*rk*x1(3);
	w=(1-alpha)*x1(2)/x1(5);
	results(i,3)=taul*w*x1(5)+tauk*rk*x1(3)+tauc*x1(1);
	x0=x1;
end

figure
plot(results(:,1),results(:,2:3))
title('Tax Revenue'); 
xlabel('Capital income tax \tau^k');
legend('Capital income tax','All taxes');
pause;



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% steady state
function [y] = steadystate(x,par)
	

	css=x(1);
	yss=x(2);
	kss=x(3);
	xss=x(4);
	nss=x(5);
	y=zeros(5,1);
    gbar = par.gbar;
    tauc = par.tauc;
    taul = par.taul;
    tauk = par.tauk;
    varphi = par.varphi;
    alpha = par.alpha;
    sigma = par.sigma;
    kappa = par.kappa;
    rbar = par.rbar;
    delta = par.delta;
    psi0 = par.psi;
    
	y(1)=yss-gbar-css-xss;
	temp=(1+tauc)/(1-taul)*(1+1/varphi)/(1-alpha);
	cy=1/(sigma*kappa*nss^(1+1/varphi))+1-1/sigma;
	cy=cy/temp;
	y(2)=cy-css/yss;
%	ky=(rbar-1)/(alpha*(1-tauk))+delta/alpha;
	ky=(rbar-1)/(alpha*(1-tauk))+delta/(alpha*(1-tauk));
	ky=1/ky;
	yn=ky^(alpha/(1-alpha));
	y(3)=yss-yn*nss;
	y(4)=ky*yss-kss;
	y(5)=xss-(psi0-1+delta)*kss; 
	
end
    

