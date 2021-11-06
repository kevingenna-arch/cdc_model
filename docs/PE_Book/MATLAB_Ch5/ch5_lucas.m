function Ch5_lucas();
% 
%	Chapter 5.5 Laffer Curve
%   Author: Burkhard Heer
%   Last Change: July 3, 2018


clc;
close all;


% Parameterization 
%

grate=0.015;	% growth rate
beta0=0.96;		% annual discount factor
delta=0.08;		% depreciation rate
gam=0.8;		% elasticity of growth rate w.r.t. learning
alpha=0.361;	% share of capital
sigma=2.0;		% IES
sigmap=-2/3;	% production elasticity
ubar=0.3;		% steady-state labor supply
gybar=0.19;		% gov consumption share
taukbar=0.41;	% tax rate capital
taulbar=0.28;   % tax rate labor 
a0=1.0;         % production parameter


% implied values
rbar=(1+grate)^(sigma)/beta0-1+delta;
rbar=rbar/(1-taukbar);

% storing calibrated values to the variable par
par.grate = grate;
par.beta0 = beta0;
par.delta = delta;
par.gam = gam;
par.alpha = alpha;
par.sigma = sigma;
par.sigmap = sigmap;
par.ubar = ubar;
par.gybar = gybar;
par.taukbar = taukbar;
par.taulbar = taulbar;
par.a0 = a0;
par.rbar = rbar;

k0=0.1;
    k1 = fsolve(@(k)findk(k,par),k0);
kbar=k1;
par.kbar = kbar;
par
pause;


wbar=wage(kbar,ubar,par);
ybar=production(kbar,ubar,par);
gbar=gybar*ybar;
cbar=ybar-gbar-(grate+delta)*kbar;

% calibration of D
% initial values for v and D
vbar=grate*ubar*gam/ ( (1+grate)^(sigma-1)/beta0-1);
Dbar=grate/(vbar^gam);
vD0=[vbar, Dbar];
    vD = fsolve(@(x)findvd(x,par),vD0);
vbar=vD(1)
Dbar=vD(2)

par.Dbar = Dbar;
par.vbar = vbar;
par.ybar = ybar;
par.gbar = gbar;

% calibration of rho (iota in the text)
rho=(1-taulbar)*wbar*(1-ubar-vbar)/cbar;
trbar=taulbar*wbar*ubar+taukbar*rbar*kbar-gbar;
par.rho = rho
par.trbar = trbar;

par.tauk = taukbar;
x0=[grate, vbar, ybar, kbar, ubar, wbar, rbar, cbar, taulbar];

disp('Steady state eqs = 0?');
steadystate(x0,par)
pause;
xinitial=x0;

% computation of the steady state effects of lower capital taxes
%
ntau=40;
taukgrid = linspace(0.41-0.01,0.02,ntau);
results=zeros(ntau,10);

for i=1:1:ntau
	tauk=taukgrid(i);
    par.tauk = tauk;
    x1 = fsolve(@(x)steadystate(x,par),x0);
	x0=x1;
	results(i,1) = tauk;
	results(i,2:10)=x1';
end


ntau=39;
taukgrid=linspace(taukbar,0.79,ntau);
results1=zeros(ntau,10);
x0=xinitial;

for i=1:1:ntau
	tauk=taukgrid(i);
    par.tauk = tauk;
    x1 = fsolve(@(x)steadystate(x,par),x0);
	x0=x1;
	results1(i,1) = tauk;
	results1(i,2:10)=x1';
end
clc;

y=vertcat(results,results1);
y1 = sortrows(y,1);
% transformation in percentage points
y=y1;
y(:,1) = y(:,1)*100;
y(:,2) = y(:,2)*100;

figure
plot(y(:,1),y(:,2));
xlabel('Capital income tax rate \tau^k');
title('Growth rate \gamma');
pause;

        figure
        subplot(3,2,1);
        plot(y(:,1),y(:,10))
        title('Labor income tax rate \tau^l'); 
        xlabel('Capital income tax rate \tau^k');
        subplot(3,2,2);
        plot(y(:,1),y(:,6)) 
        title('Working hours u'); 
        xlabel('Capital income tax rate \tau^k');
        subplot(3,2,3);
        plot(y(:,1),y(:,3))
        title('Learning v'); 
        xlabel('Capital income tax rate \tau^k');
        subplot(3,2,4);
        plot(y(:,1),y(:,5))
        title('Capital K/H'); 
        xlabel('Capital income tax rate \tau^k');
        subplot(3,2,5);
        plot(y(:,1),y(:,4))
        title('Output Y/H');
        xlabel('Capital income tax rate \tau^k');
        subplot(3,2,6);
        plot(y(:,1),y(:,9))
        title('Consumption C/H');
        xlabel('Capital income tax rate \tau^k');
       

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% implicit function to find k from interest rate
function [y] = findk(k,par)
	rbar = par.rbar;
    ubar = par.ubar;
	y=rbar-interest(k,ubar,par);
end

%% CES production function 
function [y] = production(k,u,par)
    a0 = par.a0;
    alpha = par.alpha;
    sigmap = par.sigmap;
	y=a0*(alpha*k^(sigmap)+(1-alpha)*u^(sigmap))^(1/sigmap);
end

%% interest rate for CES production function
function [y] = interest(k,u,par)
    a0 = par.a0;
    alpha = par.alpha;
    sigmap = par.sigmap;
	y=a0*alpha*(alpha*k^(sigmap)+(1-alpha)*u^(sigmap))^(1/sigmap-1)*k^(sigmap-1);
end

%% wage for CES production function
function [w] = wage(k,u,par)
	a0 = par.a0;
    alpha = par.alpha;
    sigmap = par.sigmap;
	w=a0*(1-alpha)*(alpha*k^(sigmap)+(1-alpha)*u^(sigmap))^(1/sigmap-1)*u^(sigmap-1);
end

%% calibration of parameter D in human capital accumulation function
function [y] = findvd(x,par)
	v=x(1);
	D=x(2);
	y=zeros(2,1);
    grate = par.grate;
    gam = par.gam;
    sigma = par.sigma;
    beta0 = par.beta0;
    ubar = par.ubar;
	y(1)=grate-D*v^gam;
	y(2)=(1+grate)^(sigma-1)/beta0-1-gam*D*v^(gam-1)*ubar;
end

%% steady state
function [y] = steadystate(x,par)
	

	y=zeros(9,1);
	g=x(1);
	v=x(2);
	y0=x(3);    % production
	k=x(4);
	u=x(5);
	w=x(6);
	r=x(7);
	c=x(8);
	taul=x(9);
    
    gbar = par.gbar;
    trbar = par.trbar;
    delta = par.delta;
    tauk = par.tauk;
    Dbar = par.Dbar;
    sigma = par.sigma;
    beta0 = par.beta0;
    gam = par.gam;
    rho = par.rho;
    
	y(1)=r-interest(k,u,par);
	y(2)=w-wage(k,u,par);
	y(3)=y0-production(k,u,par);
	y(4)=y0-c-gbar-(delta+g)*k;
	y(5)=taul*w*u+tauk*r*k-gbar-trbar;
	y(6)=g-Dbar*v^(gam);
	y(7)=(1+g)^(sigma-1)/beta0-1-gam*Dbar*v^(gam-1)*u;
	y(8)=(1+g)^(sigma)/beta0-1+delta-(1-tauk)*r;
	y(9)=(1-taul)*w-rho*c/(1-u-v);
	
end
    

