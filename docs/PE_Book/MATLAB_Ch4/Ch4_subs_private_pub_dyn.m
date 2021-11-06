% ----------------------------  Ch4_subs_private_pub_dyn.m --------------------------------
%   
%   computes steady state (4.26) and tranitionin Chapter 4. of Heer (2018)
%
%   author: Burkhard Heer
%
%   last change: June 14, 2018
%
% -------------------------------------------------------------------------



clear; clc;
%parameters

par.alpha=0.36;
par.beta1=0.96;
par.GY=0.2;     % share of public consumption in output
par.delta=0.10;	% depreciation rate
par.rho=0.6;	% initial value substitution elasticity labor/total consumption
par.rhoc=0.3;	% initial value substitution elasticity private/public consumption
par.sigma=2;
par.phi=3/4;	% weight of private consumptiono in consumption aggregator
par.kappa=0.01; % weight of leisure in CES-utility aggregator
par.lbar=0.3;   % steady-state labor supply
par.nt=40;      % number of transition periods

% eq. (4.26a)
KL=(par.alpha/(1/par.beta1-1+par.delta))^(1/(1-par.alpha));


% initial guess
lbar=par.lbar;
kbar=KL*lbar;
ybar=production_eq_4_8(kbar,lbar,par.alpha);
gbar=par.GY*ybar;
cpbar=ybar-gbar-par.delta*kbar; % eq. (4.11)
zt=par.phi* cpbar^(1-1/par.rhoc) + (1-par.phi) * gbar^(1-1/par.rhoc); % eq. (4.25)
aggct=zt^(1/(1-1/par.rhoc)); % eq. (4.21)
aggctbar=aggct;
xt=aggct^(1-1/par.rho)+par.kappa*(1-lbar)^(1-1/par.rho);    % eq. (4.22)
kappa=par.phi*(1-lbar)^(1/par.rho) *aggct^(-1/par.rho)* zt^(1/ (1-1/par.rhoc) - 1) * cpbar^(-1/par.rhoc); 
kappa=kappa/wage_eq_4_9(kbar,lbar,par.alpha);   % eq. (4.26)
kappa0=kappa;
kappa0

% steady state solution
x0=[kbar, lbar, cpbar];
xinitial=x0;
 
xkappa0=[kbar, kappa0, cpbar];
disp('steadystate1(x0): ');

steady_state1_eq_4_26(par,xkappa0)


% solution of non-linear equation system (4.26) with endogenous kappa
% k, and C^p for exogenous l=0.3
kappa_solution = fsolve(@(x)steady_state1_eq_4_26(par,x),xkappa0);
kappa_solution

par.kappa=kappa_solution(2)    


% initial solution guess for steady state (4.26) with endogenous 
% k, l, and C^p for exogenonus kappa
x0=kappa_solution;
x0(2)=par.lbar;
xinitial0=x0;

steady_state0_eq_4_26(par,x0)


%
%   final steady state for higher government consumption
%
gtilde=gbar*1.01;
par.gtilde=gtilde;
x0 = fsolve(@(x)steady_state2_eq_4_26(par,x),xinitial0);


xtilde=x0;
cptilde=x0(3)
ltilde=x0(2)
ktilde=x0(1)
ytilde=production_eq_4_8(ktilde,ltilde,par.alpha);
zttilde=par.phi* cptilde^(1-1/par.rhoc) + (1-par.phi) * gtilde^(1-1/par.rhoc);
aggcttilde=zttilde^(1/(1-1/par.rhoc));
xttilde=aggcttilde^(1-1/par.rho)+par.kappa*(1-ltilde)^(1-1/par.rho);

%
%   computation of the transition
%

par.ltilde=ltilde;
par.ktilde=ktilde;
par.cptilde=cptilde;
par.gtilde=gtilde;
par.aggcttilde=aggcttilde;
par.case=1;             % case 1 -- permanent higher government consumption

govt=gtilde*ones(par.nt+1,1);   % government consumption during transition
par.lfinal=ltilde;
par.kfinal=ktilde;
par.cpfinal=cptilde;
par.gfinal=gtilde;
par.aggctfinal=aggcttilde;


% check findlcpk
par1=par;
par1.govt=govt;
par1.lold=ltilde;
par1.kold=ktilde;
par1.cpold=cptilde;
par1.k0 = ktilde;       % needs to be saved in findkinitial
par1.kbar = kbar;
par1.lbar = lbar;
par1.cpbar = cpbar;
par1.aggctbar = aggctbar;
par1.gbar = gbar;
par1.ybar = ybar;
par1.i = 10;

% steady state solution
x0=[ktilde, ltilde, cptilde];
x1=[ltilde, cptilde];
xinitial=x0;
findlcpk(par1,x0)
findlcp(par1,x1)

% find initial value for findkinitial

ngrid = 100;
kgrid = linspace(0.999*ktilde,ktilde,ngrid);

ymin = 100;
for i=1:1:ngrid
    clc;
    y0= findkinitial(par1,kgrid(i));
    if abs(y0)<ymin
        kinit=kgrid(i);
        ymin = abs(y0);
    end
end


%kinit = 1.2946273;
y0=findkinitial(par1,kinit)

x1 = fsolve(@(k0)findkinitial(par1,k0),kinit);
x1	

transpath=transitionpathch44(par1,x1);

figure
subplot(3,2,1);
plot(transpath(:,7),(transpath(:,1)/par1.kbar-1)*100)
title('Capital k'); 
xlabel('Period');
subplot(3,2,2);
plot(transpath(:,7),(transpath(:,2)/par1.lbar-1)*100) 
title('Labor l'); 
xlabel('Period');
subplot(3,2,3);
plot(transpath(:,7),(transpath(:,3)/par1.ybar-1)*100)
title('Production y'); 
xlabel('Period');
subplot(3,2,4);
plot(transpath(:,7),(transpath(:,4)/par1.cpbar-1)*100)
title('Consumption c^p'); 
xlabel('Period');
subplot(3,2,5);
plot(transpath(:,7),(transpath(:,5)/par1.aggctbar-1)*100)
title('Total consumption c'); 
xlabel('Period');
subplot(3,2,6);
plot(transpath(:,7),(transpath(:,6)/par1.gbar-1)*100) 
title('Government consumption g'); 
xlabel('Period');


