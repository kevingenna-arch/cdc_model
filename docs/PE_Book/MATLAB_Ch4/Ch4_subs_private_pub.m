% ----------------------------  Ch4_subs_private_pub.m --------------------------------
%   
%   computes steady state (4.26) in Chapter 4.2 of Heer (2018)
%
%   author: Burkhard Heer
%
%   last change: June 8, 2018
%
% -------------------------------------------------------------------------



clear; clc;
%parameters

par.alpha=0.36;
par.beta1=0.96;
par.GY=0.2;     % share of public consumption in output
par.delta=0.10;	% depreciation rate
par.rho=0.6;	% initial value substitution elasticity labor/total consumption
par.rhoc=0.5;	% initial value substitution elasticity private/public consumption
par.sigma=2;
par.phi=3/4;	% weight of private consumptiono in consumption aggregator
par.kappa=0.01; % weight of leisure in CES-utility aggregator
par.lbar=0.3;   % steady-state labor supply

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
xt=aggct^(1-1/par.rho)+par.kappa*(1-lbar)^(1-1/par.rho);    % eq. (4.22)
kappa=par.phi*(1-lbar)^(1/par.rho) *aggct^(-1/par.rho)* zt^(1/ (1-1/par.rhoc) - 1) * cpbar^(-1/par.rhoc); 
kappa=kappa/wage_eq_4_9(kbar,lbar,par.alpha);   % eq. (4.26)
kappa0=kappa;
kappa0
pause;
% steady state solution
x0=[kbar, lbar, cpbar];
xinitial=x0;
 
xkappa0=[kbar, kappa0, cpbar];
disp('steadystate1(x0): ');

steady_state1_eq_4_26(par,xkappa0)
pause;

% solution of non-linear equation system (4.26) with endogenous kappa
% k, and C^p for exogenous l=0.3
kappa_solution = fsolve(@(x)steady_state1_eq_4_26(par,x),xkappa0);
kappa_solution

par.kappa=kappa_solution(2);    

% initial solution guess for steady state (4.26) with endogenous 
% k, l, and C^p for exogenonus kappa
x0=kappa_solution;
x0(2)=par.lbar;
xinitial0=x0;

steady_state0_eq_4_26(par,x0)

%
% computation of Fig. 4.9, kappa(rho^c)
%

ngrid=50;
rhocgrid1=linspace(par.rhoc,0.99,ngrid);
rhocgrid2=linspace(par.rhoc,0.01,ngrid);
rhocgrid=[rhocgrid1 rhocgrid2];
kappagrid=rhocgrid;

xinitial=xinitial0;
for i=1:2*ngrid
    par.rhoc=rhocgrid(i);
    if (i==ngrid) 
        xinitial=xinitial0;
    end
    
    x0 = fsolve(@(x)steady_state1_eq_4_26(par,x),xinitial);
    xinitital=x0;
    kappagrid(i)=x0(2);
    i
    kappagrid(i)
end

solution=[rhocgrid; kappagrid];
solution=solution';
solution=sort(solution,1);

figure
plot(solution(:,1),solution(:,2)); hold on
xlabel('\rho_c')
ylabel('\kappa')

%
% computation of the effects of a 1% increase of government consumption
% for all rho_c in (0,1)
%

gtilde=gbar*1.01;
par.gtilde=gtilde;
result=zeros(2*ngrid,4);

xinitial=xinitial0;
for i=1:2*ngrid
    par.rhoc=rhocgrid(i);
    if (i==ngrid) 
        xinitial=xinitial0;
    end
    
    par.kappa=solution(i,2);
    par.rhoc=solution(i,1);
    x0 = fsolve(@(x)steady_state2_eq_4_26(par,x),xinitial);
    xinitital=x0;
    i
    x0
    
	result(i,1)=par.rhoc;
	result(i,2)=(x0(2)-lbar)/lbar*100;
	result(i,3)=(production_eq_4_8(x0(1),x0(2),par.alpha)-ybar)/ybar*100;
	result(i,4)=(x0(3)-cpbar)/cpbar*100;
end 

result1=sort(result,1);

figure(1)
subplot(2,1,1),plot(result(:,1),result(:,2),'k','LineWidth',3), title('Labor L'); hold on
xlabel('\rho_c')
subplot(2,1,2),plot(result(:,1),result(:,3),'k','LineWidth',3), title('Production Y'); hold on
xlabel('\rho_c')
