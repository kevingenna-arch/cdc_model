function Ch5_welfare_taul();
% 
%	Chapter 5.3 Capital Income Tax
%   Author: Burkhard Heer
%   Last Change: June 22, 2018


clc;
close all;


% Parameterization 
%

par.depreciation=1;		% if ==1, depreciation is tax deductible
						% if ==0, depreciation is not tax-deductible
par.alpha=0.36;
par.beta1=0.96;		% annual real interest rate of 4%
par.lstar=0.30;		% steady-state labor supply
par.sigma=2.0;		% 1/IES
par.taul=0.23;
par.tauk = 0.41;
%par.taul=0.40;     % tax wedge US from Prescott (2004)
%par.taul=0.59;     % tax wedge Germany
par.taul0=par.taul; % taul in initial steady state
par.delta=0.1;
tauk=par.tauk;
delta=par.delta;
beta1=par.beta1;
alpha=par.alpha;
lstar=par.lstar;
taul=par.taul;

if par.depreciation==1
	rstar=( 1/beta1-1 )/(1-tauk) +delta;
else
	rstar=( 1/beta1-1 +delta)/(1-tauk);
end

klstar=(alpha/rstar)^(1/(1-alpha));
kstar=klstar*lstar;
wstar=wage(kstar,lstar,par);
ystar=production(kstar,lstar,par);
par.kstar=kstar;
if par.depreciation==1
	gstar=taul*wstar*lstar+tauk*(rstar-delta)*kstar;
else
	gstar=taul*wstar*lstar+tauk*rstar*kstar;
end
par.gstar=gstar;

disp('government share:  gstar/ystar');
gstar/ystar

cstar=ystar-gstar-delta*kstar;

temp=(1-taul)*wstar*(1-lstar)/cstar;
gam=1/(1+temp)
par.gam=gam;


disp('parameterization');
par

% test of function steadystate
xinitial=[lstar, cstar, kstar , wstar ,rstar, taul];
disp('steady state conditions==0? ');
y=steadystate(xinitial,par)
pause;


%
%   steady state effects of tauk
%

ngrid=80;
taukgrid = linspace(-0.1,0.8,ngrid);
solution=zeros(ngrid,9);

for i=1:1:ngrid
	tauk=taukgrid(i);
	par.tauk=tauk;
    sol = fsolve(@(x)steadystate(x,par),xinitial);
	solution(i,1)=tauk;
	solution(i,2:7)=sol';
	% consumption equivalent change
	ce1=(util(cstar,lstar,par)/ util(sol(2),sol(1),par) )^(-1/(par.gam*(1-par.sigma)))-1;
	solution(i,8)=production(sol(1),sol(3),par);
	solution(i,9)=ce1;
	
end

% plotting of the solution
        figure
        subplot(3,2,1);
        plot(solution(:,1)*100,solution(:,4))
        title('Capital k'); 
        xlabel('Period');
        subplot(3,2,2);
        plot(solution(:,1)*100,solution(:,2)) 
        title('Labor l'); 
        xlabel('Period');
        subplot(3,2,3);
        plot(solution(:,1)*100,solution(:,8))
        title('Production y'); 
        xlabel('Period');
        subplot(3,2,4);
        plot(solution(:,1)*100,solution(:,3))
        title('Consumption C'); 
        xlabel('Period');
        subplot(3,2,5);
        plot(solution(:,1)*100,solution(:,9))
        title('Welfare'); 
        xlabel('Period');
        subplot(3,2,6);
        plot(solution(:,1)*100,solution(:,7))
        title('Labor income tax rate \tau^L'); 
        xlabel('Period');


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% utility function
function [y] = util(c,labor,par)
	y=c^(par.gam)*(1-labor)^(1-par.gam);
	y=y^(1-par.sigma)/(1-par.sigma);
end

%% steady state
function [y] = steadystate(x,par)
	
	labor=x(1);
	c=x(2);
	k=x(3);
	w=x(4);
	r=x(5);
	taul=x(6);
	y=zeros(6,1);
    
    gam=par.gam;
    delta=par.delta;
    beta1=par.beta1;
    tauk=par.tauk;
    gstar=par.gstar;
    
	y(1)=(1-gam)/gam*c-(1-labor)*(1-taul)*w;
	y(2)=c+delta*k+gstar-production(k,labor,par);
	y(3)=r-interest(k,labor,par);
	y(4)=w-wage(k,labor,par);
	if par.depreciation==1
		y(5)=1/beta1-1+(1-tauk)*(delta-r);
		y(6)=gstar-k*(r-delta)*tauk-w*labor*taul;
    else
		y(5)=1/beta1-1+delta-(1-tauk)*r;
		y(6)=gstar-k*r*tauk-w*labor*taul;
    end
	
end
    


% production function F
function [y] = production(k,l,par)
    alpha = par.alpha;
	y=  (k^alpha)*l^(1-alpha);
end

% wage
function [y] = wage(k,l,par)
    alpha = par.alpha;
	y= (1-alpha)*(k^alpha)*l^(-alpha);
end

% interest rate
function [y] = interest(k,l,par)
    alpha = par.alpha;
	y= alpha*k^(alpha-1)*l^(1-alpha);
end


