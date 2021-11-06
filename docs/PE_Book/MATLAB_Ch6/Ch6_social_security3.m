function Ch6_social_security3();
% 
%	Chapter 6.3.4 PAYG System with Contribution-Based Benefits
%   Author: Burkhard Heer
%   Last Change: July 5, 2018


clc;
close all;


% Parameterization 
%

ls=0.3;			% labor supply
alpha=0.36;		% production elasticity of capital
beta0=0.40;		% discount factor 
n=0.1;			% population growth rate
theta=0.5;      % parameter on contributions in pension formula
nu0=257.15;
nu1=3.33;
% case 1: no pension
tau=0.3;

par.ls = ls;
par.alpha = alpha;
par.beta0 = beta0;
par.n = n;
par.nu0 = nu0
par.nu1 = nu1;
par.tau = tau;
par.theta = theta;



% initial values from Numerical Example 2, social_security2.m
kss=0.0182;
lss=0.3;
wss=wage(kss,lss,par);
rss=interest(kss,lss,par);
c1ss=0.0499;
penminss=tau*wss*lss*(1+n-theta);
c2ss=c1ss*beta0*(1+rss);

x0=[kss, lss, wss, rss, c1ss, c2ss, penminss];
xss = fsolve(@(x)steadystate(x,par),x0);
disp('steadystate(xss)=0? ');
steadystate(xss,par)
pause;
xss
pause;


	k=xss(1);
	l=xss(2);
	w=xss(3);
	r=xss(4);
	c1=xss(5);
	c2=xss(6);
	penmin=xss(7);
	
	y0=k^(alpha)*l^(1-alpha);
	
	
	tau=0;
    par.tau = 0;
    xss0 = fsolve(@(x)steadystate(x,par),xss);
disp('steadystate(xss0)=0? ');
steadystate(xss0,par)
pause;
xss0
pause;


	k0=xss0(1);
	l0=xss0(2);
	w0=xss0(3);
	r0=xss0(4);
	c10=xss0(5);
	c20=xss0(6);
	penmin0=xss0(7);
	
	
	y00=k0^(alpha)*l0^(1-alpha);
	
	utility0=log(c10)+beta0*log(c20)-nu0*l0^(1+nu1)/(1+nu1);
	utility1=log(c1)+beta0*log(c2)-nu0*l^(1+nu1)/(1+nu1);
	
	cec=(utility1-utility0)/(1+beta0);
	disp('consumption equivalent change: ');
	cec=exp(cec)-1;
	cec



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% interest rate 
function [y] = interest(k,l,par)
    alpha = par.alpha;
	y=alpha*k^(alpha-1)*l^(1-alpha);
end

%% wage 
function [w] = wage(k,l,par)
    alpha = par.alpha;
	w=(1-alpha)*k^alpha*l^(-alpha);
end


% steady state (6.36)
function [y] = steadystate(x,par)

	k=x(1);
	l=x(2);
	w=x(3);
	r=x(4);
	c1=x(5);
	c2=x(6);
	penmin=x(7);
	
    theta = par.theta;
    n = par.n;
    tau = par.tau;
    nu0 = par.nu0;
    nu1 = par.nu1;
    beta0 = par.beta0;
    
	y=zeros(7,1);
	y(1)=penmin-tau*w*l*(1+n-theta);
	y(2)=w-wage(k,l,par);
	y(3)=r-interest(k,l,par);
	y(4)=nu0*l^(nu1)*c1-(1-tau)*w-theta*tau*w/(1+r);
	y(5)=c1*beta0*(1+r)-c2;
	y(6)=(1-tau)*w*l+(penmin+theta*tau*w*l)/(1+r)-c1-c2/(1+r);
	y(7)=(1+n)*k-(1-tau)*w*l+c1;

end
