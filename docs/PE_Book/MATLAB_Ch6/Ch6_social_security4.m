function Ch6_social_security4();
% 
%	Chapter 6.3.5 PAYG System and Growth
%   Author: Burkhard Heer
%   Last Change: July 5, 2018


clc;
close all;


% Parameterization 
%

growth=0.80;    % growth rate of labor-augmenting technical change
ls=0.3;			% labor supply
alpha=0.36;		% production elasticity of capital
beta0=0.40;		% discount factor 
n=0.1;			% population growth rate
nu0=257.15;
nu1=3.33;
% case 1: no pension
tau=0;
g=0;

par.growth = 0.80;
par.ls = ls;
par.alpha = alpha;
par.beta0 = beta0;
par.n = n;
par.nu0 = nu0;
par.nu1 = nu1;
par.tau = tau;
par.g = g;

% initial value from social-security2.m
lss=0.3;
kss=0.0182;		% solution from social_security2.g
rss=interest(kss,lss,par);
wss=wage(kss,lss,par);
c1ss=(1-tau)*wss/(nu0*lss^(nu1));
c2ss=beta0*(1+rss)*c1ss;
dss=0;


x0=[kss, lss, wss, rss, c1ss, c2ss, dss];

xss = fsolve(@(x)ksteady(x,par),x0);
disp('ksteady(x)=0? ');
ksteady(xss,par)
pause;
disp('steady state values k,l: ');
x0=xss;
xss
pause;


% tatonnement over g
ng=100;
gseq=linspace(0,growth,ng);

for i=1:1:ng
	g=gseq(i);
    par.g = g;
    x = fsolve(@(x)ksteady(x,par),x0);
	x0=x;
end


disp('steady state solution for g=growth');
x0
pause;

% calibration of nu0
g=growth;
x0(2)=nu0;
    x = fsolve(@(x)calibnu0(x,par),x0);
x0=x;
nu0=x(2);
par.nu0 = nu0;
nu0
pause;

x0(2)=0.3;
disp('test: ksteady(x0)=0? ');;
ksteady(x0,par)
pause;
x0;
pause;

kss=x0(1);
lss=x0(2);
c1ss=x0(5);
c2ss=x0(6);
utilss=log(c1ss)+beta0*log(c2ss)-nu0*lss^(1+nu1)/(1+nu1);
yss=kss^(alpha)*lss^(1-alpha);
yss
pause;


% case tau=0.3
tau=0.3;
par.tau = tau;
x1 = fsolve(@(x)ksteady(x,par),x0);


kss1=x1(1);
lss1=x1(2);
c1ss1=x1(5);
c2ss1=x1(6);
utilss1=log(c1ss1)+beta0*log(c2ss1)-nu0*lss1^(1+nu1)/(1+nu1);
yss1=kss1^(alpha)*lss1^(1-alpha);
yss1
pause;

cec=exp( (utilss1-utilss)/(1+beta0) )-1;
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



% computes steady state for pension case and growth
function [y] = ksteady(x,par)

    k=x(1);
	l=x(2);
	w=x(3);
	r=x(4);
	c1=x(5);
	c2=x(6);
	d=x(7);
    n = par.n;
    tau = par.tau;
    nu0 = par.nu0;
    nu1 = par.nu1;
    beta0 = par.beta0;
    g = par.g;
	
	y=zeros(7,1);
	
	y(1)=c1*(1+r)*beta0-(1+g)*c2;
	y(2)=nu0*l^(nu1)*c1-(1-tau)*w;
	y(3)=w-wage(k,l,par);
	y(4)=r-interest(k,l,par);
	y(5)=k*(1+g)*(1+n)-beta0/(1+beta0)*(1-tau)*w*l+1/(1+beta0)*(1+g)/(1+r)*d;
	y(6)=d-(1+n)*tau*w*l;
	y(7)=c1+c2*(1+g)/(1+r)-(1-tau)*w*l-(1+g)/(1+r)*d;
	
end


% calibrates nu_0
function [y] = calibnu0(x,par)

    k=x(1);
	nu0=x(2);
	w=x(3);
	r=x(4);
	c1=x(5);
	c2=x(6);
	d=x(7);
    n = par.n;
    tau = par.tau;
    l = par.ls;
    nu1 = par.nu1;
    beta0 = par.beta0;
    g = par.g;
	
	y=zeros(7,1);
	
	y(1)=c1*(1+r)*beta0-(1+g)*c2;
	y(2)=nu0*l^(nu1)*c1-(1-tau)*w;
	y(3)=w-wage(k,l,par);
	y(4)=r-interest(k,l,par);
	y(5)=k*(1+g)*(1+n)-beta0/(1+beta0)*(1-tau)*w*l+1/(1+beta0)*(1+g)/(1+r)*d;
	y(6)=d-(1+n)*tau*w*l;
	y(7)=c1+c2*(1+g)/(1+r)-(1-tau)*w*l-(1+g)/(1+r)*d;
	
end
