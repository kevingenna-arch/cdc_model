function problem_P6_1();
% 
%	Chapter 6.3.2 PAYG System
%   Author: Burkhard Heer
%   Last Change: June 5, 2019
%
%   Sensitivity Analysis: Computes the case for an intertemporal elasticity
%                           of substitution IES=0.5
%   Exercise, see youtube video tutorial
%


clc;
close all;

%
% Step 1: Parameterization 
%

ls=0.3;			% labor supply
alpha=0.36;		% production elasticity of capital
beta0=0.40;		% discount factor 
n=0.1;			% population growth rate
nt=20;          % number of transition periods
sigma0 = 2.0;   % coefficient of relative risk aversion
% disutility parameters: see Example 2
nu0=257.15;
nu1=0.3;
delta=0.5;
%
% Step 2:   Computaion of the steady state for the
%           case 1: no pension
%
tau=0;

par.ls = ls;
par.alpha = alpha;
par.beta0 = beta0;
par.n = n;
par.nu0 = nu0
par.nu1 = nu1;
par.tau = tau;
par.nt = nt;
par.sigma0 = sigma0;
par.delta=delta;


% compute a grid for the solution for various k
% eq. (6.19)
kmin=0.001;
kmax=5.0;
nk=1000;
k=linspace(kmin,kmax,nk);
yk=zeros(nk,1);

for i=1:1:nk
    par.i = i;
	yk(i)=ksteady(k(i),par);
	i
    yk(i)
end

clc;
yk = abs(yk);
% minimum of yk as starting value
[M, I] = min(yk)

kinitial = k(I);
kss = fsolve(@(x)ksteady(x,par),kinitial);

disp('ksteady(kss,par) ');
ksteady(kss,par)
kss
yss=kss^(alpha)*ls^(1-alpha)

rss = interest(kss,ls,par);
wss = wage(kss,ls,par);
c1ss = ( (1-tau)+ (1+n)/(1+rss)*tau) *wss*ls / ( 1 +beta0^(1/sigma0) * (1+rss)^(1/sigma0-1)  ) ;
c2ss=beta0^(1/sigma0)*c1ss*(1+rss)^(1/sigma0);
utilss=util(c1ss,c2ss,ls,par)


%
% Step 3:   Computaion of the steady state for the
%           case 2: PAYG pension 
% 
tau=0.3;
par.tau = tau;

for i=1:1:nk
    par.i = i;
	yk(i)=ksteady(k(i),par);
	i
    yk(i)
end

clc;
yk = abs(yk);
% minimum of yk as starting value
[M, I] = min(yk)

kinitial = k(I);
kssd = fsolve(@(x)ksteady(x,par),kinitial);

pause;

disp('ksteady(kssd,par) ');
ksteady(kssd,par)
kssd
yssd=kssd^(alpha)*ls^(1-alpha)

rssd = interest(kssd,ls,par);
wssd = wage(kssd,ls,par);
c1ssd = ( (1-tau)+ (1+n)/(1+rssd)*tau) *wssd*ls / ( 1 +beta0^(1/sigma0) * (1+rssd)^(1/sigma0-1)  ) ;
c2ssd = beta0^(1/sigma0)*c1ssd*(1+rssd)^(1/sigma0);
utilssd=util(c1ssd,c2ssd,ls,par);


if sigma0==1
    cec=exp( (utilssd-utilss)/(1+beta0) )-1
else
    cec=(utilssd-utilss)/( c1ss^(1-sigma0)/(1-sigma0) + beta0 * c2ss^(1-sigma0)/(1-sigma0) );
    cec=cec+1;    
    cec = cec^(1/(1-sigma0))-1
end
pause;

%
% Step 4:
% computation of the dynamics
%
kt=zeros(nt+1,1);
utilt=kt;
periods=linspace(0,nt,21);
kt(1)=kssd;
utilt(1)=utilssd;
cect=zeros(nt+1,1);

cec=(utilssd-utilss)/( c1ss^(1-sigma0)/(1-sigma0) + beta0 * c2ss^(1-sigma0)/(1-sigma0) );
cec=cec+1;    
cec = cec^(1/(1-sigma0))-1;
cect(1)=cec;


tau0=0.3;
par.tau0 = tau0;    % social security rate in period t
tau1=0.3;           % social security rate in period t+1
par.tau1 = tau1;
for i=1:1:nt
    if i==2 
        tau1=0; % first generation still has to pay contributions for pension of the old
                % but will not receive a pension itself
        par.tau1 = tau1;
    elseif i==3
        tau0=0;     % generation born in period t=2 (corresponding to i=3) is the first generation
                    % that does not pay social security contributions
        par.tau0=tau0;
    end	 
	k0=kt(i);
    par.k0 = k0;
    k1 = fsolve(@(x)kdyn(x,par),k0); 
	d0 = tau0*wage(k0,ls,par)*ls;   % social security contribution in period t=i-1
    d1 = tau1*wage(k1,ls,par)*ls;   % social security contribution in period t+1 = i
    
    c1 = (1-tau0)*wage(k0,ls,par) *ls+ (1+n)*d1/(1+interest(k1,ls,par));
    c1 = c1 / (1 + beta0^(1/sigma0)*(1+interest(k1,ls,par))^(1/sigma0-1));
        
	c2=beta0^(1/sigma0)*c1*(1+interest(k1,ls,par))^(1/sigma0);
    
    
	utilt(i+1) = util(c1,c2,ls,par);
    
    cec=(util(c1,c2,ls,par)-utilss)/( c1ss^(1-sigma0)/(1-sigma0) + beta0 * c2ss^(1-sigma0)/(1-sigma0) );
    cec=cec+1;    
    cec = cec^(1/(1-sigma0))-1;
    
	cect(i+1) = cec;
    kt(i+1)=k1;
end

clc;


save Ch6_social_security_IES05;

figure
plot(periods,kt);
xlabel('Period t');
ylabel('Capital K_t');
pause;

figure
plot(periods,cect*100);
ylabel('Welfare \Delta_t');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% interest rate 
function [y] = interest(k,l,par)
    alpha = par.alpha;
	y=alpha*k^(alpha-1)*l^(1-alpha)-par.delta;
end

%% wage 
function [w] = wage(k,l,par)
    alpha = par.alpha;
	w=(1-alpha)*k^alpha*l^(-alpha);
end


% procedure that computes k_t+1 given k_t*/
% input: k1 - next period capital stock 
% output: the value of the equilibrium condition
function [y] =  kdyn(x,par)
    beta0 = par.beta0;
	n = par.n;
    ls = par.ls;
    tau0 = par.tau0;
    tau1 = par.tau1;     
    sigma0 = par.sigma0;
    k1 = x;
    k0 = par.k0;
    w0=wage(k0,ls,par);
    r1=interest(k1,ls,par);
	d1=tau1*wage(k1,ls,par)*ls;
	y=k1*(1+n)-(1-tau0)*w0*ls*beta0^(1/sigma0)*(1+r1)^(1/sigma0-1)/(1+beta0^(1/sigma0)*(1+r1)^(1/sigma0-1));
    y = y + (1+n)*d1 / (1+r1+beta0^(1/sigma0)*(1+r1)^(1/sigma0) );
    
end

%% ksteady: eq. (6.19) 
function [y] = ksteady(k0,par)
    
    tau = par.tau;
    ls = par.ls;
    i = par.i;
    n = par.n;
    beta0 = par.beta0;
    sigma0 = par.sigma0;
    
	w=wage(k0,ls,par);
	r=interest(k0,ls,par);
	
	y=k0*(1+n)-(1-tau)*w*ls*beta0^(1/sigma0)*(1+r)^(1/sigma0-1)/ ( 1+ beta0^(1/sigma0)*(1+r)^(1/sigma0-1));
	y = y + (1+n)*tau*w*ls/ (1+ r+beta0^(1/sigma0)*(1+r)^(1/sigma0) );
    
end

%% util: lifetime utility 
function[y] = util(c1,c2,l1,par)
	sigma0 = par.sigma0;
    beta0 = par.beta0;
    nu0 = par.nu0;
    nu1 = par.nu1;
    if nu1==1
        y = ln(c1)+beta0*ln(c2)- nu0 * l1^(1+1/nu1)/(1+1/nu1);
    else
	y = c1^(1-sigma0)/(1-sigma0) + beta0 * c2^(1-sigma0)/(1-sigma0) - nu0 * l1^(1+1/nu1)/(1+1/nu1);
    end
end