function Ch7_debt1();
% 
%	Chapter 7.4.2 OLG and Public Debt
%   Author: Burkhard Heer
%   Last Change: July 6, 2018


clc;
close all;


% Parameterization 
%

alpha=0.36;         % production elasticity of capital
beta1=0.96^(30);	% discount factor 
eta=2.0;            % intertemporal elasticity of substitution
a0=1;               % productivity
delta=1;            % deprecation rate
ga=0.81;            % growth rate
GY=0.20;            %government consumption share
GY=0;
BY=1.00;            % annual debt-GDP ratio
BY=BY/30;           % 30-year value
BY=0;
n=0.2;              % population growth rate

par.eta = eta;
par.a0 = a0;
par.alpha = alpha;
par.beta1 = beta1;
par.delta = delta;
par.n = n;
par.ga = ga;
par.BY = BY;
par.GY = GY;


% calibration
% initial guess for k
nk=1000;
kgrid=linspace(0,1,nk);
ky=zeros(nk,1);
for i=1:1:nk
	ky(i)=calib(kgrid(i),par);
	i
    ky(i)
end

% selection of rows that are non-zero
x=horzcat(ky,kgrid');
y= x(x(:,1)>0,:);

kinit=y(1,2);

kss = fsolve(@(x)calib(x,par),kinit);
disp('calibrated steady state value of k: ');
kss 
pause;
gss=GY*a0*kss^(alpha);
bss=0;
par.gss = gss;
par.bss = bss;


% debt-output ratio
nb=75;
bygrid=linspace(0,0.0001*nb,nb);
kb=zeros(nb,1); % results for capital k and consumption c
cb=zeros(nb,1);
yb=cb;


BY=bygrid(1);
par.BY = BY;

kss = fsolve(@(x)ksteady(x,par),kinit);
disp('ksteady(kss) = 0?');
ksteady(kss,par)
pause;


for i=1:1:nb
	BY=bygrid(i);
    par.BY = BY;
    kss = fsolve(@(x)ksteady(x,par),kinit);
    i
    kss
	kb(i)=kss;
	kinit=kss;
	b0=BY*a0*kss^alpha;
	r=interest(kss,par);
	w=wage(kss,par);
	tr=(n+ga+ga*n-r)*b0-gss;
	cb(i)=(w+tr)/ (  1+ beta1^(1/eta)*(1+r)^(1/eta-1)  );
end

figure
plot(bygrid*3000,kb);
xlabel('Debt-output ratio B/Y');
ylabel('Capital stock k');


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% interest rate 
function [y] = interest(k,par)
    alpha = par.alpha;
    a0 = par.a0;
    delta = par.delta;
	y=alpha*a0*k^(alpha-1)-delta;
end

%% wage 
function [w] = wage(k,par)
    alpha = par.alpha;
	w=(1-alpha)*k^alpha;
end



% computes steady state for pension case and growth
function [y] = ksteady(x,par)
    
    a0 = par.a0;
    BY = par.BY;
    alpha = par.alpha;
    GY = par.GY;
    n = par.n;
    ga = par.ga;
    beta1 = par.beta1;
    eta = par.eta;
    gss = par.gss;

    k0 = x;
	b0=BY*a0*k0^alpha;
	r=interest(k0,par);
	
	tr=(n+ga+ga*n-r)*b0-gss;
	y=(1+n)*(1+ga)*(b0+k0)-(wage(k0,par)+tr)*(1-1/ ( 1+beta1^(1/eta)*(1+r)^(1/eta-1) ) );
end




% procedure that computes the steady state
% input: k0 - steady state capital stock 
% output: the value of the equilibrium condition
function [y] = calib(k0,par)
        
    alpha = par.alpha;
    BY = par.BY;
    GY = par.GY;
    a0 = par.a0;
    n = par.n;
    ga = par.ga;
    beta1 = par.beta1;
    eta = par.eta;
    
	b0=BY*a0*k0^alpha;
	r=interest(k0,par);
	g0=GY*a0*k0^alpha;
	tr=(n+ga+ga*n-r)*b0-g0;
	y=(1+n)*(1+ga)*(b0+k0)-(wage(k0,par)+tr)*(1-1/ ( 1+beta1^(1/eta)*(1+r)^(1/eta-1) ) );
end
