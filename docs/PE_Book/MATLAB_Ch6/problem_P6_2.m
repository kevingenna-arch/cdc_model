function problem_P6_2();
% 
%   Problem 6.2 in Heer, Public Economics. The Macroeconomic Perspective,
%   2019, Springer
%
%   Author: Burkhard Heer
%   Last Change: May 8, 2020
%
%   See also the solutions material on my web page:
%   https://assets.uni-augsburg.de/media/filer_public/5b/64/5b648da2-25f5-43f1-bea7-b1a8b8748084/book_heer_solutions_to_problems.pdf
%

clc;
close all;

%
% Step 1: Parameterization 
%



ls=0.3;			% labor supply
alpha=0.36;		% production elasticity of capital
beta0=0.90;		% discount factor 
sigma0=2.0;		% 1/IES
iota=2.0;       % weight of leisure in utility
delta=0.5;      % deprecation rate
tau=0.10;       % initial social security tax
rhopen=0.5;     % pension schedule
penmin=0;       % needs to be calibrated

par.ls = ls;
par.alpha = alpha;
par.beta0 = beta0;
par.iota = iota;
par.sigma0 = sigma0;
par.rhopen = rhopen;
par.tau = tau;
par.delta=delta;

% educative guess:
%
% 1. marginal utility the same at age 1 and 2
% => lambda^1 = lambda^2
% 1/beta = 1+r
%
% 2. l=0.3 => L=1/3 * (0.3+0.3) = 0.2;
rbar=1/beta0-1;
lbar=0.2;
l1bar=0.3;
l2bar=0.3;
kbar=(alpha/(rbar+delta))^(1/(1-alpha))*lbar
wbar=wage(kbar,lbar,par);
penmin=(1-rhopen)*tau*wbar*(l1bar+l2bar);
% educative guess for k^2 and k^3
% k increases until retirement, therefore, k^3=2 k^2
k3bar=2*kbar;
k2bar=kbar;
c1bar=(1-tau)*wbar*l1bar-k2bar;
c2bar=(1-tau)*wbar*l2bar+(1+rbar)*k2bar-k3bar;
c3bar=(1+rbar)*k3bar+penmin+rhopen*tau*wbar*(l1bar+l2bar);


xinitial=[c1bar,c2bar,c3bar,l1bar,l2bar,lbar,k2bar,k3bar,kbar,penmin];

% compute a grid for the solution for various k
% eq. (6.19)
%kmin=0.001;
%kmax=5.0;
%nk=1000;
%k=linspace(kmin,kmax,nk);
%yk=zeros(nk,1);


clc;
%kss = fsolve(@(x)ksteady(x,par),kinitial);

disp('steadystate(xinitial,par) ');
steadystate(xinitial,par)

 yss = fsolve(@(x)steadystate(x,par),xinitial)
steadystate(yss,par)
disp('Did fsolve find a solution to steadystate? No.');
disp('equations are not all equal to zero.');
pause

%
% solution for the steady state not found
%
% therefore: tatonnement over K,L
%


update=0.9;
 
xinitial = [c1bar, c2bar, c3bar, l1bar, l2bar, k2bar, k3bar, penmin];

 for i=1:1:100
    par.kbar = kbar;
    par.lbar = lbar;
    yss = fsolve(@(x)foc(x,par),xinitial)
    xinitial=yss;
    % update of K and L
    kbarnew = 1/3*(yss(6)+yss(7))
    kbar=update*kbar + (1-update)*kbarnew
    lbarnew = 1/3*(yss(4)+yss(5))
    lbar = update*lbar+(1-update)*lbarnew

end


foc(yss,par)
display('solution found? yes. All eqs=0.')
pause


zinitial=[yss(1:5), lbar, yss(6:7), kbar, yss(8)]
zss=steadystate(zinitial,par);

    solution = fsolve(@(x)steadystate(x,par),zinitial)
    
    zss=steadystate(solution,par)
"solution vector for steady state and value of non-linear eqs: ";    
yssinitial=solution;
pause

x=solution;
    c1=x(1);
    c2=x(2);
    c3=x(3);
    l1=x(4);
    l2=x(5);
    bigl=x(6);
    k2=x(7);
    k3=x(8);
    bigk=x(9);
    penmin=x(10);

w=wage(bigk,bigl,par);
disp('pension replacement rate: ');
pen_repl=(penmin+rhopen*tau*w*(l1+l2))/(w*l2);
pen_repl
par.pen_repl=pen_repl;
disp('benchmark values: ');
production=bigk^alpha*bigl^(1-alpha);
production
bigl
bigk
Ubench=lifetimeutil(c1,c2,c3,l1,l2,par);
disp('lifetime utility U: '); 
Ubench
penbench=penmin+rhopen*tau*w*(l1+l2);
disp('pension level in benchmark: ');
penbench
pause


tau=0;
par.tau=tau;

    solution0 = fsolve(@(x)steadystate(x,par),solution)
    
    zss=steadystate(solution0,par)
"solution vector for steady state, tau=0, and value of non-linear eqs: ";    
pause

x=solution0;
    c10=x(1);
    c20=x(2);
    c30=x(3);
    l10=x(4);
    l20=x(5);
    bigl0=x(6);
    k20=x(7);
    k30=x(8);
    bigk0=x(9);
    penmin0=x(10);

disp('aggregate variables with tau=0: ');
production=bigk0^alpha*bigl0^(1-alpha);
production
bigl0
bigk0
U0=lifetimeutil(c10,c20,c30,l10,l20,par);
disp('lifetime utility U: ');
U0
disp('consumption equivalent change: ');
delta0=(U0/Ubench)^(1/(1-sigma0))-1;
delta0
pause


% defined benefits
% pension are porportional to pre-retirement earnings, 
% however, pensions do not depend on contributions
tau=0.10;   % initial guess, tau is now endogenous
solutiondb = [solution(1:9), tau];
y=steadystatedb(solutiondb,par);  % tests steadystatedb



    solution1 = fsolve(@(x)steadystatedb(x,par),solutiondb)
    
    zss=steadystatedb(solution1,par)
"solution vector for steady state, tau=10, defined contributions and value of non-linear eqs: ";    
pause


x=solution1;
    c11=x(1);
    c21=x(2);
    c31=x(3);
    l11=x(4);
    l21=x(5);
    bigl1=x(6);
    k21=x(7);
    k31=x(8);
    bigk1=x(9);
    tau1=x(10);

disp('aggregate variables with defined benefits: ');
production=bigk1^alpha*bigl1^(1-alpha);
production
bigl1
bigk1
U1=lifetimeutil(c11,c21,c31,l11,l21,par);
disp('lifetime utility U: '); 
U1
disp('consumption equivalent change: ');
delta1=(U1/Ubench)^(1/(1-sigma0))-1;
delta1
disp('new social security tax tau: '); 
tau1
pen1=pen_repl*wage(bigk1,bigl1,par)*l21;
disp('pension level in this case: ');
pen1


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

%% steady state in the case of defined contributions 
function [y] = steadystate(x,par)

    
	sigma0 = par.sigma0;
    beta0 = par.beta0;
    alpha = par.alpha;
    iota = par.iota;
    rhopen = par.rhopen;
    tau = par.tau;
    
    c1=x(1);
    c2=x(2);
    c3=x(3);
    l1=x(4);
    l2=x(5);
    bigl=x(6);
    k2=x(7);
    k3=x(8);
    bigk=x(9);
    penmin=x(10);
    
    
	w=wage(bigk,bigl,par);
	r=interest(bigk,bigl,par);
	
    lambda1 = c1^(-sigma0)*(1-l1)^(iota*(1-sigma0));
    lambda2 = c2^(-sigma0)*(1-l2)^(iota*(1-sigma0));
    lambda3 = c3^(-sigma0);
    
    y=zeros(10,1);
    
    % foc labor
    y(1) = iota*c1^(1-sigma0)*(1-l1)^(iota*(1-sigma0)-1)-lambda1*(1-tau)*w-beta0^2*lambda3*rhopen*tau*w;
    y(2) = iota*c2^(1-sigma0)*(1-l2)^(iota*(1-sigma0)-1)-lambda2*(1-tau)*w-beta0^2*lambda3*rhopen*tau*w;
    % Euler eqs.
    y(3) = lambda1-beta0*lambda2*(1+r);
    y(4) = lambda2-beta0*lambda3*(1+r);
	%  budget constraints
    y(5) = (1-tau)*w*l1-k2-c1;
    y(6) = (1-tau)*w*l2+(1+r)*k2-k3-c2;
    y(7) = (1+r)*k3+penmin+rhopen*tau*w*(l1+l2)-c3;
    % aggregate consistency conditions
    y(8) = bigk-1/3*k2-1/3*k3;
    y(9) = bigl-1/3*l1-1/3*l2;
    y(10) = penmin-(1-rhopen)*tau*w*(l1+l2);
	
end


% procedure that computes the equilibrium for given K and L
%
    function[y] = foc(x,par)
    
    
	sigma0 = par.sigma0;
    beta0 = par.beta0;
    alpha = par.alpha;
    iota = par.iota;
    rhopen = par.rhopen;
    bigk = par.kbar;
    bigl = par.lbar;
    tau = par.tau;
    
	w=wage(bigk,bigl,par);
	r=interest(bigk,bigl,par);
    
    c1=x(1);
    c2=x(2);
    c3=x(3);
    l1=x(4);
    l2=x(5);
    k2=x(6);
    k3=x(7);
    penmin=x(8);
    
    
    
	
    lambda1 = c1^(-sigma0)*(1-l1)^(iota*(1-sigma0));
    lambda2 = c2^(-sigma0)*(1-l2)^(iota*(1-sigma0));
    lambda3 = c3^(-sigma0);
    
    y=zeros(8,1);
    
    % foc labor
    y(1) = iota*c1^(1-sigma0)*(1-l1)^(iota*(1-sigma0)-1)-lambda1*(1-tau)*w-beta0^2*lambda3*rhopen*tau*w;
    y(2) = iota*c2^(1-sigma0)*(1-l2)^(iota*(1-sigma0)-1)-lambda2*(1-tau)*w-beta0^2*lambda3*rhopen*tau*w;
    % Euler eqs.
    y(3) = lambda1-beta0*lambda2*(1+r);
    y(4) = lambda2-beta0*lambda3*(1+r);
	%  budget constraints
    y(5) =(1-tau)*w*l1-k2-c1;
    y(6) =(1-tau)*w*l2+(1+r)*k2-k3-c2;
    y(7) = (1+r)*k3+penmin+rhopen*tau*w*(l1+l2)-c3;
    % aggregate consistency conditions
    y(8) = penmin-(1-rhopen)*tau*w*(l1+l2);
	
end

% procedure that computes the steady state with defined benefits
function[y] = steadystatedb(x,par)

	sigma0 = par.sigma0;
    beta0 = par.beta0;
    alpha = par.alpha;
    iota = par.iota;
    pen_repl = par.pen_repl;
    

    c1=x(1);
    c2=x(2);
    c3=x(3);
    l1=x(4);
    l2=x(5);
    bigl = x(6);
    k2=x(7);
    k3=x(8);
    bigk=x(9);
    tau=x(10);
    
    
	w=wage(bigk,bigl,par);
	r=interest(bigk,bigl,par);
	
    pen=pen_repl*w*l2;  % defined benefits, exogenous to workers
    
    lambda1 = c1^(-sigma0)*(1-l1)^(iota*(1-sigma0));
    lambda2 = c2^(-sigma0)*(1-l2)^(iota*(1-sigma0));
    lambda3 = c3^(-sigma0);
    
    y=zeros(10,1);
    
    % foc labor
    y(1) = iota*c1^(1-sigma0)*(1-l1)^(iota*(1-sigma0)-1)-lambda1*(1-tau)*w; % needs to be adjusted because of defined benefits
    y(2) = iota*c2^(1-sigma0)*(1-l2)^(iota*(1-sigma0)-1)-lambda2*(1-tau)*w;
    % Euler eqs.
    y(3) = lambda1-beta0*lambda2*(1+r);
    y(4) = lambda2-beta0*lambda3*(1+r);
	%  budget constraints
    y(5) =(1-tau)*w*l1-k2-c1;
    y(6) =(1-tau)*w*l2+(1+r)*k2-k3-c2;
    y(7) = (1+r)*k3+pen-c3;
    % aggregate consistency conditions
    y(8) = bigk-1/3*k2-1/3*k3;
    y(9) = bigl-1/3*l1-1/3*l2;
    y(10) = 1/3*pen-tau*w*bigl;
	
end

% lifetime utility
function[y] = lifetimeutil(c1,c2,c3,l1,l2,par)
    beta0 = par.beta0;
    y=util(c1,l1,par)+beta0*util(c2,l2,par)+beta0^2*util(c3,0,par);
end

function[y] = util(c,l,par)

    beta0 = par.beta0;
    iota = par.iota;
    sigma0 = par.sigma0;
    
    y=(c^(1-sigma0)*(1-l)^(iota*(1-sigma0)))/(1-sigma0);
    
end
