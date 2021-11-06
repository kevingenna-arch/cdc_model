% ----------------------------  Problem_P3_5.m ----------------------------
% computes the solution to Problem 3.5 in Heer (2019), Public Economics
%
% builds upon Ch3_turnpike.m
%
% computes the dynamics of the central planner OLG economy using
%
% 1. solving a non-linear eqs system in (k_1,...,k_20) simulataneously
% 2. reverse shooting: NEW,
%    starting with a guess for k_20 given initial and final values k_0 and k_21
%
% author: Burkhard Heer
%
% this version: May 5, 2020
%
% -------------------------------------------------------------------------



clear; clc;
%parameters


par.alpha=0.36;		% production elasticity of capital
par.beta=0.4;			%  discount factor 

par.n=0.1;			% population growth rate

k=(par.alpha/par.n)^(1/(1-par.alpha));	% steady state value of capital stock 

par.k0=0.7*k;		% initial capital stock 

par.kfinal=0.7*k;   % final capital stock

par.bigt=20;        % number of transition periods

tt=0:1:par.bigt+2;	% periods

kt=zeros(par.bigt+2,1);		% time series for capital stock

x=kdyn(par,k);

kguess=1.1*par.k0;          % initial guess for solution

ksol = fsolve(@(k)kdyn(par,k),kguess);


    kt2=zeros(par.bigt+2,1);
    kt2(1)=par.k0;
	kt2(2)=ksol;
	for i=3:par.bigt+2,
		kt2(i)=kt2(i-2)+kt2(i-2)^(par.alpha) -(1+par.n)*kt2(i-1) ;
        kt2(i)=(kt2(i-1)+kt2(i-1)^(par.alpha))/(1+par.n) - (1+par.alpha*kt2(i-1)^(par.alpha-1))/((1+par.n)^2) * kt2(i);
    end;


figure
plot([0:par.bigt+1],kt2); hold on
xlabel('Period t')
ylabel('Capital stock k_t')
title('Solution found with forward shooting')

%
%   Part II: Reverse Shooting
%



k20initial=par.kfinal*1.01;        % initial guess for k_20

kdyn_reverse(par,k20initial)

k20sol = fsolve(@(k)kdyn_reverse(par,k),k20initial)

    kt3=zeros(par.bigt+2,1);
	kt3(par.bigt+2)=par.kfinal;
    kt3(par.bigt+1)=k20sol;
	for i=par.bigt:-1:1,
        
        par.k1=kt3(i+1);
        par.k2=kt3(i+2);
        % non-linear equation in k_t given k_t+1 and k_t+2
        kguess=kt3(i+1);
        ksol = fsolve(@(k)kdyn1(par,k),kguess);
        
        kt3(i)=ksol;
    end;

figure
plot([0:par.bigt+1],kt3); hold on
xlabel('Period t')
ylabel('Capital stock k_t')
title('Solution found with reverse shooting')

function [y] = kdyn_reverse(par,x)

    alpha1=par.alpha;
    n=par.n;
    
    kt2=zeros(par.bigt+2,1);
    kt2(1)=par.k0;
	kt2(par.bigt+2)=par.kfinal;
    kt2(par.bigt+1)=x;
	for i=par.bigt:-1:1,
        
        par.k1=kt2(i+1);
        par.k2=kt2(i+2);
        % non-linear equation in k_t given k_t+1 and k_t+2
        kguess=kt2(i+1);
        ksol = fsolve(@(k)kdyn1(par,k),kguess);
        
        kt2(i)=ksol;
    end;
	y=kt2(1)-par.k0; 
end

function [y] = kdyn1(par,x)
% Dynamics of the capital stock
% Equation (3.35) in Chapter 3.3.4 in Heer (2018)
	
    alpha1=par.alpha;
    n=par.n;
    k0=x;
    k1=par.k1;
    k2=par.k2;
    y=k2-(k1+k1^(alpha1))/(1+n)+(1+alpha1*k1^(alpha1-1))/((1+n)^2)*(k0+k0^(alpha1)-(1+n)*k1);
end

function [kdiff] = kdyn(par,x)
% Dynamics of the capital stock
% Equation (3.38) in Chapter 3.3.4 in Heer (2018)
	
    alpha1=par.alpha;
    n=par.n;
    
    kt2=zeros(par.bigt+2,1);
    kt2(1)=par.k0;
	kt2(2)=x;
	for i=3:par.bigt+2,
		kt2(i)=kt2(i-2)+kt2(i-2)^(alpha1) -(1+n)*kt2(i-1) ;
        kt2(i)=(kt2(i-1)+kt2(i-1)^(alpha1))/(1+n) - (1+alpha1*kt2(i-1)^(alpha1-1))/((1+n)^2) * kt2(i);
    end;
	kdiff=kt2(par.bigt+2)-par.kfinal; 

end
