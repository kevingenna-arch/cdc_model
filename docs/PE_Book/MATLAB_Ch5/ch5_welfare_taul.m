function Ch5_welfare_taul();
% 
%	Chapter 5.3 Labor Income Tax, transition dynamics
%   Author: Burkhard Heer
%   Last Change: June 21, 2018


%
%  Step 1: House Cleaning
%

clc;
close all;

%
% Step 2: Parameterization/Calibration 
%
par.alpha=0.36;
par.beta1=0.96;		% annual real interest rate of 4%
par.lstar=0.30;		% steady-state labor supply
par.sigma=2.0;		% 1/IES
par.taul=0.23;
%par.taul=0.40;     % tax wedge US from Prescott (2004)
%par.taul=0.59;     % tax wedge Germany
par.taul0=par.taul; % taul in initial steady state
par.delta=0.1;
par.nt=40;          % number of transition periods
par.display=0;      % set to 1 below so that the transition path is displayed in function Ch53_findk
par.utility=0;      % set to 1 below so that the welfare loss during transition is computed

rstar=1/par.beta1-1+par.delta;
klstar=(par.alpha/rstar)^(1/(1-par.alpha));
kstar=klstar*par.lstar;
par.kstar=kstar;
par.kfinal = par.kstar;
wstar=(1-par.alpha)*klstar^(par.alpha);
ystar=kstar^(par.alpha)*par.lstar^(1-par.alpha);

disp('parameterization');
par

gstar=par.taul*wstar*par.lstar;
disp('government share: ');
gstar/ystar
cstar=ystar-gstar-par.delta*kstar;
par.cstar = cstar;
par.gstar = gstar;

temp=(1-par.taul)*wstar*(1-par.lstar)/cstar;
gam=1/(1+temp);
par.gam=gam;
disp('calibration of gamma '); 
gam
pause;

xinitial= [par.lstar, cstar];
y=foc1(xinitial,par)
pause;

% Step 3: computation of partial equilibrium effect;
% computation of partial equilibrium effect;
% k and w constant
taulnew=par.taul+0.01;
par.taul1 = taulnew;
par.taul=taulnew;      

sol = fsolve(@(x)foc1(x,par),xinitial)
disp('solution : ');
disp('c and labor: ');
sol

lnew=sol(1);
cnew=sol(2);
disp('new equilibrium values: ');
cnew
lnew
disp('percentage changes: '); 
(cnew-cstar)/cstar
(lnew-par.lstar)/par.lstar

('consumption equivalent change: ');
ce0=(util(cstar,par.lstar,par)/ util(cnew,lnew,par) )^(-1/(par.gam*(1-par.sigma)))-1
pause;

%
% Step 4: general equlibrium effects: comparative steady states
par.taul=par.taul0;
% test: =0?";
xinitial=[par.lstar, cstar, par.kstar, wstar, rstar];
y=steadystate(xinitial,par)
pause;

par.taul=par.taul1;

disp('new allocation in general equilibrium');
disp('taul: ');
par.taul
xf = fsolve(@(x)steadystate(x,par),xinitial)


lnew1=xf(1);
cnew1=xf(2);
knew1=xf(3);
wnew1=xf(4);
rnew1=xf(5);

gnew1=par.taul*lnew1*wnew1;

disp('percentage changes c and l'); 
(cnew1-cstar)/cstar
(lnew1-par.lstar)/par.lstar

disp('consumption equivalent change: ');
ce1=(util(cstar,par.lstar,par)/ util(cnew1,lnew1,par) )^(-1/(par.gam*(1-par.sigma)))-1
pause;

%
% Step 5: Computation of the Transition
%
%      outer loop: over k_{nt} in the final period of the transition -> function findk()
%                  return: k_0 - kstar, the implied capital stock in period 0 minus 
%                                  the initial steady state value of capital
%
%      inner loop: within &findkinitial, we have to compute the dynamics
%                      Step 5.1: for given k_{nt}, we compute L_{nt} and C_{nt} in function findlcp()
%                      Step 5.2: Iteration over t=nt-1,...0
%                                  For given k_{t+1},L_{t+1},C_{t+1}, we have to find k_t,L_t,C_t 
%                                  in function findlcpk()
%
% Capital converges to the new steady state value from above, therefore k_{nt} has to lie between the old and
%      the new steady state value, kstar and knew1, and closer to the new steady state knew1

par.kfinal=knew1;
par.lfinal=lnew1;
par.cfinal=cnew1;
par.gfinal=taulnew*wage(knew1,lnew1,par)*lnew1;

kinitial=1.001*knew1;
y=Ch53_findk(kinitial,par)
pause;

clc;
k_solution = fsolve(@(k0)Ch53_findk(k0,par),kinitial);
clc;
disp('solution capital stock in period nt: ');
k_solution
pause;
par.display=1;
par.utility=1;
y=Ch53_findk(k_solution,par)



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% first-order condition labor 
function [y] = findlcpk(x,par)

	k=x(1);
	labor=x(2);
	c=x(3);
	k1=par.kold;
	labor1=par.lold;
	c1=par.cold;
	y=zeros(3,1);
	
	w0=wage(k,labor,par);
	r1=interest(k1,labor1,par);
	g0=par.taul*w0*labor;
	
	y(1)=w0*(1-par.taul)*(1-labor)-(1-par.gam)/par.gam*c;
	y(2)=(1-par.delta)*k+production(k,labor,par)-k1-c-g0;
	y(3)=(c1/c)^(1-par.gam*(1-par.sigma))*( (1-labor1)/(1-labor) )^(-(1-par.gam)*(1-par.sigma))-par.beta1*(1+r1-par.delta);

end


%% first-order condition labor 
function [y] = foc1(x,par)
	y=zeros(2,1);
    labor = x(1);
    c = x(2);
    k = par.kfinal;
    wstar = wage(k,labor,par);
    rstar = interest(k,labor,par);
	y(1)=(1-par.gam)/par.gam*c-(1-labor)*(1-par.taul)*wstar;
	y(2)=c+par.delta*k-(1-par.taul)*labor*wstar-rstar*k;
end

%% first-order condition labor 
function [y] = Ch53_findk(k0,par)
	kyoung=k0;      % last period nt of transition
	kold=par.kfinal;	% new steady state in period nt+1
	lold=par.lfinal;
	cold=par.cfinal;
	nt = par.nt;
	
	xtimepath=zeros(nt+2,6);    % time path for endogenous variables k,l,y,c,g 
	xtimepath(nt+2,1)=par.kfinal;
	xtimepath(nt+2,2)=par.lfinal;
	xtimepath(nt+2,3)=production(par.kfinal,par.lfinal,par);
	xtimepath(nt+2,4)=par.cfinal;
	xtimepath(nt+2,5)=par.gfinal;
	xtimepath(nt+2,6)=nt+1;
	
	xtimepath(1,1)=par.kstar;
	xtimepath(1,2)=par.lstar;
	xtimepath(1,3)=production(par.kstar,par.lstar,par);
	xtimepath(1,4)=par.cstar;
	xtimepath(1,5)=par.gstar;
	xtimepath(1,6)=0;
	
	% computation of l and cp in period nt
	x0=[lold, cold];
    y1 = foc1(x0,par);
    sol = fsolve(@(x)foc1(x,par),x0);
	lyoung=sol(1); 
	cyoung=sol(2);
	
	xtimepath(nt+1,1)=kyoung;
	xtimepath(nt+1,2)=lyoung;
	xtimepath(nt+1,3)=production(kyoung,lyoung,par);
	xtimepath(nt+1,4)=cyoung;
	xtimepath(nt+1,5)=par.taul*wage(kyoung,lyoung,par)*lyoung;
	xtimepath(nt+1,6)=nt;
	
	
	
	x0=[kyoung, lyoung, cyoung];
	x=x0;
	
	for i=nt-1:-1:1
		kold=x(1);
		lold=x(2);
		cold=x(3);
        par.kold=kold;
        par.lold=lold;
        par.cold=cold;
		
        x = fsolve(@(x)findlcpk(x,par),x0);
		x0=x;
		
		kyoung=x(1);
		lyoung=x(2); 
		cyoung=x(3);
	
		xtimepath(i+1,1)=kyoung;
		xtimepath(i+1,2)=lyoung;
		xtimepath(i+1,3)=production(kyoung,lyoung,par);
		xtimepath(i+1,4)=cyoung;
		xtimepath(i+1,5)=par.taul*wage(kyoung,lyoung,par)*lyoung;
		xtimepath(i+1,6)=i;
		
		
    end
	
    if par.display==1
        figure
        subplot(3,2,1);
        plot(xtimepath(:,6),xtimepath(:,1))
        title('Capital k'); 
        xlabel('Period');
        subplot(3,2,2);
        plot(xtimepath(:,6),xtimepath(:,2)) 
        title('Labor l'); 
        xlabel('Period');
        subplot(3,2,3);
        plot(xtimepath(:,6),xtimepath(:,3))
        title('Production y'); 
        xlabel('Period');
        subplot(3,2,4);
        plot(xtimepath(:,6),xtimepath(:,4))
        title('Consumption c'); 
        xlabel('Period');
        subplot(3,2,5);
        plot(xtimepath(:,6),xtimepath(:,5))
        title('Tax revenue'); 
        xlabel('Period');
        pause;
    end
    
% compuation of welfare loss during transition
    if par.utility==1
        utillife0=util(par.cstar,par.lstar,par)/(1-par.beta1);
        utillife1=0;
        for i=1:1:nt
            c0=xtimepath(i+1,4);
            l0=xtimepath(i+1,2);
            utillife1=utillife1+par.beta1^(i-1)*util(c0,l0,par);
        end
		utilnew=util(par.cfinal,par.lfinal,par);
		utillife1=utillife1+par.beta1^(nt)*utilnew/(1-par.beta1);

        disp('consumption equivalent change: ');
        ce2=(utillife0/ utillife1 )^(-1/(par.gam*(1-par.sigma)))-1;
        ce2*100  
        pause;
    end
    
	y=kyoung-par.kstar;
end

%% utility function
function [y] = util(c,labor,par)
	y=c^(par.gam)*(1-labor)^(1-par.gam);
	y=y^(1-par.sigma)/(1-par.sigma);
end

%% steady state
function [y] = steadystate(x,par)
	
    labor = x(1);
	c=x(2);
	k=x(3);
	w=x(4);
	r=x(5);
	y=zeros(5,1);
	y(1)=(1-par.gam)/par.gam*c-(1-labor)*(1-par.taul)*w;
	y(2)=c+par.delta*k-(1-par.taul)*labor*w-r*k;
	y(3)=r-interest(k,labor,par);
	y(4)=w-wage(k,labor,par);
	y(5)=1/par.beta1-1+par.delta-r;
	
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



