function [ y ] = findlcpk( par1, x )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%   input: 
%   x -- initial value for labor, capital and consumption
%   par1 -- initial parameters par plus next-period values of k, xi, cp
%           and period i
%   output: y -- eq. (4.75)

    y = zeros(3,1);
    k = x(1);
    labor = x(2);
	cp = x(3);
    
    
	k1 = par1.kold;
	labor1 = par1.lold;
	cp1 = par1.cpold;
    govt = par1.govt;
    it = par1.it;     % period
	
	zt=par1.phi*cp^(1-1/par1.rhoc)+(1-par1.phi)*govt(it)^(1-1/par1.rhoc);
	aggct=zt^(1/(1-1/par1.rhoc));
	xt=aggct^(1-1/par1.rho)+par1.kappa*(1-labor)^(1-1/par1.rho);
	zt1=par1.phi*cp1^(1-1/par1.rhoc)+(1-par1.phi)*govt(it+1)^(1-1/par1.rhoc);
	aggct1=zt1^(1/(1-1/par1.rhoc));
	xt1=aggct1^(1-1/par1.rho)+par1.kappa*(1-labor1)^(1-1/par1.rho);
	w0=wage_eq_4_9(k,labor,par1.alpha);
	r1=interest_eq_4_9(k1,labor1,par1.alpha);
    y0=production_eq_4_8(k,labor,par1.alpha);

    it
    govt(it)
    
	y(1) = w0/par1.kappa*par1.phi-((1-labor)/aggct)^(-1/par1.rho)* zt^( 1-1/(1-1/par1.rhoc)) * cp^(1/par1.rhoc); 
    y(2) = (xt1/xt)^(1-(1-par1.sigma)/(1-1/par1.rhoc))*(aggct1/aggct)^(1/par1.rho) * (zt1/zt)^(1-1/(1-1/par1.rhoc))*(cp1/cp)^(1/par1.rhoc);
	y(2) = y(2)-(1+r1-par1.delta)*par1.beta1;
	y(3) = k1-(1-par1.delta)*k-y0+govt(it)+cp;

end

