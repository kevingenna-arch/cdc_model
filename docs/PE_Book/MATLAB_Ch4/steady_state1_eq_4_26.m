function [y] = steady_state1_eq_4_6(par,x)
% computes the steady state for exogenous l=0.3
% and endogenous kappa,k,C^p
% in Chapter 4.2.1 in Heer (2018)
	
	y=zeros(3,1);
	k=x(1);
	kappa=x(2);
	cp=x(3);
	l=par.lbar;
    
    KL=(par.alpha/(1/par.beta1-1+par.delta))^(1/(1-par.alpha));
	w=wage_eq_4_9(k,l,par.alpha);
	r=interest_eq_4_9(k,l,par.alpha);
	y0=production_eq_4_8(k,l,par.alpha);
	g=par.GY*y0;
	
	zt=par.phi* cp^(1-1/par.rhoc) + (1-par.phi) * g^(1-1/par.rhoc);
	aggct=zt^(1/(1-1/par.rhoc));
	xt=aggct^(1-1/par.rho)+kappa*(1-l)^(1-1/par.rho);
	
	y(1)=cp+g+par.delta*k-y0;	% resource constraint
	y(2)=w*par.phi*((1-l)^(1/par.rho)) - kappa*aggct^(1/par.rho) *  zt^(1-1/ (1-1/par.rhoc)) * cp^(1/par.rhoc);
	y(3)=k/l-KL;
	

end
