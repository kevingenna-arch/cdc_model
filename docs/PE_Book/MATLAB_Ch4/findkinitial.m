function [ y ] = findkinitial( par2, k0)
% computes the transition dynamics
% initial value problem k0
% eq. (4.75) in Heer (2018)

    kyoung=k0;      % last period nt of transition
	kold=par2.kfinal;	% new steady state in period nt+1
	lold=par2.lfinal;
	cpold=par2.cpfinal;
	govt = par2.govt;
	nt = par2.nt;

	xtimepath=zeros(par2.nt+2,7);   % time path for endogenous variables k,l,y,cp,c,g
	
    % final steady state after transition is complete
    xtimepath(nt+2,1) = par2.kfinal;
	xtimepath(nt+2,2) = par2.lfinal;
	xtimepath(nt+2,3) = production_eq_4_8(kold,lold,par2.alpha);
	xtimepath(nt+2,4) = par2.cpfinal;
	xtimepath(nt+2,5) = par2.aggctfinal;
	xtimepath(nt+2,6) = par2.gfinal;
	xtimepath(nt+2,7) = par2.nt+1;
	
    % initial steady state in period t=-1
	xtimepath(1,1) = par2.kbar;
	xtimepath(1,2) = par2.lbar;
	xtimepath(1,3) = production_eq_4_8(par2.kbar,par2.lbar,par2.alpha);
	xtimepath(1,4) = par2.cpbar;
	xtimepath(1,5) = par2.aggctbar;
	xtimepath(1,6) = par2.gbar;
	xtimepath(1,7) = 0;
	
	
	% computation of l and cp in period nt
	x0 = [lold, cpold];
    par2.k0 = k0;

    y0 = findlcp(par2,x0)
    x1 = fsolve(@(x)findlcp(par2,x),x0);
	lyoung = x1(1); 
	cpyoung =x1(2);
	
	xtimepath(nt+1,1) = kyoung;
	xtimepath(nt+1,2) = lyoung;
	xtimepath(nt+1,3) = production_eq_4_8(kyoung,lyoung,par2.alpha);
	xtimepath(nt+1,4) = cpyoung;
	xtimepath(nt+1,5) = ( par2.phi*cpyoung^(1-1/par2.rhoc)+(1-par2.phi)*(govt(nt))^(1-1/par2.rhoc) )^(1/(1-1/par2.rhoc));
	xtimepath(nt+1,6) = govt(nt);
	xtimepath(nt+1,7) = nt;
	
	
	
	x0=[kyoung, lyoung, cpyoung];
	x=x0;
	
	
	% computation of k, l and cp in period nt-1, nt-2, ..., 0
	
    for i=nt-1:-1:1
		kold = x(1);
		lold = x(2);
		cpold = x(3);
        
        par2.kold = kold;
        par2.lold = lold;
        par2.cpold = cpold;
        par2.it = i;
		
        x1 = fsolve(@(x)findlcpk(par2,x),x0);
		x0=x1;
        x=x1;
		
		kyoung = x1(1);
		lyoung = x1(2); 
		cpyoung = x1(3);
	
		xtimepath(i+1,1) = kyoung;
		xtimepath(i+1,2) = lyoung;
		xtimepath(i+1,3) = production_eq_4_8(kyoung,lyoung,par2.alpha);;
		xtimepath(i+1,4) = cpyoung;
		xtimepath(i+1,5) = ( par2.phi*cpyoung^(1-1/par2.rhoc)+(1-par2.phi)*govt(i)^(1-1/par2.rhoc) )^(1/(1-1/par2.rhoc));
		xtimepath(i+1,6) = govt(i);
		xtimepath(i+1,7) = i;
		
		
    end
	
    par2.xtimepath = xtimepath;
	
	y= x(1)- par2.kbar;


end

