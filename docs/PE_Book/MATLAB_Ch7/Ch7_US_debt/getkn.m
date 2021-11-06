function [ ktnew,ntnew,tauntnew,tauptnew,pentnew,ksharetnew,ctnew,trtnew,btnew ] = getkn()
% computes the new transition path of {K,N,tau^n,tau^p,pen} 
% from individual decision
%

%   local ktnew,ntnew,tauntnew,tauptnew,pentnew,ctnew,ksharetnew,meanlabtnew,x0,xf,jcode,asset,labors,c,i,nage0;
%	local taultnew, trtnew, btnew;
%	local tp, tax, labortax, beqtnew, tp0;
%	local income0t,income1t;
%	local wbar1, dbar1, taun1, taup1, pen1;
%	local tr1, debt1;

    def_global_USdebt
    def_global_USdebt_transition
    
    tt=nt+1; 
    ktnew=zeros(nt,1);
    ntnew=zeros(nt,1);
    tauntnew=zeros(nt,1);
    tauptnew=zeros(nt,1);
	taultnew=zeros(nt,1);
    pentnew=zeros(nt,1);
	ctnew=zeros(nt,1);
	ksharetnew=zeros(nt,1);
	meanlabtnew=zeros(nt,1);
	beqtnew=zeros(nt,1);
	trtnew=zeros(nt,1);
	btnew=zeros(nt,1);
	
    % computation of the optimal allocation of an agent born in period tt */
	nage0=nage-1;
	x0=xstartfinal(1:nage+nw-1);
	x0=xagginitial(1:nage+nw-1);
	wbar1=wt(nt);
	dbar1=dbart(nt);
	pen1=pent(nt);
	taun1=tauntold(nt);
	taup1=taupt(nt);
	tr1=trtold(nt);
	debt1=btold(nt);
	
    for tt=nt:-1:-nage0+1
		disp('tt~q~kritt~ (ndebt):  ');
        if policy<=1
            [tt, q, kritt]
        else
            [tt, q, kritt, ndebt]
        end
        % wage, interest rate and pensions over the lifetime of the individual 
        if tt>nt-nage0      % agent is alive after period nt 
            wseq(1:nt-tt+1) = wt(tt:nt);
            wseq(nt-tt+2:nage) = ones(nage-nt+tt-1,1)*wbar1;    % economy in steady state after period nt
            dseq(1:nt-tt+1) = dbart(tt:nt);
            dseq(nt-tt+2:nage) = ones(nage-nt+tt-1,1)*dbar1;
            penseq(1:nt-tt+1) = pent(tt:nt);
            penseq(nt-tt+2:nage) = ones(nage-nt+tt-1,1)*pen1;
            taunseq(1:nt-tt+1) = tauntold(tt:nt);
            taunseq(nt-tt+2:nage) = ones(nage-nt+tt-1,1)*tauntold(nt);
            taupseq(1:nt-tt+1) = taupt(tt:nt);
            taupseq(nt-tt+2:nage) = ones(nage-nt+tt-1,1)*taup1;
            trseq(1:nt-tt+1) = trtold(tt:nt);
            trseq(nt-tt+2:nage) = ones(nage-nt+tt-1,1)*tr1;
        elseif tt<1         % agent is born prior to transition 
            wseq(1:1-tt) = ones(-tt+1,1)*wbar0;
            wseq(2-tt:nage) = wt(1:nage0+tt);
            dseq(1:1-tt) = ones(-tt+1,1)*dbar0;
            dseq(2-tt:nage) = dbart(1:nage0+tt);
            penseq(1:1-tt) = ones(-tt+1,1)*pen0;
            penseq(2-tt:nage) = pent(1:nage0+tt);
            taunseq(1:1-tt) = ones(-tt+1,1)*tauntold(1);
            taunseq(2-tt:nage) = tauntold(1:nage0+tt);
            taupseq(1:1-tt) = ones(-tt+1,1)*taup0;
            taupseq(2-tt:nage) = tauptold(1:nage0+tt);
            trseq(1:1-tt) = ones(-tt+1,1)*trbar0;
            trseq(2-tt:nage) = trtold(1:nage0+tt);
        else                % agent is born and dies during transition
            wseq(1:nage) = wt(tt:tt+nage0);
            dseq(1:nage) = dbart(tt:tt+nage0);
            penseq(1:nage) = pent(tt:tt+nage0);
            taunseq(1:nage) = tauntold(tt:tt+nage0);
            taupseq(1:nage) = taupt(tt:tt+nage0);
            trseq(1:nage) = trtold(tt:tt+nage0);
        end

		taulseq=taunseq-taupseq;
	
		
        % computation of allocation of the individual born in period tt 
        [xf, Fval] = fsolve(fhandle_function1,x0);
		x0=xf;
        
		%if maxc(rftr(xf))==miss(1,1) or maxc(rftr(xf))==miss(1,1)>0.001;
		%	x0=xstartfinal[1:nage+nw-1];
		%	{xf,jcode}=FixVMN1(x0,&rftr);
		%	"maxc(rftr(xf)): " maxc(rftr(xf));
		%endif;
		if max(abs(Fval))>0.0001 
            disp('could not find solution for rftr');
            pause; 
        end
			
		asset=[0; xf(1:nage-1); 0];		% assets of the cohorts
		labors=xf(nage:nage+nw-1);
	
		% consumption
		c=zeros(nage,1);	
		income0t=zeros(nage,1);
		income1t=zeros(nage,1);
		
		
		for i=1:1:nw	
			c(i) = (1-taulseq(i)-taupseq(i)) * wseq(i)*ef(i)*labors(i);
            c(i) = c(i)	+(1+(1-tauk)*(dseq(i)-delta))*asset(i)-ygrowth*asset(i+1)+trseq(i);
			c(i) = c(i)/(1+tauc);
            if i+tt<2       % consumption pre-determined in period t=0,-1,-2
				c(i)=cinitial(i);
            end
			income0t(i) = (1-taulseq(i) -taupseq(i)) * wseq(i)*ef(i)*labors(i)+(1-tauk)*(dseq(i)-delta)*asset(i)+trseq(i);
        end

	
		for i=1:1:nr-1
			c(i+nw) = penseq(i+nw)+(1+(1-tauk)*(dseq(i+nw)-delta))*asset(nw+i)-ygrowth*asset(nw+i+1)+trseq(nw+i);
			c(i+nw) = c(i+nw)/(1+tauc);
			% consumption predetermined
			if i+nw+tt<2
				c(i+nw) = cinitial(i+nw);
            end		
			income0t(i+nw) = penseq(i+nw)+(1-tauk)*(dseq(i+nw)-delta)*asset(nw+i)+trseq(nw+i);
			income1t(i+nw) = penseq(i+nw);
        end
        
		c(nr+nw) = penseq(nr+nw)+( 1 + (1-tauk)*(dseq(nr+nw)-delta))*asset(nw+nr)+trseq(nw+nr);
		c(nr+nw) = c(nr+nw)/(1+tauc);
		% consumption predetermined
		if nr+nw+tt<2
			c(nr+nw) = cinitial(nr+nw);
        end
        
		income0t(nr+nw) = penseq(nr+nw)+(1-tauk)*(dseq(nr+nw)-delta)*asset(nw+nr)+trseq(nw+nr);
		income1t(nr+nw) = penseq(nr+nw);

		utilityt(tt+nage0,1) = utility(c,labors);   % life-time utility of agent born in tt
		utilityt(tt+nage0,2) = tt;
		consumprofilet(1:nage,tt+nage0) = c;
		laborprofilet(1:nw,tt+nage0) = labors;
		incomeprofilet(1:nage,tt+nage0) = income0t;
		income1profilet(1:nage,tt+nage0) = income1t;
		assetprofilet(1:nage,tt+nage0) = asset(1:nage);
		

		% Aggregation of K, N, C, and Beq
        tp=tt-1; i=0;
        while (tp<nt) && (i<nage)
            tp=tp+1; 
            i=i+1;
			if (tp>0) && (tp<nt)
				if (i>1) && (i<nage)
                    % choice of survival probability
                    % if life period is prior to 2010: sp from 2010
                    tp0=max(tp,1);
                    % if life period is after the final year: sp from yearfinal
                    tp0=min(tp,nsp);
					beqtnew(tp+1) = beqtnew(tp+1)+massvec(i,tp)*(1-spvec(i,tp0))*(1+(1-tauk)*(dseq(i+1)-delta))*asset(i+1);
                end
            end
			
            if tp>0
                ktnew(tp) = ktnew(tp) + massvec(i,tp)*asset(i);
                ctnew(tp) = ctnew(tp) + massvec(i,tp)*c(i);
                if i<=nw
                    ntnew(tp) = ntnew(tp) + massvec(i,tp) * labors(i) * ef(i);
					meanlabtnew(tp) = meanlabtnew(tp) + massvec(i,tp) * labors(i);
                end
            end
        end     % loop over age i



    end     % loop over period tt=nt,nt-1,..,-nage0+1

	ktnew=ksharetold.*ktnew;
	beqtnew(1)=bequestsinitial;
    
	% Dynamics of Debt, Pensions: Fiscal and Social Security Budgets	
	for tp=1:1:nt
		tp0 = min(tp,nsp);
		taxes = tauk*(dbart(tp)-delta)*ktnew(tp)+tauc*ctnew(tp);
				
        % policy:
		% ========
		% 0 -- benchmark, transfers are constant and labor income tax adjusts
		% 1 -- transfers adjust, 
		% 2 -- during the first ndebt periods, debt adjusts. Afterwards, labor income tax rate
		% 3 -- during the first ndebt periods, debt adjusts. Afterwards, transfers
				
		btnew(1) = bigb;		
		if (policy==1) || (policy==3)   % transfers adjust to balance budget
			taxes = taxes+(tauntold(tp)-tauptold(tp)) * wt(tp) * ntnew(tp);
			taultnew(tp)=taul0		% labor income tax constant
			if policy==1            % transfers adjust, constant debt
				trtnew(tp) = taxes-bigg+((1+popgrowthvec(tp0)*ygrowth)-(1+(1-tauk)*(dbart(tp)-delta)))*btnew(tp)+beqtnew(tp);
				btnew(tp) = bigb;   % constant debt under policy 1
			elseif policy==3		
				if tp<=ndebt        % debt increases during the first ndebt periods, constant transfers
					trtnew(tp) = trbar0;
					btnew(tp+1)=(1+(1-tauk)*(dbart(tp)-delta))*btold(tp)+bigg+trtnew(tp)-beqtnew(tp)-taxes;
					btnew(tp+1)=btnew(tp+1) / (1+popgrowthvec(tp0)*ygrowth);    % per capita debt
                else
					trtnew(tp) = taxes-bigg+((1+popgrowthvec(tp0)*ygrowth)-(1+(1-tauk)*(dbart(tp)-delta)))*btold(tp)+beqtnew(tp);
					if tp>ndebt+1
						btnew(tp)=btnew(ndebt+1);
                    end
                end
            end
        end

		if (policy==0) || (policy==2)   % constant transfers, labor income taxes adjust
			trtnew(tp) = trbar0;
			
			if policy==0
				labortax = bigg+bigb*(1+(1-tauk)*(dbart(tp)-delta)-(1+popgrowthvec(tp0))*ygrowth)+trbar0;
                labortax = labortax-taxes-beqtnew(tp);
				taultnew(tp) = labortax /(wt(tp)*ntnew(tp));
				btnew(tp) = bigb;
			elseif policy==2;
				if tp>ndebt
					labortax = bigg+btold(ndebt)*(1+(1-tauk)*(dbart(tp)-delta)-(1+popgrowthvec(tp0))*ygrowth);
                    labortax = labortax +trbar0-taxes-beqtnew(tp);
					taultnew(tp) = labortax/(wt(tp)*ntnew(tp));
					
					if tp>ndebt+1
						btnew(tp)=btnew(ndebt+1);
                    end
					
                else
					
					taultnew(tp) = taul0;
					taxes = taxes+(tauntold(tp)-tauptold(tp)) * wt(tp) *ntnew(tp);
					btnew(tp+1) = (1+(1-tauk)*(dbart(tp)-delta))*btold(tp)+bigg+trtnew(tp)-beqtnew(tp)-taxes;
					btnew(tp+1) = btnew(tp+1)/ (1+popgrowthvec(tp0)*ygrowth);					
                end
            end
        end
		
		ksharetnew(tp) = ktnew(tp) / (ktnew(tp)+btnew(tp));
		
		if case_pen==1
			if case_level_pen1==0
				pentnew(tp) = replacement_ratio*wt(tp)*meanlabtnew(tp)/sum(massvec(1:nw,tp));
			elseif case_level_pen1==1
				pentnew(tp) = pen0;
            end
		
			tauptnew(tp) = pentnew(tp) * sum(massvec(nw+1:nage,tp))/(wt(tp)*ntnew(tp));
        else			
			disp('wrong parameter for case_pen');
			pause;
        end	

     end
	
     tauntnew=taultnew+tauptnew;
      
end

