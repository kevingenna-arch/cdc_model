function [ ktnew, bigltnew, tauwtnew, taubtnew, trbartnew, averagehourstnew ] = gettranspath()
% computes the transition dynamics 

    def_global
    ktnew=zeros(ntrans,1);
    bigltnew=zeros(ntrans,1);
    tauwtnew=zeros(ntrans,1);
    taubtnew=zeros(ntrans,1);
	trbartnew=zeros(ntrans,1);
    beqtnew=zeros(ntrans,1);
	averagehourstnew=zeros(ntrans,1);
	penbartnew=zeros(ntrans,1);
	
	wseq1=zeros(nage,1);
	rseq1=zeros(nage,1);
	penseq1=zeros(nage,2);
	trseq1=zeros(nage,1);
	tauwseq1=zeros(nage,1);
	taubseq1=zeros(nage,1);
	
	gainitial0=zeros(nage,neps,neta,na);
	
    tt=ntrans+1; 
	
	
    % computation of the optimal allocation of an agent born in period tt 
	nage0=nage-1; 
    clc;
	
    for tt=ntrans:-1:-nage0+1
        
        tt
        itrans
        krittrans
		if itrans>2				
			krittranst(3:itrans)'
        end	
		
			
        % wage, interest rate and taxes over the lifetime of the individual 
        if tt>ntrans-nage0 % agent is alive after period ntrans 
            wseq1(1:ntrans-tt+1)=waget(tt:ntrans);      % CHECK: Column vectors?
            wseq1(ntrans-tt+2:nage)=ones(nage-ntrans+tt-1,1)*wbarfinal;
            rseq1(1:ntrans-tt+1)=rbart(tt:ntrans);
            rseq1(ntrans-tt+2:nage)=ones(nage-ntrans+tt-1,1)*rbarfinal;
            penseq1(1:ntrans-tt+1,:)=penbart(tt:ntrans,:);
			if policy_change==0
				penfinal=pen_repl *(1-taubfinal-tauwfinal)* wbarfinal* averagehoursfinal;
            else
				penfinal=pen_repl_new *(1-taubfinal-tauwfinal)* wbarfinal* averagehoursfinal;
            end
			
			if case_pension==1
				% net replacement ratio				
				penseq1(ntrans-tt+2:nage,:)=ones(nage-ntrans+tt-1,1)*penfinal*eps';
			elseif case_pension==2
				% net replacement ratio
				penseq1(ntrans-tt+2:nage,:)=ones(nage-ntrans+tt-1,1)*penfinal*ones(1,neps);
            end
			
            trseq1(1:ntrans-tt+1) = trtold(tt:ntrans);
            trseq1(ntrans-tt+2:nage) = ones(nage-ntrans+tt-1,1)*trbarfinal;
            tauwseq1(1:ntrans-tt+1) = tauwtold(tt:ntrans);
            tauwseq1(ntrans-tt+2:nage) = ones(nage-ntrans+tt-1,1)*tauwfinal;
            taubseq1(1:ntrans-tt+1) = taubtold(tt:ntrans);
            taubseq1(ntrans-tt+2:nage) =ones(nage-ntrans+tt-1,1)*taubfinal;
			
			
        elseif tt<1
			
            wseq1(1:1-tt) = ones(-tt+1,1)*wbarinitial;
            wseq1(2-tt:nage) = waget(1:nage0+tt);
			
            rseq1(1:1-tt)=ones(-tt+1,1)*rbarinitial;
            rseq1(2-tt:nage)=rbart(1:nage0+tt);
			
			peninitial=pen_repl_old *(1-taubinitial-tauwinitial)* wbarinitial* averagehoursinitial;

			if case_pension==1
				% net replacement ratio				
				penseq1(1:1-tt,:) = ones(-tt+1,1)*peninitial*eps';
			elseif case_pension==2
				% net replacement ratio
				penseq1(1:1-tt,:) = ones(-tt+1,1)*peninitial*ones(1,neps);
            end
            penseq1(2-tt:nage,:) = penbart(1:nage0+tt,:);
			
            trseq1(1:1-tt) = ones(-tt+1,1)*trbarinitial;
            trseq1(2-tt:nage) = trtold(1:nage0+tt);
			
            tauwseq1(1:1-tt) = ones(-tt+1,1)*tauwinitial;
            tauwseq1(2-tt:nage) = tauwtold(1:nage0+tt);
			
            taubseq1(1:1-tt) = ones(-tt+1,1)*taubinitial;
            taubseq1(2-tt:nage) = taubtold(1:nage0+tt);
			
			
        else
			
            wseq1(1:nage) = waget(tt:tt+nage0);
            rseq1(1:nage) = rbart(tt:tt+nage0);
            trseq1(1:nage) = trtold(tt:tt+nage0);
            penseq1(1:nage,:) = penbart(tt:tt+nage0,:);
            tauwseq1(1:nage) = tauwtold(tt:tt+nage0);
            taubseq1(1:nage) = taubtold(tt:tt+nage0);
			
			
        end
		
		
		
%	 	survival probabilities of the household born in period tt
        for isp = 1:1:nage
			sp1(isp) = surviveprob(isp,tt+isp-1);
        end	
		
        % prepare all input vectors to be column vectors
        if size(wseq1,2)>1
            wseq1 = wseq1';
        end
        
        if size(rseq1,2)>1
            rseq1 = rseq1';
        end
        
        if size(penseq1,2)>2
            penseq1 = penseq1';
        end
        
        if size(trseq1,2)>1
            trseq1 = trseq1';
        end
        
        if size(tauwseq1,2)>1
            tauwseq1 = tauwseq1';
        end
        
        if size(taubseq1,2)>1
            taubseq1 = taubseq1';
        end
        
		% computation of the optimal policy function for the generation born in period tt
		[ lopt1,copt1,aopt1,v ] = getoptpolicy(wseq1,rseq1,penseq1,trseq1,tauwseq1,taubseq1);


		lifetimeutiltemp = 1/4*(v(1,1,1,1)+v(1,1,2,1)+v(1,2,1,1)+v(1,2,2,1)); % Check: Scalar?
        
		if tt>=1
			lifetimeutilt(tt+nage) = lifetimeutiltemp;
			if policy_change==1         
				deltatnew(tt+nage)= (lifetimeutilt(tt+nage)/lifetimeutiltbenchmark(tt+nage))^(-1/gam)-1;	% consumption equivalent change
				deltat(tt+nage,itrans) = deltatnew(tt+nage);
            end
        else
			lifetimeutilt(nage+tt) = lifetimeutiltemp;
			if policy_change==1	% household is born in -tt and has age 1-tt in period 1
				deltatnew(nage+tt) = (lifetimeutiltemp/lifetimeutiltbenchmark(nage+tt))^(-1/gam)-1;	% consumption equivalent change
				deltat(nage+tt,itrans) = deltatnew(nage+tt);
            end		
        end

%		// 1. check if optimal policies in last period of transition
%		// are equal to those in the final steady state
%		//
%		//if tt==ntrans;
%		//save aopt1, lopt1, copt1;
%		//save v, ve, vu;
%		//"optimal policy functions saved";
%		//wait;
%		//endif;
%
%		// 2. check if optimal policy functions in final steady state imply
%
%		// run the following with the optimal policy functions in the final steady state
%		
%		//lopt1=loptfinal;
%		//copt1=coptfinal;
%		//aopt1=aoptfinal;


%
%       Aggregation 
%

		ga0=zeros(nage,neps,neta,na);
		
		% initialization at age 1
		% equal distribution of abilities
		
		if tt>0
			ga0(1,1,1,1) = 1/2*1/2*massvec(1,tt);		% mass of the 1-year-old in period tt
			ga0(1,2,1,1) = 1/2*1/2*massvec(1,tt);
			ga0(1,1,2,1) = 1/2*1/2*massvec(1,tt);
			ga0(1,2,2,1) = 1/2*1/2*massvec(1,tt);
			
        else
			ga0(1,1,1,1) = 1/2*1/2*massvec(1,1);		% mass of the 1-year-old in period tt
			ga0(1,2,1,1) = 1/2*1/2*massvec(1,1);
			ga0(1,1,2,1) = 1/2*1/2*massvec(1,1);
			ga0(1,2,2,1) = 1/2*1/2*massvec(1,1);
        end
		
		
		for tp = 1:1:nage-1     % age of the generation born in period tt
						
			for iap = 1:1:na    % asset holding in period t */
				asset0p=agrid(iap);
				for iepsp = 1:1:neps 
					for ietap = 1:1:neta
						measure=ga0(tp,ietap,iepsp,iap);	% measure of household it,ieta,ieps,ia
					
						if tp+tt-1<=1		%	in the first transition period at age tp, transition period tt+tp-1=1
											%  the distribution of the assets is equal to the initial distribution
											% also sets the mass of the population equal to one in transition period 0
							measure=gainitial(tp,ietap,iepsp,iap);
                        end
							
						if tp<=t
							l0=lopt1(tp,ietap,iepsp,iap);	% working hours at age it
                        end
					
					
						a1=aopt1(tp,ietap,iepsp,iap);		% optimal next-period assets

						if (tt+tp>0) && (tt+tp<=ntrans)		% bequest of individual who diese at age tp 
															% period is during the transition?			
							% bequests of cohort in tt+tp-1 
							beqtnew(tt+tp) = beqtnew(tt+tp)+(1-sp1(tp))*(1+(1-taur)*rbart(tt+tp))*a1*measure;
                        end
							
						if (tt+tp-1>0) && (tt+tp-1<=ntrans)		% individual at age tp is alive during the transition
							ktnew(tt+tp-1) = ktnew(tt+tp-1)+asset0p*measure;
							massvec1(tt+tp-1) = massvec1(tt+tp-1)+measure;
							massvec0(tp,tt+tp-1) = massvec0(tp,tt+tp-1)+measure;
							if tp<=t
								averagehourstnew(tt+tp-1) = averagehourstnew(tt+tp-1)+measure*l0;
								bigltnew(tt+tp-1) = bigltnew(tt+tp-1)+l0*eta1(ietap)*eps1(iepsp)*efage(tp)*measure;
							elseif tp>t
								penbartnew(tt+tp-1) = penbartnew(tt+tp-1)+penbart(tt+tp-1,iepsp)*measure;
                            end
				
				
                        end		% cohort tp in transition period tt=1,..,ntrans?	

						

					% computation of next-periods distribution: bilinear interpolation 
					%  for a[ia1]<= a1 <= a[ia1+1] 
					%  and x[ix1]<= averagex1 <= x[ix1+1]
					%
					% ga(it,ieta,ieps,ia)
					% check if a1<=0
					
					
					
						if a1<=assetmin
							ga0(tp+1,1,iepsp,1) = ga0(tp+1,1,iepsp,1)+pi_eta(ietap,1)*sp1(tp)*measure;
							ga0(tp+1,2,iepsp,1) = ga0(tp+1,2,iepsp,1)+pi_eta(ietap,2)*sp1(tp)*measure;
						
						elseif a1>=assetmax
						
							ga0(tp+1,1,iepsp,na) = ga0(tp+1,1,iepsp,na)+pi_eta(ietap,1)*sp1(tp)*measure;
							ga0(tp+1,2,iepsp,na) = ga0(tp+1,2,iepsp,na)+pi_eta(ietap,2)*sp1(tp)*measure;
						
                        else
							ia1=sum(agrid<a1);
							lambda1=(agrid(ia1+1)-a1) / (agrid(ia1+1)-agrid(ia1) );
						
							ga0(tp+1,1,iepsp,ia1) = ga0(tp+1,1,iepsp,ia1)	+lambda1 * pi_eta(ietap,1)*sp1(tp)*measure;
							ga0(tp+1,2,iepsp,ia1) = ga0(tp+1,2,iepsp,ia1)	+lambda1 * pi_eta(ietap,2)*sp1(tp)*measure;
	
							ga0(tp+1,1,iepsp,ia1+1) = ga0(tp+1,1,iepsp,ia1+1)	+ (1-lambda1) * pi_eta(ietap,1)*sp1(tp)*measure;
							ga0(tp+1,2,iepsp,ia1+1) = ga0(tp+1,2,iepsp,ia1+1)+ (1-lambda1) * pi_eta(ietap,2)*sp1(tp)*measure;
							
                        end
                    end	% ietap
                end	% iepsp
            end	% iap
        end % tp


	% aggregation for the last period
	tp=nage;
	for iap=1:1:na
		asset0p=agrid(iap);
		for iepsp=1:1:neps
			for ietap = 1:1:neta
				
				measure=ga0(tp,ietap,iepsp,iap);	% measure of household tp,ietap,iepsp,iap
					
				if tp+tt-1==1		%	in the first transition period at age tp, transition period tt+tp-1=1
									%  the distribution of the assets is equal to the initial distribution
									% also sets the mass of the population equal to one in transition period 0
					measure=gainitial(tp,ietap,iepsp,iap);				
                end
						
				if (tt+tp-1>0) && (tt+tp-1<=ntrans)		% individual at age tp is alive during the transition
					massvec1(tt+tp-1) = massvec1(tt+tp-1)+measure;
					massvec0(tp,tt+tp-1) = massvec0(tp,tt+tp-1)+measure;
					ktnew(tt+tp-1) = ktnew(tt+tp-1)+asset0p*measure;
					penbartnew(tt+tp-1)=penbartnew(tt+tp-1)+penbart(tt+tp-1,iepsp)*measure;
                end
					
            end	% ietap
        end	% iepsp
     end % iap
	

    end     % tt - transition periods
	

	ktnew=ktnew./massyear;
	bigltnew=bigltnew./massyear;
	beqtnew=beqtnew./massyear;
	penbartnew=penbartnew./massyear;
	
	averagehourstnew=averagehourstnew./(sum(massvec(1:t,:)));	% correct for the share of the labor force in population	

	if endog_tr==0      % tauw adjust to balance gov budget;
		trbartnew=trtold;
		taxest=ones(ntrans,1)*gbar;
		taxest=taxest+trbartnew-taur*rbart.*ktnew-beqtnew;
		tauwtnew=taxest./(waget.*bigltnew);
		taubtnew=penbartnew./(waget.*bigltnew);
	
	elseif endog_tr==1  % tr adjust to balance gov budget
		tauwtnew=tauwtold;
		trbartnew=tauw*waget.*bigltnew+ taur*rbart.*ktnew+beqtnew-ones(ntrans,1)*gbar;
		taubtnew=penbartnew./(waget.*bigltnew);	
		
    end
	
end

