function [ y1,y2,y3,y4,y5,y6 ] = getvaluess(kbar,lbar,averagehours,trbar0,taub,tauw0)
% computes the steady state


%	local knew,lnew,hoursnew,trnew,taubnew, tauwnew;
%	local wbar, rbar, penbar, bigk;
%	local wseq1, rseq1, penseq1, trseq1, tauwseq1, taubseq1, hickstrseq1;
%	local bigl, biga, bigtax, bigcontrib, bigbequests, bigpensions, xbar0, bigc;
%	local measure, c, p, prob, in1, lambda, l0, a1, wageincome, netincome, w0, grossincome;
%	local ia1, lambda1, pensionincome, bigtax0, p1,d1,p2,d2, wtaxes, totalmeasure;

    def_global
	wbar=wagerate(kbar,lbar);
	rbar=interest(kbar,lbar);
	
    % pension dependent on efficiency type?
	if case_pension==1
		% net replacement ratio
		penbar=pen_repl *(1-taub-tauw0)* wbar * averagehours*eps;
	elseif case_pension==2
		 penbar=pen_repl  *(1-taub-tauw0)* wbar * averagehours*ones(neps,1);
    end
	
    
	wseq1=ones(nage,1)*wbar;
	rseq1=ones(nage,1)*rbar;
	penseq1=ones(nage,1)*penbar';
	trseq1=ones(nage,1)*trbar0;
	tauwseq1=ones(nage,1)*tauw0;
	taubseq1=ones(nage,1)*taub;
	
	[lopt,copt,aopt,v]=getoptpolicy(wseq1,rseq1,penseq1,trseq1,tauwseq1,taubseq1);
	
    %
	%  aggregation	
	%	
	% computation of the aggregate capital stock and employment nbar 
	ga=zeros(nage,neps,neta,na);    % distribution function
	totalmeasure=0;

	% age-profiles
	agen=zeros(nage,2);		% assets workers
	cgen=zeros(nage,2);		% consumption epsilon=1,2
	lgen=zeros(t,1);		% working hours workers
	lgeneps=zeros(t,2);		% working hours, for each epsilon
	hoursnew=0;
	bigl=0;					% aggregate labor 
	biga=0;					% aggregate assets
	bigtax=0;				% aggregate taxes
	bigcontrib=0;			% aggregate contributions to pension system
	bigbequests=0;			% aggregate accidental bequests
	bigpensions=0;			% total pensions
	xbar0=0;				% average contributions to pension system
	bigc=0;					% aggregate consumption

	% distribution functions for Gini coefficients
	fa=zeros(na,1);         % distribution of wealth a
	fnetincome=zeros(nin,1);    % net income
	fgrossincome=zeros(nin,1);	% gross income
	fwage=zeros(nin,1);			% wage income
	fpension=zeros(nin,1);		%  pension income
	fconsumption=zeros(nin,1);	% consumption

	% initialization at age 1
	% equal distribution of abilities
	ga(1,1,1,1)=1/2*1/2*mass(1);
	ga(1,2,1,1)=1/2*1/2*mass(1);
	ga(1,1,2,1)=1/2*1/2*mass(1);
	ga(1,2,2,1)=1/2*1/2*mass(1);


    for it=1:1:nage-1
		for ia=1:1:na   % asset holding in period t 
			asset0=agrid(ia);
			for ieps=1:1:neps
				for ieta=1:1:neta
					measure=ga(it,ieta,ieps,ia);    % measure of household it,ieta,ieps,ia
					
					% if ismissing(measure) 
					%	pause;
                    % end
							
					c=copt(it,ieta,ieps,ia);
					bigc=bigc+c*measure;
					
					if c<=0
						fconsumption(1)=fconsumption(1)+measure;
					elseif c>=incomemax
						fconsumption(nin)=fconsumption(nin)+measure;
                    else
						in1=sum(c>incomegrid);
						lambda=(incomegrid(in1+1)-c) / (incomegrid(in1+1)-incomegrid(in1) );
						fconsumption(in1)=fconsumption(in1)+lambda*measure;
						fconsumption(in1+1)=fconsumption(in1+1)+(1-lambda)*measure;
                    end
					
					
					fa(ia)=fa(ia)+measure;	% wealth distribution
					if it<=t
						l0=lopt(it,ieta,ieps,ia);	% working hours at age it	
                    end
					
					a1=aopt(it,ieta,ieps,ia);		% optimal next-period assets
					
					cgen(it,ieps)=cgen(it,ieps)+measure*c;
					agen(it,ieps)=agen(it,ieps)+measure*asset0;
					bigtax=bigtax+taur*rbar*asset0*measure;		% interest rate tax worker
					totalmeasure=totalmeasure+measure;
						
						
					if it<=t
						lgen(it)=lgen(it)+measure*l0;
						lgeneps(it,ieps)=lgeneps(it,ieps)+l0*measure;
						bigl=bigl+l0*eta1(ieta)*eps1(ieps)*efage(it)*measure;		% effective labor supply				
						bigtax=bigtax+tauw0*l0*efage(it)*eps1(ieps)*eta1(ieta)*wbar*measure;		% wage income tax
						bigcontrib=bigcontrib+taub*l0*efage(it)*eps1(ieps)*eta1(ieta)*wbar*measure;	% pension contributions			
						wageincome=l0*efage(it)*eps1(ieps)*eta1(ieta)*wbar;
						netincome=(1-tauw0-taub)*wageincome+trbar0+(1-taur)*rbar*asset0;
						
						% wage income distribution: only employed workers
						if wageincome<=0
							fwage(1)=fwage(1)+measure;
						elseif wageincome>=incomemax
							fwage(nin)=fwage(nin)+measure;
                        else
							in1=sum(wageincome>incomegrid);
							lambda=(incomegrid(in1+1)-wageincome) / (incomegrid(in1+1)-incomegrid(in1) );
							fwage(in1)=fwage(in1)+lambda*measure;
							fwage(in1+1)=fwage(in1+1)+(1-lambda)*measure;
                        end
							
						grossincome=wageincome+trbar0+rbar*asset0;
					
                    else
						bigpensions=bigpensions+penbar(ieps)*measure;
						pensionincome=penbar(ieps);
						grossincome=penbar(ieps)+trbar0+rbar*asset0;
						netincome=penbar(ieps)+trbar0+(1-taur)*rbar*asset0;	
                    end
						
					
					biga=biga+agrid(ia)*measure;

					% income distributions
					
					if grossincome<=0
						fgrossincome(1)=fgrossincome(1)+measure;
					elseif grossincome>=incomemax
						fgrossincome(nin)=fgrossincome(nin)+measure;
                    else
						in1=sum(grossincome>incomegrid);
						lambda=(incomegrid(in1+1)-grossincome) / (incomegrid(in1+1)-incomegrid(in1) );
						fgrossincome(in1)=fgrossincome(in1)+lambda*measure;
						fgrossincome(in1+1)=fgrossincome(in1+1)+(1-lambda)*measure;
                    end

					if netincome<=0
						fnetincome(1)=fnetincome(1)+measure;
					elseif netincome>=incomemax
						fnetincome(nin)=fnetincome(nin)+measure;
                    else
						in1=sum(netincome>incomegrid);
						lambda=(incomegrid(in1+1)-netincome) / (incomegrid(in1+1)-incomegrid(in1) );
						fnetincome(in1)=fnetincome(in1)+lambda*measure;
						fnetincome(in1+1)=fnetincome(in1+1)+(1-lambda)*measure;
                    end

					if it>t
						if pensionincome<=0
							fpension(1)=fpension(1)+measure;
						elseif pensionincome>=incomemax
							fpension(nin)=fpension(nin)+measure;
                        else
							in1=sum(pensionincome>incomegrid);
							lambda=(incomegrid(in1+1)-pensionincome) / (incomegrid(in1+1)-incomegrid(in1) );
							fpension(in1)=fpension(in1)+lambda*measure;
							fpension(in1+1)=fpension(in1+1)+(1-lambda)*measure;
                        end
                    end
							

					% computation of next-periods distribution: linear interpolation 
					%  for a(ia1)<= a1 <= a(ia1+1) 
					
					if a1<=assetmin
						ga(it+1,1,ieps,1)=ga(it+1,1,ieps,1)+pi_eta(ieta,1)*sp1(it)/(1+gn)*measure;
						ga(it+1,2,ieps,1)=ga(it+1,2,ieps,1)+pi_eta(ieta,2)*sp1(it)/(1+gn)*measure;
					elseif a1>=assetmax
						ga(it+1,1,ieps,na)=ga(it+1,1,ieps,na)+pi_eta(ieta,1)*sp1(it)/(1+gn)*measure;
						ga(it+1,2,ieps,na)=ga(it+1,2,ieps,na)+pi_eta(ieta,2)*sp1(it)/(1+gn)*measure;	
                    else
						ia1=sum(agrid<a1);
						lambda1=(agrid(ia1+1)-a1) / (agrid(ia1+1)-agrid(ia1) );
						ga(it+1,1,ieps,ia1)=ga(it+1,1,ieps,ia1)	+lambda1 * pi_eta(ieta,1)*sp1(it)/(1+gn)*measure;
						ga(it+1,2,ieps,ia1)=ga(it+1,2,ieps,ia1)	+lambda1 * pi_eta(ieta,2)*sp1(it)/(1+gn)*measure;
						ga(it+1,1,ieps,ia1+1)=ga(it+1,1,ieps,ia1+1)	+ (1-lambda1) * pi_eta(ieta,1)*sp1(it)/(1+gn)*measure;
						ga(it+1,2,ieps,ia1+1)=ga(it+1,2,ieps,ia1+1)+ (1-lambda1) * pi_eta(ieta,2)*sp1(it)/(1+gn)*measure;
                    end
                end % ieta
            end % ieps
        end	% ia
    end % it


	% aggregation for the last period
	it=nage;
	for ia=1:1:na
		asset0=agrid(ia);
		for ieps= 1:1:neps
			for ieta=1:1:neta
				measure=ga(it,ieta,ieps,ia);	% measure of household it,ieta,ieps,ia,ix
				totalmeasure=totalmeasure+measure;		
				c=copt(it,ieta,ieps,ia);	
				bigc=bigc+c*measure;
					
				if c<=0
					fconsumption(1)=fconsumption(1)+measure;
				elseif c>=incomemax
					fconsumption(nin)=fconsumption(nin)+measure;
                else
					in1=sum(c>incomegrid);
					lambda=(incomegrid(in1+1)-c) / (incomegrid(in1+1)-incomegrid(in1) );
					fconsumption(in1)=fconsumption(in1)+lambda*measure;
					fconsumption(in1+1)=fconsumption(in1+1)+(1-lambda)*measure;
                end
							
				fa(ia)=fa(ia)+measure;	% wealth distribution
	
				grossincome=penbar(ieps)+trbar0+rbar*asset0;
				netincome=penbar(ieps)+trbar0+(1-taur)*rbar*asset0;	
				pensionincome=penbar(ieps);

				if grossincome<=0
					fgrossincome(1)=fgrossincome(1)+measure;
				elseif grossincome>=incomemax
					fgrossincome(nin)=fgrossincome(nin)+measure;
                else
					in1=sum(grossincome>incomegrid);
					lambda=(incomegrid(in1+1)-grossincome) / (incomegrid(in1+1)-incomegrid(in1) );
					fgrossincome(in1)=fgrossincome(in1)+lambda*measure;
					fgrossincome(in1+1)=fgrossincome(in1+1)+(1-lambda)*measure;
                end

				if netincome<=0
					fnetincome(1)=fnetincome(1)+measure;
				elseif netincome>=incomemax
					fnetincome(nin)=fnetincome(nin)+measure;
                else
					in1=sum(netincome>incomegrid);
					lambda=(incomegrid(in1+1)-netincome) / (incomegrid(in1+1)-incomegrid(in1) );
					fnetincome(in1)=fnetincome(in1)+lambda*measure;
					fnetincome(in1+1)=fnetincome(in1+1)+(1-lambda)*measure;
                end
							
							
				if pensionincome<=0
					fpension(1)=fpension(1)+measure;
				elseif pensionincome>=incomemax
					fpension(nin)=fpension(nin)+measure;
                else
					in1=sum(pensionincome>incomegrid);
					lambda=(incomegrid(in1+1)-pensionincome) / (incomegrid(in1+1)-incomegrid(in1) );
					fpension(in1)=fpension(in1)+lambda*measure;
					fpension(in1+1)=fpension(in1+1)+(1-lambda)*measure;
                end
							
				cgen(it,ieps)=cgen(it,ieps)+measure*c;
				agen(it,ieps)=agen(it,ieps)+measure*asset0;
				bigtax=bigtax+taur*rbar*asset0*measure;		% interest rate tax 75-year old
					
				bigpensions=bigpensions+penbar(ieps)*measure;
				biga=biga+asset0*measure;
					
            end	% ieta
        end	% ieps
    end % ia
	
	disp('aggregation complete'); 
	
	% computation of the Ginis
	fwage=fwage/sum(fwage);     % normalization of mass of workers to one
	fpension=fpension/sum(fpension);
	giniwage=ginid(incomegrid,fwage);
	ginigrossincome=ginid(incomegrid,fgrossincome);
	gininetincome=ginid(incomegrid,fnetincome);
	giniwealth=ginid(agrid,fa);
	ginipension=ginid(incomegrid,fpension);
	giniconsumption=ginid(incomegrid,fconsumption);
	progressivityindex=1-ginipension/giniwage;
					
	% total bequests: check if equilibrium condition is fine
	biga=sum(sum(agen));
	bigk=biga;
	bigbequests=(1-sp1)'*agen;
	bigbequests=sum(bigbequests');
	bigbequests=(1+(1-taur)*rbar)*bigbequests;
	bigy=kbar^(alpha1)*lbar^(1-alpha1);
	hoursnew=sum(lgen)/sum(mass(1:t));

	bigtax0=taur*rbar*biga+tauw0*wbar*bigl;


 biga
 bigl
 bigy
 disp('totalmeasure=1?: ');
 totalmeasure
 bigpensions
 bigtax
 bigtax0
 disp('bigtax/bigy: '); 
 bigtax0/bigy
 giniwealth
 ginigrossincome
 gininetincome
 giniwage
 ginipension
 giniconsumption
 progressivityindex
 averagehours
 rbar
 trbar0


% average value of newborn
welfare=1/4*(v(1,1,1,1)+v(1,1,2,1)+v(1,2,1,1)+v(1,2,2,1));	
welfare 
 
totalpension=bigpensions;
taubnew=bigpensions/(wbar*bigl);
if calib==1
	gbar=gy*bigy;
end
	
if endog_tr==1		% endogenous transfers
	trnew=bigtax0+bigbequests-gbar;
	tauwnew=tauw0;
else				% endogenous wage tax
	wtaxes=gbar+trbar-bigbequests-taur*rbar*biga;
	tauwnew=wtaxes/(wbar*bigl);
	trnew=trbar0;
end
	
	
	
% if _ssfinal==1; save bigbequestsfinal=bigbequests; endif;	

% return steady-state state variables K, L, l, tr, taub, tauw
y1 = bigk;
y2 = bigl;
y3 = hoursnew;
y4 = trnew;
y5 = taubnew;
y6 = tauwnew;
end

