function [ y ] = ss_aggregate(x)
% steady state with endogenous labor supply 

    def_global_USdebt
	y=x;
	asset=[0; x(1:nage-1); 0];		% assets of the cohorts
	labors=x(nage:nage+nw-1);       % labor supply
	kbar=x(nage+nw);
	nbar=x(nage+nw+1);
	debt=x(nage+nw+2);
	if case_pen==1
		taupss=x(nage+nw+3);
    else
		disp('wrong parameter case_pen');
		pause;
    end
	trbarss=x(nage+nw+4);
	taunss=taulbar-taupss;
	
    wbarss=wagerate(kbar,nbar); 
	dbarss=interest(kbar,nbar);
	ybar=production(kbar,nbar);
	anew=mass'*asset(1:nage);       % aggregate wealth
	knew=anew-debt;
	ksharess=kbar/anew;
	bigbnew=debtoutputratio*production(kbar,nbar);
	
	
	k=ksharess*asset;
	b=(1-ksharess)*asset;
	if nr>0
		labors=[labors; zeros(nr,1)];
    end

    % update of N, K, transfer, kshare, and pen
	nopt1=ef.*labors(1:nw);
	nnew=mass(1:nw)'*nopt1;
	y(nage+nw+3)=nnew-nbar;
	
	
	mean_labor=nopt1'*mass(1:nw);
	mean_labor=mean_labor/sum(mass(1:nw));
	
	if case_pen==1
		penss=replacement_ratio*wbarss*mean_labor;
		taupnew=penss*sum(mass(nw+1:nage))/(wbarss*nbar);

		y(nage+nw+4) = taupnew-taupss;

	elseif case_pen==0
        
		penss = taupss*nbar*wbarss/sum(mass(nw+1:nage));
    else
		disp('wrong parameter for case_pen');
		pause;
    end	
	
	% consumption
	c=zeros(nage,1);	
	
	for i=1:1:nw	
		c(i)=(1-taunss-taupss)*wbarss*ef(i)*labors(i)+(1+(1-tauk)*(dbarss-delta))*(k(i)+b(i))-ygrowth*(k(i+1)+b(i+1))+trbarss;
		c(i)=c(i)/(1+tauc);
    end

	if nr>0
		for i=1:1:nr-1
			c(i+nw) = penss +(1+(1-tauk)*(dbarss-delta))*(k(nw+i)+b(nw+i))-ygrowth*(k(nw+i+1)+b(nw+i+1))+trbarss;
			c(i+nw) = c(i+nw)/(1+tauc);
        end
		c(nr+nw) = penss+(1+(1-tauk)*(dbarss-delta))*(k(nw+nr)+b(nw+nr))+trbarss;
		c(nr+nw) = c(nr+nw)/(1+tauc);
    end
	
	% intertemporal first-order conditions
	for i=1:1:nage-1
		y(i) = ygrowth^(eta1)*uc(c(i),labors(i))/uc(c(i+1),labors(i+1))-sp(i)*beta1*(1+(1-tauk)*(dbarss-delta));
    end
	
	% optimal labor supply
	for i=1:1:nw
		y(i+nage-1) = (1-taunss-taupss)/(1+tauc)*wbarss*ef(i)*(1-kappa*(1-eta1)*labors(i)^(1+1/varphi));
		y(i+nage-1) = y(i+nage-1) -kappa*eta1*(1+1/varphi)*c(i)*labors(i)^(1/varphi);
    end



    % computation of the aggregate capital stock and employment nbar 

    anew=mass'*asset(1:nage);       % aggregate wealth
	knew=ksharess*anew;
	y(nage+nw)=kbar-knew;
	
	bigc=mass'*c;                   % total consumption
	taxes=taunss*wbarss*nbar+tauk*(dbarss-delta)*anew+tauc*bigc;
	temp=(1-sp(1:nage-1)).*mass(1:nage-1);
	bequests=temp'*asset(2:nage);
	bequests=bequests*(1+(1-tauk)*(dbarss-delta)); 
	ybarss=production(kbar,nbar);
	transfernew=taxes+bequests+debt*((1+popgrowth)*ygrowth-(1+dbarss-delta) )-bigg;
	y(nage+nw+1) = trbarss-transfernew;
%
%	bigb=debtratio*ybar;
%
	ksharenew=kbar/(kbar+bigbnew);
	y(nage+nw+2)=debt-bigbnew;

end

