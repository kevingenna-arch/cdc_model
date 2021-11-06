function [ c,sy,penss,taxes ] = ss_computec( x )
%
% computes c, s/y, pen and taxes in steady state
%

    def_global_USdebt
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
	bigbnew=debt;
	
	
	k=ksharess*asset;
	b=(1-ksharess)*asset;
	if nr>0
		labors=[labors; zeros(nr,1)];
    end

    % update of N, K, transfer, kshare, and pen
	nopt1=ef.*labors(1:nw);
	nnew=mass(1:nw)'*nopt1;
	
	mean_labor=nopt1'*mass(1:nw);
	mean_labor=mean_labor/sum(mass(1:nw));
	
	if case_pen==1
		penss=replacement_ratio*wbarss*mean_labor;
		taupnew=penss*sum(mass(nw+1:nage))/(wbarss*nbar);

	elseif case_pen==0
        
		penss = taupss*nbar*wbarss/sum(mass(nw+1:nage));
    else
		disp('wrong parameter for case_pen');
		pause;
    end	
	
	% consumption
	c=zeros(nage,1);	
    sy = c;
    
	for i=1:1:nw	
		c(i)=(1-taunss-taupss)*wbarss*ef(i)*labors(i)+(1+(1-tauk)*(dbarss-delta))*(k(i)+b(i))-ygrowth*(k(i+1)+b(i+1))+trbarss;
		c(i)=c(i)/(1+tauc);
        % savings rate
		sy(i)=(1-taunss-taupss)*wbarss*ef(i)*labors(i)+(1-tauk)*(dbarss-delta)*(k(i)+b(i))+trbarss-(1+tauc)*c(i);
		sy(i)=sy(i)/((1-taunss-taupss)*wbarss*ef(i)*labors(i)+(1-tauk)*(dbarss-delta)*(k(i)+b(i))+trbarss);
    end

	if nr>0
		for i=1:1:nr-1
			c(i+nw) = penss +(1+(1-tauk)*(dbarss-delta))*(k(nw+i)+b(nw+i))-ygrowth*(k(nw+i+1)+b(nw+i+1))+trbarss;
			c(i+nw) = c(i+nw)/(1+tauc);
			sy(i+nw) = penss+(1-tauk)*(dbarss-delta)*(k(i+nw)+b(i+nw))+trbarss-(1+tauc)*c(i+nw);
			sy(i+nw) = sy(i+nw)/(pen+(1-tauk)*(dbarss-delta)*(k(i+nw)+b(i+nw))+trbarss);
        end
		c(nr+nw) = penss+(1+(1-tauk)*(dbarss-delta))*(k(nw+nr)+b(nw+nr))+trbarss;
		c(nr+nw) = c(nr+nw)/(1+tauc);
        
		sy(nr+nw) = penss+(1-tauk)*(dbarss-delta)*(k(nr+nw)+b(nr+nw))+trbarss-(1+tauc)*c(nr+nw);
		sy(nr+nw) = sy(nr+nw)/(penss+(1-tauk)*(dbarss-delta)*(k(nr+nw)+b(nr+nw))+trbarss);
    end
	

    % computation of the aggregate capital stock and employment nbar 

    anew=mass'*asset(1:nage);       % aggregate wealth
	
	bigc=mass'*c;                   % total consumption
	taxes=taunss*wbarss*nbar+tauk*(dbarss-delta)*anew+tauc*bigc;
    
end

