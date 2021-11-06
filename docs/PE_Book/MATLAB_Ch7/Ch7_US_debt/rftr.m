function [ y ] = rftr(x)
% first-order condtions of the household for wages wseq, interest rates rseq, pensions penseq */
%
%	local y,asset,labors,c,i,i0,l0,l1;

    def_global_USdebt
    def_global_USdebt_transition
    
    y=x;
	asset=[0; x(1:nage-1); 0];		% assets of the cohorts
	labors=x(nage:nage+nw-1);
	
	% consumption
	c=zeros(nage,1);	
	
	for i=1:1:nw	
		c(i)= (1-taulseq(i)-taupseq(i))*wseq(i)*ef(i)*labors(i);
        c(i) = c(i) +(1+(1-tauk)*(dseq(i)-delta))*asset(i)-ygrowth*asset(i+1)+trseq(i);
		c(i) = c(i)/(1+tauc);
    end

	if nr>0;
		for i=1:1:nr-1
			c(i+nw) = penseq(i+nw)+(1+(1-tauk)*(dseq(i+nw)-delta))*asset(nw+i)-ygrowth*asset(nw+i+1)+trseq(nw+i);
			c(i+nw) = c(i+nw)/(1+tauc);
        end
		c(nr+nw) = penseq(nr+nw)+( 1 + (1-tauk)*(dseq(nr+nw)-delta))*asset(nw+nr)+trseq(nw+nr);
		c(nr+nw) = c(nr+nw)/(1+tauc);
    end
	
	% intertemporal first-order conditions
	for i=1:1:nage-1
		% choice of survival probability
		% if life period is prior to 2010: sp from 2010
		% if life period is after the final year: sp from yearfinal
		i0=min(tt+i,nsp);
		i0=max(i0,1);
		if i<=nw
			l0=labors(i);
        else
			l0=0;
        end
	
		if i+1<=nw
			l1=labors(i+1);
        else
			l1=0;
        end
        
		y(i) = ygrowth^(eta1)*uc(c(i),l0)/uc(c(i+1),l1) - spvec(i,i0) * beta1*(1+(1-tauk)*(dseq(i)-delta));
		if i+tt<2       % asset demand pre-determined in period t=0,-1,-2r
			y(i) = asset(i+1) - asset0(i+1);
        end
		
    end
	
	% optimal labor supply
	for i=1:1:nw
		y(i+nage-1) = (1-taulseq(i)-taupseq(i)) / (1+tauc) * wseq(i)*ef(i)*(1-kappa*(1-eta1)*labors(i)^(1+1/varphi));
		y(i+nage-1) = y(i+nage-1) - kappa*eta1*(1+1/varphi)*c(i)*labors(i)^(1/varphi);
		if i+tt<2      % labor supply pre-determined in period t=0,-1,-2
			y(i+nage-1) = labors(i)-labors0(i);
        end
    end

end

