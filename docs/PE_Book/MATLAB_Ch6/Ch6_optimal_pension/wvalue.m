function [ y ] =  wvalue(a1)
% computes the value if the agent
% is a worker in the current period
%
% Input: 	a1 -- next period asset
% 		  
% Output: right-hand-side of the bellman equation
		
    def_global
	w0=wseq(it)*efage(it)*eps1(ieps)*eta1(ieta);
	
% computation of labor supply 
	if case_endogenouslabor==1
		% first-order conditions with respect to c and l	
		c=gam*( (1+(1-taur)*rseq(it))*asset0+trseq(it)+(1-tauwseq(it)-taubseq(it))*w0-a1*(1+garate) );
        if c>0; 
            l0=1-c/ ( (1-tauwseq(it)-taubseq(it))*w0 ) * (1-gam)/gam;
            if l0<0
                l0=0;
                c= (1+(1-taur)*rseq(it))*asset0+trseq(it)-a1*(1+garate);
            end
        else
            l0=0;
        end
    else
		l0=laborexogenous;	
		c=(1+(1-taur)*rseq(it))*asset0+trseq(it)+(1-tauwseq(it)-taubseq(it))*w0*l0-a1*(1+garate);
    end

    
	if c<=0
        y=neg;
    else
%        x1 = interp1(agrid,v11e,a1,VI_method);
%        x2 = interp1(agrid,v12e,a1,VI_method);
% interpolation too slow?

        if a1==assetmin 
            x1= v11e(1); 
            x2= v12e(1); 
  
        elseif a1>=assetmax  
            x1=v11e(na);
            x2=v12e(na);
        else
            ia0=sum(a1>agrid);
            lambda=(agrid(ia0+1)-a1) / (agrid(ia0+1)-agrid(ia0) );

            if ia0==na
                x1=v11e(ia0);
                x2=v12e(ia0);
            else
                x1=lambda*v11e(ia0)+(1-lambda)*v11e(ia0+1);
                x2=lambda*v12e(ia0)+(1-lambda)*v12e(ia0+1);
            end
        end

		expected_value=pi_eta(ieta,1)*x1+pi_eta(ieta,2)*x2;
		y=u(c,l0);	
		y=y+beta1*sp1(it)*(1+garate)^(gam*(1-sigma))*expected_value;
    end     % c>0?
end

