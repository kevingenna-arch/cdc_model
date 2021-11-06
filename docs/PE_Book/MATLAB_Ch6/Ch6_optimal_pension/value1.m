function [ y ] = value1(x)
% computes the RHS of the Bellman equation for the retired
    def_global
	c=(1+(1-taur)*rseq(it+t))*agrid(ia)+penseq(it+t,ieps)+trseq(it+t)-x*(1+garate);
	if c<=0 
        y=neg; 
    else
        y=u(c,0)+beta1*sp1(t+it)*(1+garate)^(gam*(1-sigma))*rvalue(x);
    end
end

