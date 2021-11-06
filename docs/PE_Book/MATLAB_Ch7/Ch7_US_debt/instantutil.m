function [ y ] = instantutil(c,n)
%instantaneous utility

    def_global_USdebt
	y=c^(1-eta1)*(1-kappa*(1-eta1)*n^(1+1/varphi))^eta1;
	y=y/(1-eta1);
end

