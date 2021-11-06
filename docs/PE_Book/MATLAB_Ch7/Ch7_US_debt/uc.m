function [ y ] = uc(c,l)
% marginal utility of consumption
%
    def_global_USdebt
	y=c^(-eta1)*(1-kappa*(1-eta1)*l^(1+1/varphi))^(eta1);

end

