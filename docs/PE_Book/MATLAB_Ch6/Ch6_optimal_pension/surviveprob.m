function [y] = surviveprob(agex,periodx)
% computes the survival probability at age agex in period periodx;
    def_global
	if periodx>rowsspall
		periodx=rowsspall;		% all survival probs are equal to those in the year 2095
	elseif periodx<1
		periodx=1;
    end
	y=spall(agex,periodx);
end
