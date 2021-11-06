function [ y ] = gnproc(periodx)
% computes the population growth rate in period periodx;
    def_global
	if periodx>rowsspall
		periodx=rowsspall;		% all survival probs are equal to those in the year 2095
	elseif periodx<1;
		periodx=1;
    end
	y=gnall(periodx);
end	
		