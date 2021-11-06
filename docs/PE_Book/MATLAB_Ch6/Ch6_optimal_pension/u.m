function [ y ] = u(c,h)
% computes utility
% c - consumption
% h - working hours

    def_global
	if sigma==1
        y = gam*log(c)+(1-gam)*log( 1-h);
    else
		y= ( (c)^gam * (1-h)^(1-gam))^(1-sigma) /(1-sigma);
    end
end

