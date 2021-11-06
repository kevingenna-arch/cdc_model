% un(c,l)
%
% marginal utility of leisure
function [y] = un(c,l)
    def_global_AK_value
    y=gam*(c+psi0)^(1-sigma).*(1-l)^(gam*(1-sigma)-1);
end