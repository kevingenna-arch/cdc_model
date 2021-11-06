
% uc(c,l)
%
% marginal utility of consumption
function [y] = uc(c,l)
    def_global_AK_value
    y=(c+psi0)^(-sigma).*(1-l)^(gam*(1-sigma));
end
