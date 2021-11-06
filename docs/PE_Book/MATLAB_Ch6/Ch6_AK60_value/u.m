
% instantaneous utility
function [y] = u(c,l)

    def_global_AK_value

    if sigma==1
        y=ln(c+psi0)+gam*ln(1-l);
    else
        y=(((c+psi0)*(1-l)^gam)^(1-sigma)-1)/(1-sigma);
    end
end