% value1(x)
% value function of the retired
% with present asset k0 and
% next-period asset x
function [y] = value1(x)
    def_global_AK_value
    c=(1+r)*k0+pen-x;
    if (c<0)
        y=neg;
    else
        y=u(c,0)+beta0*rvalue(x,period+1);
    end
end
