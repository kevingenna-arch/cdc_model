% rvalue(x,it)
% returns the interpolated value of
% the retired agent's value function
% at age t+it with wealth x
function [y] = rvalue(x,it)
    def_global_AK_value
    if (x<assetmin)
        y=neg;
    elseif (x>assetmax)
        y=vr(na,it);
    else    % interpolation
        vr1 = vr(:,it);
        y = interp1(agrid,vr1,x,VI_method);
    end
end