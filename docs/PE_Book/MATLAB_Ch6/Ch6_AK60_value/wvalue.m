% wvalue(x,it)
% returns the interpolated value of
% the worker's value function
% at age it with wealth x
function [y] = wvalue(x,it)
    def_global_AK_value
    if (x<assetmin)
        y=neg;
    elseif (x>assetmax)
        y=vw(na,it);
    else    % interpolation
        vw1 = vw(:,it);
        y = interp1(agrid,vw1,x,VI_method);
    end
end

