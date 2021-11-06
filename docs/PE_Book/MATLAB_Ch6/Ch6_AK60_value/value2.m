% value2(x): 
%
% returns the Bellman equation
%
% value function of the worker
% at age period, wealth level k0
% and next-period wealth x

function [y] = value2(x)           
    def_global_AK_value
    k1 = x;         % global variable, to be used in rf()
    
    n= 1/(1+gam)*(1-gam/((1-tau)*w)*(psi0+(1+r)*k0-k1));
    if n<0
        n=0;
    elseif n>1
        n=1;
    end
            
    c=(1-tau)*w*n+(1+r)*k0-k1;
    
    if c<0
        y=neg;
    else
        if period==t    % next period at age 40 the household retires
                        % value function of the 1-year retired agent,
                        % rvalue(:,1)
            y=u(c,n)+beta0*rvalue(x,1);
        else        % worker next period
            y=u(c,n)+beta0*wvalue(x,period+1);
        end
    end
end