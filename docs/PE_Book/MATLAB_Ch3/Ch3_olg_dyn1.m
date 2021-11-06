function Ch3_olg_dyn1()
%Code computes the dynamics of the Example economy

% %author: Burkhard Heer

%this version: June 6, 2018

clc;
close all;
dbstop if error


%% PARAMETERS
par.alpha   = 0.36;
par.beta    = 0.40;
par.TT      = 20;	% number of transition periods
par.nn      = 0.1;  % populatoin growth rate

%% STEADY STATE

par.kss     = (par.beta/(1+par.beta)*(1-par.alpha)/(1+par.nn) )^(1/(1-par.alpha));
par.yss     = production(par,par.kss)


%% Check stability of system in k_t+1,  k_t , k_t-1

kt=zeros(1,par.TT+1);

kt(1)=par.kss/3;


for j=1:par.TT;

    kt(j+1)= par.beta/(1+par.beta) * (1-par.alpha)/(1+par.nn) * kt(j)^(par.alpha);

    kt(j+1) = eq_3_21(par,kt(j));
end    
  

figure
plot([0:par.TT],kt); hold on
xlabel('Period t')
ylabel('Capital stock')


x=par.kss;

% numerical computation of the Jacobian: central differences

    dx= max(1e-8,abs(x)/1e8);
        x1 = x + dx;    
    y1  = eq_3_21(par,x1);
   
    x2 = x - dx;    
    y2  = eq_3_21(par,x2);
    
    Jac_3_21  = (y1-y2)./(x1-x2) 


save Ch3_olg_dyn1;

end

%% Auxiliary functions
%%
%% Production function
function [yy] = production(par,kk)

yy = kk.^par.alpha;

end


%% Function 3.21 in k_t 
function [f] = eq_3_21(par,kx)

f = par.beta/(1+par.beta) * (1-par.alpha)/(1+par.nn) * kx^(par.alpha);


end
