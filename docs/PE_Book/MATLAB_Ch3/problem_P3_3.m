function Ch3_olg_dyn1()
%Code computes the dynamics of the Example economy

% %author: Burkhard Heer

%this version: June 6, 2018

clc;
close all;
dbstop if error


%% PARAMETERS
par.alpha1   = 0.36;
par.beta1    = 0.40;
par.n      = 0.1;  % populatoin growth rate
par.tr = 0;
    

%% STEADY STATE

par.kss     = (par.beta1/(1+par.beta1)*(1-par.alpha1)/(1+par.n) )^(1/(1-par.alpha1));
kss=par.kss;
yss= kss^(par.alpha1);
rss= par.alpha1*kss^(par.alpha1-1);
ntr=200;    % rows of transfer grid

% maximum transfers from old to the young when the transfer program is initiated:
% size all the savings from the old (including interest)
trmax=kss*(1+rss);

trt=linspace(0,trmax,ntr);

kt=zeros(ntr,1);       % solution capital stock
c1t=kt; c2t=kt; utilt=kt;
kguess=kss;

for i=1:ntr,
    par.tr=trt(i);
    ksol = fsolve(@(k)findk_P3_3(par,k),kguess);
    kguess=ksol;
    kt(i)=ksol;
    s=(1+par.n)*ksol; % savings
    w=wage(par,ksol);
    r=interest(par,ksol);
    c1=par.tr+w-s;
    c2=par.beta1*(1+r)*c1;
    u0=log(c1)+par.beta1*log(c2);
    c1t(i)=c1;
    c2t(i)=c2;
    utilt(i)=u0;
end;





figure
plot(trt,utilt); hold on
xlabel('transfers')
ylabel('welfare')




save problem_P3_3;

end

%% Auxiliary functions
%%
%% Wage function
function [yy] = wage(par,kk)

yy = (1-par.alpha1)*kk^(par.alpha1);

end

%% Interest Rate function
function [yy] = interest(par,kk)

yy = (par.alpha1)*kk^(par.alpha1-1);

end


function [y] = findk_P3_3(par,k)
% non-linear equation in Problem 3.3
% finds the solution to the steady-state capital stock
% see also solutions_manual (download from my homepage)
	
    alpha1=par.alpha1;
    n=par.n;
    beta1=par.beta1;
    tr=par.tr;
    w=wage(par,k);
    r=interest(par,k);
   
    y = (1+n)*k-beta1/(1+beta1)*w-tr/(1+beta1)*(1+beta1+beta1*r+n)/(1+r);
    
end