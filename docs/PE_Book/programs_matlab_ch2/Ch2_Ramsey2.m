function Ch1_Ramsey2()
%Code adapted from Gauss code by Heer, S. Duineveld, 10/11/2015

clc;
close all;
dbstop if error


%% PARAMETERS
par.delta   = 0.08;
par.alpha   = 0.36;
par.nn      = 0;
par.beta    = 0.96;
par.sigma   = 2.0;

par.lss     = 0.3; %steady state labour supply

par.TT      = 40;	% number of transition periods
par.Tz_sta   = 10;  % start of shock
par.Tz_end   = 12;  % end of shock

%% STEADY STATE 
par.zss     = 1;
par.kss     = (par.alpha/(1/par.beta-1+par.delta))^(1/(1-par.alpha))*par.lss; 
par.yss     = production(par,par.zss,par.kss,par.lss);
par.css     = par.yss-(par.nn+par.delta)*par.kss;
% calibration of gamma
gam0        = (1-par.lss)*(1-par.alpha)*par.kss^(par.alpha)*par.lss^(-par.alpha)/par.css;
par.gam     = 1/(1+gam0);
clear gam0;


%% Solve for expected shock and unexpectged shock (1.35)
[irf.kt_exp,irf.lt_exp] = solve_1_35_exp(par);

[irf.kt_une,irf.lt_une] = solve_1_35_une(par);


%% Generate IRF consumption & production
irf.zt                          = par.zss*ones(1,par.TT+2);
irf.zt(1,par.Tz_sta+1:par.Tz_end+1) = 1.1*par.zss; %column = t+1

irf.ct_exp = consumption(par,irf.zt,irf.kt_exp,irf.lt_exp,[irf.kt_exp(1,2:end),par.kss]);
irf.yt_exp = production(par,irf.zt,irf.kt_exp,irf.lt_exp);
    
irf.ct_une = consumption(par,irf.zt,irf.kt_une,irf.lt_une,[irf.kt_une(1,2:end),par.kss]);
irf.yt_une = production(par,irf.zt,irf.kt_une,irf.lt_une);


%% Plot varaibles
figure
subplot(2,2,1)
plot([0:par.TT+1],irf.yt_exp)
hold all
plot([0:par.TT+1],irf.yt_une,'--');

xlabel('period')
title('production');
legend('expected shock','unexpected shock');
axis([0 par.TT+1 -inf inf]);
axis 'auto y';


subplot(2,2,2)
plot([0:par.TT+1],irf.ct_exp)
hold all
plot([0:par.TT+1],irf.ct_une,'--');

xlabel('period')
title('consumption');
legend('expected shock','unexpected shock');
axis([0 par.TT+1 -inf inf]);
axis 'auto y';


subplot(2,2,3)
plot([0:par.TT+1],irf.kt_exp)
hold all
plot([0:par.TT+1],irf.kt_une,'--');

xlabel('period')
title('capital stock');
legend('expected shock','unexpected shock');
axis([0 par.TT+1 -inf inf]);
axis 'auto y';


subplot(2,2,4)
plot([0:par.TT+1],irf.lt_exp)
hold all
plot([0:par.TT+1],irf.lt_une,'--');

xlabel('period')
title('labor');
legend('expected shock','unexpected shock');
axis([0 par.TT+1 -inf inf]);
axis 'auto y';


%% SAVE DATA

save Ch1_ramsey2;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solve equation 1.35
function [kt,lt] = solve_1_35_exp(par)

options     = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

zt                          = par.zss*ones(1,par.TT);
zt(1,par.Tz_sta:par.Tz_end) = 1.1*par.zss;

f_res_135_exp   = @(XX)res_1_35_exp(par,zt,XX(1,1:par.TT),XX(1,par.TT+1:par.TT*2));
XX0             = [repmat(par.kss,1,par.TT),repmat(par.lss,1,par.TT)];

XX              = fsolve(f_res_135_exp,XX0,options);

kt = [par.kss,XX(1,1:par.TT),par.kss]; %period 0 and period par.TT + 1 added
lt = [par.lss,XX(1,par.TT+1:2*par.TT),par.lss]; %period 0 and period par.TT + 1 added

end


%% Residual function 1.35
function [res] = res_1_35_exp(par,zt,kt,lt)

res = NaN(2*par.TT,1);
for it = 1:par.TT;
    if it == 1;
        k_p = par.kss; % k_t-1
        z_p = par.zss;
        l_p = par.lss;
    else
        k_p = kt(1,it-1); % k_t-1
        z_p = zt(1,it-1);
        l_p = lt(1,it-1);
    end
    if it == par.TT;
        k_n = par.kss; % k_t+1
    else
        k_n = kt(1,it+1);% k_t+1
    end 
    cc          = consumption(par,zt(1,it),kt(1,it),lt(1,it),k_n); 
    c_p         = consumption(par,z_p,k_p,l_p,kt(1,it));
    [NR]        = net_return(par,zt(1,it),kt(1,it),lt(1,it));
    res(it,1)   = par.beta*NR - (cc/c_p)^par.sigma;
    
    res(par.TT+it,1) = (1-par.gam)/par.gam * cc/(1-lt(1,it)) - ...
                        (1-par.alpha)* zt(1,it)*kt(1,it)^par.alpha * lt(1,it)^-par.alpha;    
end

end


%% Solve 1.35 for unexpected shock
function [kt,lt] = solve_1_35_une(par)

T_lng               = par.TT - par.Tz_sta;

options = optimset('Display','off','TolFun',1e-12,'TolX',1e-12);

zt                          = par.zss*ones(1,par.TT+1);
zt(1,par.Tz_sta:par.Tz_end) = 1.1*par.zss;

f_res_135_une   = @(XX)res_1_35_une(par,zt,XX(1,1:T_lng),XX(1,T_lng+1:2*T_lng+1));
XX0             = [repmat(par.kss,1,T_lng),repmat(par.lss,1,T_lng+1)];

XX              = fsolve(f_res_135_une,XX0,options);

kt = [repmat(par.kss,1,par.Tz_sta+1),XX(1,1:T_lng),par.kss];%add t=0:t=Tz_sta, t=TT+1 
lt = [repmat(par.lss,1,par.Tz_sta),XX(1,T_lng+1:2*T_lng+1),par.lss];%add t=0:t=Tz_sta-1, t=TT+1  

end

%% Residual function 1.35 for unexpected shock
function [res] = res_1_35_une(par,zt_f,kt,lt)

if size(zt_f,2) ~= par.TT+1
    error('length of vector zt_f and time periods not consistent')
end

T_lng   = size(kt,2);         % length of period under consideration
T_st    = par.TT - T_lng + 1;  % starting period, for capital;

kt_f    = [repmat(par.kss,1,par.TT-T_lng),kt,par.kss];  %kt over full length
lt_f    = [repmat(par.lss,1,par.TT-T_lng-1),lt,par.kss]; %lt over full length
if size(kt_f,2) ~= par.TT+1 ||  size(lt_f,2) ~= par.TT+1;
    error('length of vectors kt_f or lt_f and time periods not consistent');
end

res = NaN(2*T_lng+1,1);

w_st                = (1-par.alpha)* zt_f(1,T_st-1)*kt_f(1,T_st-1)^par.alpha * lt_f(1,T_st-1)^-par.alpha; 
c_st                = consumption(par,zt_f(1,T_st-1),kt_f(1,T_st-1),lt_f(1,T_st-1),kt_f(1,T_st));
res(T_lng+1,1)      = w_st*(1-lt_f(1,T_st-1))*par.gam-(1-par.gam)*c_st;

cnt = 1;
for it = T_st:par.TT;
    
    cc                  = consumption(par,zt_f(1,it),kt_f(1,it),lt_f(1,it),kt_f(1,it+1));
    c_p                 = consumption(par,zt_f(1,it-1),kt_f(1,it-1),lt_f(1,it-1),kt_f(1,it));
    NR                  = net_return(par,zt_f(1,it),kt_f(1,it),lt_f(1,it));
    
    res(cnt,1)          = par.beta*NR - (cc/c_p)^par.sigma;
    
    res(T_lng+cnt+1,1)  = (1-par.gam)/par.gam * cc/(1-lt_f(1,it)) - ...
                        (1-par.alpha)* zt_f(1,it)*kt_f(1,it)^par.alpha * lt_f(1,it)^-par.alpha; 
                    
    cnt                 = cnt + 1;    
end


end


%% Net return
function [NR] = net_return(par,zz,kk,ll)

NR = par.alpha * zz .* kk.^(par.alpha-1) .* ll.^(1 - par.alpha) + 1 - par.delta;

end


%% Consumption function
function [cc] = consumption(par,zz,kk,ll,k_n)

cc = production(par,zz,kk,ll) + (1-par.delta)*kk - (1+par.nn)*k_n;

end


%% PRoduction function
function [yy] = production(par,zz,kk,ll)

yy = zz .* kk.^par.alpha .* ll.^(1-par.alpha);

end
