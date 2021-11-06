function Ch1_Ramsey1()
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
par.TT      = 40;	% number of transition periods
par.Tz_sta   = 10;  % start of shock
par.Tz_end   = 12;  % end of shock


%% STEADY STATE

par.kss     = (par.alpha/(1/par.beta-1+par.delta))^(1/(1-par.alpha));
par.yss     = production(par,1,par.kss);
par.css     = par.yss-(par.nn+par.delta)*par.kss;
par.zss     = 1;


%% Check stability of system in k_t+1,  k_t , k_t-1

Jac_115 = jac_cd(@(X)eq_1_15(par,X(1,1),X(1,2)),[par.kss,par.kss])
EG_115  = eig(Jac_115)

Jac_117 = jac_cd(@(X)eq_1_17(par,X(1,1),X(1,2)),[par.kss,par.css])
EG_117  = eig(Jac_117)

Jac_ana = jac_117_ana(par,par.kss,par.css)
EG_ana  = eig(Jac_ana)

[SS,TT] = schur(Jac_ana)


%% Solve with expected shock
irf.kt_exp = solve_1_14_exp(par);

%% Solve with unexpected shock
irf.kt_une = solve_1_14_une(par);

irf.zt                              = par.zss*ones(1,par.TT+2);
irf.zt(1,par.Tz_sta+1:par.Tz_end+1) = 1.1*par.zss;


%% Plot capital stock

figure
plot([0:par.TT+1],irf.kt_exp);
hold all
plot([0:par.TT+1],irf.kt_une,'--');

xlabel('period')
ylabel('capital stock');
legend('expected shock','unexpected shock');
axis([0 par.TT+1 -inf inf]);
axis 'auto y';


%% Plot production

irf.y_exp = production(par,irf.zt,irf.kt_exp);
irf.y_une = production(par,irf.zt,irf.kt_une);

figure
plot([0:par.TT+1],irf.y_exp);
hold all
plot([0:par.TT+1],irf.y_une,'--');

xlabel('period')
ylabel('production');
legend('expected shock','unexpected shock');
axis([0 par.TT+1 -inf inf]);
axis 'auto y';


%% Plot consumption (TIMING: time=column-1)

irf.c_exp = consumption(par,irf.zt,irf.kt_exp,[irf.kt_exp(2:end),par.kss]);
irf.c_une = consumption(par,irf.zt,irf.kt_une,[irf.kt_une(2:end),par.kss]);

figure
plot([0:par.TT+1],irf.c_exp);
hold all
plot([0:par.TT+1],irf.c_une,'--');

xlabel('period')
ylabel('consumption');
legend('expected shock','unexpected shock');
axis([0 par.TT+1 -inf inf]);
axis 'auto y';


%% CHECK LAW OF MOTION CAPITAL

irf.kt_exp_chck(1,1) = par.kss;
irf.kt_une_chck(1,1) = par.kss;
for it=2:par.TT+2;
    irf.kt_exp_chck(1,it) = irf.y_exp(1,it-1) - irf.c_exp(1,it-1) + (1-par.delta)*irf.kt_exp_chck(1,it-1);
    irf.kt_une_chck(1,it) = irf.y_une(1,it-1) - irf.c_une(1,it-1) + (1-par.delta)*irf.kt_une_chck(1,it-1);
end
clear it;


%% LINEARIZATION, eq. 1.21(using method of Klein)

AA     = [1,0;0,1];
par.rss = par.alpha*par.kss^(par.alpha-1) + 1 - par.delta; 
par.aux = par.css/(1+par.nn) * par.beta/par.sigma *par.alpha*(par.alpha-1)*par.kss^(par.alpha-2);
BB    = [1/(1+par.nn) * par.rss, -1/(1+par.nn);...
           par.aux*par.rss,1-par.aux];
        
       
[pol_c,pol_k] = solab(AA,BB,1); %  a*x(t+1) = b*x(t)

%% Calculate linearization for K_t starting from k_13
irf.kt_lin = NaN(1,par.TT+2);
%start in t=13  (column 14)
irf.kt_lin(1,par.Tz_end+2) = irf.kt_une(1,par.Tz_end+2); %initial value
for ic = par.Tz_end+3:par.TT+2; 
    irf.kt_lin(1,ic) = par.kss + pol_k*(irf.kt_lin(1,ic-1)-par.kss);    
end
clear ic;


%% Plot direct method & linearization

figure
plot([par.Tz_end+1:par.TT+1],irf.kt_une(1,par.Tz_end+2:par.TT+2));
hold all
plot([par.Tz_end+1:par.TT+1],irf.kt_lin(1,par.Tz_end+2:par.TT+2),'--');
hold all;
plot([par.Tz_end+1:par.TT+1],repmat(par.kss,1,par.TT-par.Tz_end+1),':k');

xlabel('period')
ylabel('capital stock');
legend('direct compuation','linear approximation','steady state');
axis([par.Tz_end+1 par.TT+1 -inf inf]);
axis 'auto y';


%% SAVE DATA
save Ch1_ramsey1;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SOlve equation 1.14 for expected shock
function [kt] = solve_1_14_exp(par)

% give pattern of jacobian, to prevent algorithm from estimating full
% jacobian
jac_pat = eye(par.TT) + diag(ones(par.TT-1,1),1) + diag(ones(par.TT-1,1),-1);

options = optimset('Display','off','TolFun',1e-12,'TolX',1e-12,'JacobPattern',sparse(jac_pat),'Algorithm','trust-region-reflective');

zt                          = par.zss*ones(1,par.TT);
zt(1,par.Tz_sta:par.Tz_end) = 1.1*par.zss;

kt = fsolve(@(kt)res_1_14_exp(par,zt,kt),repmat(par.kss,1,par.TT),options);

kt = [par.kss,kt,par.kss]; %period 0 and period par.TT + 1 added

end

%% Residual function 1.14 for expected shock
function [res] = res_1_14_exp(par,zt,kt)

res = NaN(par.TT,1);
for it = 1:par.TT;    
    if it == 1;
        k_p = par.kss; % k_t-1
        z_p = par.zss;
    else
        k_p = kt(1,it-1); % k_t-1
        z_p = zt(1,it-1);
    end
    if it == par.TT;
        k_n = par.kss; % k_t+1
    else
        k_n = kt(1,it+1);% k_t+1
    end
    
    cc      = consumption(par,zt(1,it),kt(1,it),k_n);
    c_p     = consumption(par,z_p,k_p,kt(1,it));

    res(it,1) = (cc / c_p)^par.sigma - par.beta*(par.alpha*zt(1,it)*kt(1,it)^(par.alpha - 1) + 1 -par.delta);
end

end

%% Solve equation 1.14 for unexpected shock
function [kt] = solve_1_14_une(par)

T_lng = par.TT - par.Tz_sta;

% give pattern of jacobian, to prevent algorithm from estimating full
% jacobian
jac_pat = eye(T_lng)+ diag(ones(T_lng-1,1),1) + diag(ones(T_lng-1,1),-1);

options = optimset('Display','off','Algorithm','trust-region-reflective','TolFun',1e-12,'TolX',1e-12,'JacobPattern',sparse(jac_pat));

zt                          = par.zss*ones(1,par.TT+1);
zt(1,par.Tz_sta:par.Tz_end) = 1.1*par.zss;

kt = fsolve(@(kt)res_1_14_une(par,zt,kt),repmat(par.kss,1,T_lng),options);

kt = [repmat(par.kss,1,par.Tz_sta+1),kt,par.kss];  %add t=0:t=Tz_sta, t=TT+1

end


%% Residual function 1.14 for unexpected shock
function [res] = res_1_14_une(par,zt_f,kt)
% zt_f should be lenght 1:par.TT+1 
%(x_f means variable x over t=1:par.TT+1, while kt is T_st:par.TT)

if size(zt_f,2) ~= par.TT+1
    error('length of vector zt_f and time periods not consistent')
end

T_lng   = size(kt,2);         % length of period under consideration
T_st    = par.TT - T_lng + 1;  % starting period;

kt_f    = [repmat(par.kss,1,par.TT-T_lng),kt,par.kss]; %kt over full length

if size(kt_f,2) ~= par.TT+1
    error('length of vectors kt_f and time periods not consistent')
end

cnt = 1;
res = NaN(T_lng,1);
for it = T_st:par.TT;
    z_p = zt_f(1,it-1); %z_t-1;
    k_p = kt_f(1,it-1); % k_t-1
    k_n = kt_f(1,it+1);% k_t+1
    
    cc      = consumption(par,zt_f(1,it),kt_f(1,it),k_n);
    c_p     = consumption(par,z_p,k_p,kt_f(1,it));
    
    res(cnt,1) = (cc / c_p)^par.sigma - par.beta*(par.alpha*kt_f(1,it)^(par.alpha - 1) + 1 -par.delta);    
    cnt = cnt  + 1;
end

end

%% Consumption function
function [cc] = consumption(par,zz,kk,k_n)

cc = production(par,zz,kk) + (1-par.delta)*kk - (1+par.nn)*k_n;

end


%% PRoduction function
function [yy] = production(par,zz,kk)

yy = zz.*kk.^par.alpha;

end

%% System 1.15 in k_t & k_t-1
function [f] = eq_1_15(par,kk,xx)

cc = consumption(par,1,xx,kk);

f(1,1) = 1/(1+par.nn) * (  production(par,1,kk) + (1-par.delta)*kk - ...
            ( par.beta*(par.alpha*kk^(par.alpha-1) + 1 - par.delta) )^(1/par.sigma)*...
            cc  );
f(2,1) = kk;
		
end

%% System 1.15 in k & c_t
function [f] = eq_1_17(par,kk,cc)

yy      = production(par,1,kk);
k_n     = 1/(1+par.nn) *(yy  + (1-par.delta)*kk - cc);

f(1,1)  = k_n;

f(2,1) = cc * (  par.beta*( par.alpha*k_n^(par.alpha-1) + 1- par.delta)  )^(1/par.sigma);

end

%% Analytical Jacobian
function [Jac] = jac_117_ana(par,kss,css)

Jac(1,1)    = 1/(1+par.nn)*(par.alpha*kss^(par.alpha-1)+1-par.delta);
Jac(1,2)    = -1/(1+par.nn);
Jac(2,1)    = css*par.beta*1/par.sigma*par.alpha*(par.alpha-1)/(1+par.nn)*kss^(par.alpha-2)*(par.alpha*kss^(par.alpha-1)+1-par.delta);
Jac(2,2)    = 1-css*par.beta/par.sigma*par.alpha*(par.alpha-1)/(1+par.nn)*kss^(par.alpha-2);

end
