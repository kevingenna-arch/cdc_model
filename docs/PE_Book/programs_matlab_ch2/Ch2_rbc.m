% Ch2_rbc.m
% Burkhard Heer, 23.10.2015
%
% computes the real business cycle model in Ch2
%

clear all;
tic

% parameters 
alpha=0.36;         % elasticity of output with respect to labor    @
eta=2.0;            % elasticity of maringal utility of consumpiton @
delta=0.02;        % rate of capital depreciation                  @
a=1.00;            % growth rate of output,4.5 percent per annum   @
nstar=0.30;         % one third of time endowment devoted to work    @
rho=0.95;           % autoregressive parameter of productivity shock@
sigma=0.0072;       % standard deviation of epsilon                 @
beta=0.99;
n=0;

% periods for impulse response functions
nimp=40;

% nofs artificial time series of length nobs
nobs=60;
nofs=500;

if beta>1
    disp('beta larger 1');
    break
end

% Stationary solution 
yk=((a^eta)-beta*(1-delta))/(beta*alpha);   % output-capital ratio      @
ck=yk+(1-delta-a);                              % consumption-capital ratio @
yc=yk/ck;         % inverse of average propensity to consume @
yi=yk/(a+delta-1);
theta=((1-alpha)*yc*(1-nstar))/nstar;               % theta implied by nstar          @
gamma=1/(1+theta);
kn=yk^(-1/(1-alpha));                               % capital-labor ratio             @
kstar=kn*nstar;                                 % stationary capital stock k=K/A  @
istar=(a+delta-1)*kstar;                        % stationary level of investment  @
cstar=yk*kstar - istar;                         % stationary level of consumption @
ystar=kstar^alpha*nstar^(1-alpha);
rstar=alpha*ystar/kstar;

nx=1;
Cu=zeros(2,2);
zeta1=(1-(1-gamma)*(1-eta))*nstar/(1-nstar)+alpha;
Cu=[-(1-gamma*(1-eta)) -(1-gamma)*(1-eta)*nstar/(1-nstar); 
gamma*(1-eta) zeta1];

Cxl=[0 1; alpha 1];

Cz=[0;1];

Dxl=[(1-alpha)*beta*rstar -1; (1+n)*kstar/ystar 0];

Fxl=[0 1; -(alpha+(1-delta)*kstar/ystar) 0];

Du=[0 (1-alpha)*beta*rstar;0 0];

Fu=[0 0; -cstar/ystar (1-alpha)];

Dz=[beta*rstar;0];

Fz=[0; 1];

[Lxx,Lxz,Llx,Llz,Lux,Luz] = SolveLA(Cu,Cxl,Cz,Dxl,Fxl,Du,Fu,Dz,Fz,rho,nx);

Lxx
Lxz
Llx
Llz
Lux
Luz


pause
% computation of the impulse response functions
% technology shock = 1 in period 2, zero thereafter


    ximp=zeros(nimp+1,1);
    zimp=ximp;
    zimp(2)=1;
    uimp=zeros(2,nimp);
    cimp=zeros(nimp,1);
    yimp=cimp;
    iimp=cimp;
    wimp=cimp;
    
    % compute time path
     for t=2:nimp
        zimp(t+1)=rho*zimp(t);
        ximp(t+1)=Lxx*ximp(t)+Lxz*zimp(t);
        uimp(:,t)= Lux*ximp(t) + Luz * zimp(t);
        yimp(t) = (1-alpha) * uimp(2,t) + alpha * ximp(t) + zimp(t);
        wimp(t) = yimp(t)-uimp(2,t);
        iimp(t) = yi*yimp(t)-(yi-1)*uimp(1,t);    % i_t = y_t -c_t
    end
    
    
    t2=(1:1:nimp)' ;
    zers=zeros(nimp,1) ;

    figure(1)

    subplot(221)
    plot(t2,[ zimp(1:nimp) zers],'k-' )
    title('z')
    set(gca,'XLim',[0 nimp]) ;

    subplot(222)
    plot(t2,[ ximp(1:nimp) zers ],'k-',t2,yimp,'k:')
    title('k and y')
    set(gca,'XLim',[0 nimp]) ;

    subplot(223)
    plot(t2,[ uimp(1,:)' zers ],'k-',t2,iimp,'k:')
    title('c and i')
    set(gca,'XLim',[0 nimp]) ;

    subplot(224)
    plot(t2,[ uimp(2,:)' zers ],'k-',t2,wimp,'k:')
    title('n and w')
    set(gca,'XLim',[0 nimp]) ;



pause

% time series simulation
% output, investment, consumption, hours, real wage

sx=zeros(1,5);  % mean standard deviation of y, i, c, n, w
rxy=zeros(1,4);
rxx=0;

eps=randn(nobs+1,nofs);
eps=sigma*eps;

for itn=1:nofs
    disp(['simulation number  ', num2str(itn)]);
    xt=zeros(nobs+1,1);
    zt=xt;
    zt(1)=eps(1,itn);
    ut=zeros(2,nobs);
    ct=zeros(nobs,1);
    yt=ct;
    it=ct;
    wt=ct;
    
    % compute time path
    for t=1:nobs
        zt(t+1)=rho*zt(t)+eps(t,itn);
        xt(t+1)=Lxx*xt(t)+Lxz*zt(t);
        ut(:,t)= Lux*xt(t) + Luz * zt(t);
        yt(t) = (1-alpha) * ut(2,t) + alpha * xt(t) + zt(t);
        wt(t) = yt(t)-ut(2,t);
        it(t) = yi*yt(t)-(yi-1)*ut(1,t);    % i_t = y_t -c_t
    end
    
%    t1=1:1:nobs;
%    figure
%    plot(t1,xt(t1));
%    pause
%    figure
%    plot(t1,ut);
%    pause
%    figure
%    plot(t1,yt);
%    pause
%    figure
%    plot(t1,it);
%    pause
    
    ct1=hpfilter(ut(1,:)',1600);
    nt1=hpfilter(ut(2,:)',1600);
    yt1=hpfilter(yt,1600);
    it1=hpfilter(it,1600);
    wt1=hpfilter(wt,1600);
    ct=ct1-ut(1,:)';
    nt=nt1-ut(2,:)';
    yt=yt1-yt;
    it=it1-it;
    wt=wt1-wt;
    st=[yt it ct nt wt];
 
%    t1=1:1:nobs;
%    figure
%    plot(t1,ct, t1,nt, t1,yt, t1,it, t1,wt);
%   pause
%    std(st)
%    corrcoef(st)

    sx=sx+std(st);
    temp=corrcoef(st);
    rxy=rxy+temp(1,2:5);
    
    
    temp=[yt(1:nobs-1) yt(2:nobs)];
    temp=corrcoef(temp);
    rxx=rxx+temp(1,2);
    
end

sx=sx/nofs;
rxy=rxy/nofs;
rxx=rxx/nofs;

disp('Standard deviation of y, i, c, n, w');
sx
disp('Correlation of i, c, n, w with output y');
rxy
disp('Autocorrelation of output');
rxx

