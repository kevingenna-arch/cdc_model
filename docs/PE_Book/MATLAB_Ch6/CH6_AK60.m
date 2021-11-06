function Ch6_AK60();
% 
%   Auerbach-Kotlikoff Model 
%
%   13.5.2020
%    
%    author: Burkhard Heer
%
%    algorithm 9.1 from Heer/Maussner, Dynamic General Equilibrium Modeling
%
%    direct computation of the OLG model in section 9.1
%

clc;
close all;

%
% Step 1: Parameterization 
%

beta0=0.96;         % discount factor 
r=0.045;            % initial value of the interest rate 
sigma=2;            % coefficient of relative risk aversion 
alpha=0.36;         % production elasticity of capital 
rep=0.3;            % replacement ratio 
delta=0.1;          % rate of depreciation 
tr=20;              % retired 
t=40;               % working time 
tau=rep/(2+rep);    % income tax rate 
gam=2;              % disutility from working 

psi0=0.001;         % parameter of utility function 
phi=0.8;            % update parameter of aggregates
tol=0.001;          % percentage deviation of final solution 
tolk=0.001;         % percentage deviation of final solution for a^1 
nq1=30;             % maximum number of iterations over a^60
nq=30;              %  maximum number of iterations over K,L

par.psi0=psi0;
par.gam=gam;
par.sigma=sigma;
par.alpha=alpha;
par.beta0=beta0;
par.delta=delta;
par.tr=tr;
par.t=t;
par.tau=tau;

% ---------------------------------
%
%   initialization of K, L, tau
%
% --------------------------------- 

nbar=0.2;
kbar=(alpha/(r+delta))^(1/(1-alpha))*nbar;
kold=100;
nold=2;


% agents' policy function  
aopt=zeros(t+tr,1);     % optimal asset 
copt=zeros(t+tr,1);     % optimal consumption 
nopt=0.3*ones(t,1);     % optimal labor supply 
par.nopt=nopt;

% --------------------------------------------------
%
%   iteration of policy function, wealth distribution
%
% -------------------------------------------------- 

q=0;
krit=1+tol;
while (q<nq) | (krit>tol)       % iteration over K,L while deviation is too large
   
    clc;
    q=q+1
    krit
    
    w=wage(kbar,nbar,par);
    r=interest(kbar,nbar,par);
    pen=rep*(1-tau)*w*nbar*3/2;
    par.w=w;
    par.r=r;
    par.pen=pen;
    
    k60q=zeros(nq1,1);
    k1q=zeros(nq1,1);

    q1=0;
    while (q1<nq1)| ((q1<5) | abs(aopt(1))>tolk)
        q1=q1+1;
        disp('iteration over q1: ');
        q1
        
        if q1==1
            k60=0.1;
        elseif q1==2
            k60=0.2;
        else    % secant method, see Heer/Maussner, 2nd ed., p. 609
            k60=k60q(q1-1)-(k60q(q1-1)-k60q(q1-2))/(k1q(q1-1)-k1q(q1-2))*k1q(q1-1); 
        end

        aopt(t+tr) = k60;
        copt(t+tr) = pen+(1+r)*k60;
        k60q(q1)=k60;
        par.aopt=aopt;

        % computation of the decision rules for the retired 
        i=tr;
        while i>1      % all periods t=T+1,T+2,..T+TR 
            i=i-1;
            par.i=i;
            disp('old worker: ');
            i+t
            x0=aopt(i+t+1); % initial guess for a^s in Euler condition
                            % given a^{s+1}, a^{s+2}
                            % simply: a^s=a^{s+1}

            a0 = fsolve(@(x)rfold(x,par),x0)
            aopt(i+t)=a0;
            copt(i+t)=pen+(1+r)*aopt(i+t)-aopt(i+t+1);
            par.aopt=aopt;
        end


        % compuation of the decsion rules for the worker 
        i=t+1;
        while i>1       % all periods t=1,2,..T 
            i=i-1;
            par.i=i;
            disp('young worker: ');
            i
            kbar
            x0=[aopt(i+1),nopt(i)];
            
            xf = fsolve(@(x)rfyoung(x,par),x0);
            aopt(i)=xf(1);
            nopt(i)=xf(2);
            copt(i)=nopt(i)*w*(1-tau)+(1+r)*aopt(i)-aopt(i+1);
            par.aopt=aopt;
            par.nopt=nopt;
        end   
    
      k1q(q1)=aopt(1);
      %  k1q(1:q1) 
        

    end     % q1 - iteration of k^1,..,k^60 */

    % computation of the aggregate capital stock and employment nbar 
    knew=mean(aopt);
    nnew=mean(nopt)*2/3;    
    krit=abs((kbar-knew)/kbar);
    krit0=abs((nbar-nnew)/nbar);
    kbar=phi*kbar+(1-phi)*knew
    nbar=phi*nbar+(1-phi)*nnew
end   %  iteration over K,L

disp('solution');
kbar
nbar
aopt(1)
pause;

periods=linspace(20,20+t+tr-1,t+tr);
periods_workers=linspace(20,20+t-1,t);

figure
plot(periods,aopt);
xlabel('Real-fife age');
ylabel('Individual wealth');
pause;


figure
plot(periods,copt);
xlabel('Real-fife age');
ylabel('Individual consumption');
pause;

figure
plot(periods_workers,nopt);
xlabel('Real-fife age');
ylabel('Individual labor supply');
pause;

save CH6_AK60.mat;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% AUXILIARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   u           -- utility function
%   uc          -- marginal utility of consumption
%   un          -- marginal utility of leisure

function [y] = u(c,l,par)

    if par.sigma==1
        y=ln(c+par.psi0)+par.gam*ln(1-l);
    else
        y=(((c+par.psi0)*(1-l)^par.gam)^(1-par.sigma)-1)/(1-par.sigma);
    end
end

function [y] = uc(c,l,par)
    y=(c+par.psi0)^(-par.sigma).*(1-l)^(par.gam*(1-par.sigma));
end

function [y] = un(c,l,par)
    y=par.gam*(c+par.psi0)^(1-par.sigma).*(1-l)^(par.gam*(1-par.sigma)-1);
end

function [y] = rfold(x,par)
    aopt=par.aopt;
    i=par.i;
    t=par.t;
    tr=par.tr;
    r=par.r;
    k0=x(1);
    k1=aopt(i+t+1);
    pen=par.pen;
    beta0=par.beta0;
    if i==tr-1
        k2=0;
    else
        k2=aopt(i+t+2);
    end
    c0=(1+r)*k0+pen-k1;
    c1=(1+r)*k1+pen-k2;
    y=uc(c0,0,par)/beta0-uc(c1,0,par)*(1+r);
end


function [y] = rfyoung(x,par)
    rf1=0;
    rf2=0;
    k0=x(1);
    n0=x(2);
    
    aopt=par.aopt;
    nopt=par.nopt;
    i=par.i;
    t=par.t;
    tr=par.tr;
    r=par.r;
    pen=par.pen;
    beta0=par.beta0;
    tau=par.tau;
    w=par.w;
    
    k1=aopt(i+1);
    k2=aopt(i+2);
    if i==t
        n1=0;
        c1=(1+r)*k1+pen-k2;
    else
        n1=nopt(i+1);
        c1=(1+r)*k1+(1-tau)*w*n1-k2;
    end
    c0=(1+r)*k0+(1-tau)*w*n0-k1;
    rf1=uc(c0,n0,par)/beta0-uc(c1,n1,par)*(1+r);
    rf2=(1-tau)*w*uc(c0,n0,par)-un(c0,n0,par);
    y=[rf1,rf2];
end



%% interest rate 
function [y] = interest(k,l,par)
    alpha = par.alpha;
	y=alpha*k^(alpha-1)*l^(1-alpha)-par.delta;
end

%% wage 
function [w] = wage(k,l,par)
    alpha = par.alpha;
	w=(1-alpha)*k^alpha*l^(-alpha);
end
