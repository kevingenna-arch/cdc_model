/*
 * Ch4newkeyneslin.mod with DYNARE:
 * 
 * computes the New Keynesian model 
 * in Chapter 4.5 of Heer (2018) "Public Economics: The Macroeconomic
 * Perspective"
 *
 * Author: Burkhard Heer
 * Last Change: June 12, 2018
 *
 * To run the program, you need to install DYNARE on your PC
 * see: http:\\www.dynare.org
 */


%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------


var y, ytil, ce, cp, i, sy, q, mc, pa, wa, l, infl, gam1, gam2, psi1, psi2, lamb, ch, bigq, w, inflp, g, geps, qeps, k, z, xi, markup, bigqp, wp, syp;
varexo e_z, e_g, e_q;

parameters	beta alpha delta rhoC sigmaC gy chi nu1 nu0
            epsn epsy phin phiy zeta delta1 delta2 delta4y 
            a1 a2
			phi rho_z rho_g gstar inflstar ystar;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
% parameters
alpha   	= 0.36;				% share of capital in ouput
beta    	= 0.99;				% discount factor
chi         = 0.65;             % habit parameter consumption
sigmaC		= 2;				% intertemporal elasticity of substitution
nu1         = 5.0;              % inverse of the Frisch labor supply elasticity
rhoC		= 0.5;				% substitution elasticity of private/public
                                % consumption

delta  		= 0.025;			% depreciation of capital
epsn        = 6.0;              % wage elasticity labor
epsy        = 6.0;              % price elasticity goods

phin        = 0.5;              % share of wages that adjust optimally
phiy        = 0.5;              % share of prices that adjust optimally
gy 			= 0.2; 			 	% Public spending in GDP
phi         = 1.0;             % relative weight of private consumption in C

% adjustment cost function
zeta        = 3.0;              % elasticity adjustment costs
a1          = delta^zeta;
a2          = ((-zeta)/(1.00-zeta))*delta;

% Taylor rule
delta2      = 1.50;             % coefficient on inflation
delta1      = 0.90;             % autoregressive parameter
delta4y     = 0.25;             % coefficient on output             


% shock process
rho_z   	= 0.95; 			% productivity 
rho_g   	= 0.90; 			% public spending

% steady states
INFL        = 1.005;            % inflation factor
L           = 0.3;
Q           = 1.0;
PA          = 1.0;
SY          = 1.0;
BIGQ        = INFL/beta;        % nominal interest rate facctor
MC          = (epsy-1)/epsy;
YK          = (1-beta*(1-delta))/(alpha*beta*MC);
KL          = YK^(1/(alpha-1));
K           = KL*L;
Y           = YK*K;
I			= delta*K;
G           = gy*Y;
CP			= (1-gy)*Y-delta*K;
W			= (1-alpha)*MC*Y/L;
WA          = W;
XI           = phi*CP^(1-1/rhoC)+(1-phi)*G^(1-1/rhoC);
CE           = XI^(1/(1-1/rhoC));
LAMB        = phi*((CE-chi*CE)^(-sigmaC))*XI^(1/(1-1/rhoC)-1)*CP^(-1/rhoC);

GAM1=(MC*LAMB*Y)/(1-beta*phiy);
GAM2=GAM1/MC;

nu0=((epsn-1)/epsn)*W*LAMB*(L^(-nu1));
PSI1=(nu0*L^(1+nu1))/(1-beta*phin);
PSI2=(LAMB*L)/(1-beta*phin);

gstar       = G;
inflstar    = INFL;
ystar       = Y;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
model;
    %% Household
    exp(xi) = phi* exp(cp)^(1-1/rhoC)+(1-phi)*exp(g)^(1-1/rhoC);
    exp(ce) = exp(xi)^(1/(1-1/rhoC));
    exp(lamb) = phi*((exp(ce)-chi*exp(ch))^(-sigmaC))*exp(xi)^(1/(1-1/rhoC)-1)*exp(cp)^(-1/rhoC);
    exp(lamb) = beta*exp(lamb(+1))*exp(bigq)/exp(infl(+1));
    exp(ch) = exp(ce(-1));

    %% Firm
    exp(ytil) = exp(sy)*exp(y);
    exp(w) = (1-alpha) * exp(mc) * exp(ytil)/exp(l);
    exp(ytil) = exp(z)*(exp(k(-1))^alpha)*(exp(l)^(1-alpha));
    exp(q) = (1/a1)*((exp(i)/exp(k(-1)))^zeta);
    exp(pa) = (epsy/(epsy-1))*(exp(gam1)/exp(gam2));
    exp(k) = ((a1/(1-zeta))*(exp(i)/exp(k(-1)))^(1-zeta)+a2)*exp(k(-1))+(1-delta)*exp(k(-1));
    (1-phiy)*(exp(pa)^(1-epsy)) = 1- phiy*((exp(inflp)/exp(infl))^(1-epsy));
    exp(wa) = (epsn/(epsn-1))*(exp(psi1)/exp(psi2));
    exp(w)^(1-epsn) = (1-phin)*(exp(wa)^(1-epsn)) + phin*((exp(inflp)/exp(infl))*exp(wp))^(1-epsn);
    exp(sy) = (1-phiy)*(exp(pa)^(-epsy)) + phiy*((exp(inflp)/exp(infl))^(-epsy))*exp(sy(-1));
    exp(y) = exp(cp) + exp(i) + exp(g);
    exp(q) = beta*(exp(lamb(+1))/exp(lamb))*(alpha*exp(mc(+1))*(exp(ytil(+1))/exp(k))-(exp(i(+1))/exp(k))+exp(q(+1))*(1-delta+((a1/(1-zeta))*(exp(i(+1))/exp(k))^(1-zeta)+a2)));
    exp(gam1) = (exp(mc)*exp(lamb)*exp(y))+(beta*phiy)*(exp(infl)/exp(infl(+1)))^(-epsy)*exp(gam1(+1));
    exp(gam2) = (exp(lamb)*exp(y))+(beta*phiy) *(exp(infl)/exp(infl(+1)))^(1-epsy)*exp(gam2(+1));
    exp(psi1) = (exp(l)^(1+nu1))*nu0*(exp(wa)/exp(w))^(-epsn*(1+nu1))+(beta*phin) *(exp(infl)*exp(wa)/(exp(infl(+1))*exp(wa(+1))))^(-epsn*(1+nu1)) *exp(psi1(+1));
    exp(psi2) = exp(lamb)*(exp(wa)/exp(w))^(-epsn)*exp(l)  + (beta*phin)* (exp(wa)/exp(wa(+1)))^(-epsn)  *(exp(infl)/exp(infl(+1)))^(1-epsn) *exp(psi2(+1));
    exp(markup) = 1/exp(mc);

    %% Government consumption
    exp(g) = gstar^(1-rho_g)*exp(g(-1))^(rho_g)*exp(geps);	

    %% Interest rate rule
    exp(bigq) = ((inflstar/beta)^(1-delta1))*(exp(bigqp)^delta1)*((exp(infl)/inflstar)^delta2)*(exp(y)/ystar)^(delta4y)*exp(qeps);
    exp(inflp) = exp(infl(-1)); 
    exp(bigq(-1)) = exp(bigqp);
    exp(wp) = exp(w(-1));
    exp(syp) = exp(sy(-1));

    % Exogenous shocks
	z = rho_z*z(-1) + e_z;
	geps = e_g;
    qeps = e_q;

end;


%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
    y = log(Y);
    ce = log(CE);
    cp = log(CP);
    ch = log(CE);
    i = log(I);
    l = log(L); 
    w = log(W);  
    g = log(G);
    k = log(K);
    xi = log(XI);
    z = 0;
    q = log(Q);
    infl = log(INFL);
    inflp = log(INFL);
    ytil = log(Y);
    sy = log(1);
    bigq = log(BIGQ);
    mc = log(MC);
    wa = log(W);
    lamb = log(LAMB);
    gam1 = log(GAM1);
    gam2 = log(GAM2);
    psi1 = log(PSI1);
    psi2 = log(PSI2);
    pa = log(PA);
    markup = log(1/MC); 
    bigqp = log(BIGQ);
    wp = log(W); 
    syp = log(1);

    geps = 0;
    qeps = 0;

end;

shocks;
var e_z;  stderr .0072;
var e_g;  stderr .01;
var e_q;  stderr .02520;
end;

%check;

stoch_simul(order=1,irf=20);
