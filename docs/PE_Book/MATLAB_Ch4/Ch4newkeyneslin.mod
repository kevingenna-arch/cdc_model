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
varexo e_z, e_g, e_Q;

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
    xi = phi* cp^(1-1/rhoC)+(1-phi)*g^(1-1/rhoC);
    ce = xi^(1/(1-1/rhoC));
    lamb = phi*((ce-chi*ch)^(-sigmaC))*xi^(1/(1-1/rhoC)-1)*cp^(-1/rhoC);
    lamb = beta*lamb(+1)*bigq/infl(+1);
    ch = ce(-1);

    %% Firm
    ytil = sy*y;
    w = (1-alpha) * mc * ytil/l;
    ytil = exp(z)*(k(-1)^alpha)*(l^(1-alpha));
    q = (1/a1)*((i/k(-1))^zeta);
    pa = (epsy/(epsy-1))*(gam1/gam2);
    k = ((a1/(1-zeta))*(i/k(-1))^(1-zeta)+a2)*k(-1)+(1-delta)*k(-1);
    (1-phiy)*(pa^(1-epsy)) = 1- phiy*((inflp/infl)^(1-epsy));
    wa = (epsn/(epsn-1))*(psi1/psi2);
    w^(1-epsn) = (1-phin)*(wa^(1-epsn)) + phin*((inflp/infl)*wp)^(1-epsn);
    sy = (1-phiy)*(pa^(-epsy)) + phiy*((inflp/infl)^(-epsy))*sy(-1);
    y = cp + i + g;
    q = beta*(lamb(+1)/lamb)*(alpha*mc(+1)*(ytil(+1)/k)-(i(+1)/k)+q(+1)*(1-delta+((a1/(1-zeta))*(i(+1)/k)^(1-zeta)+a2)));
    gam1 = (mc*lamb*y)+(beta*phiy)*(infl/infl(+1))^(-epsy)*gam1(+1);
    gam2 = (lamb*y)+(beta*phiy) *(infl/infl(+1))^(1-epsy)*gam2(+1);
    psi1 = (l^(1+nu1))*nu0*(wa/w)^(-epsn*(1+nu1))+(beta*phin) *(infl*wa/(infl(+1)*wa(+1)))^(-epsn*(1+nu1)) *psi1(+1);
    psi2 = lamb*(wa/w)^(-epsn)*l  + (beta*phin)* (wa/wa(+1))^(-epsn)  *(infl/infl(+1))^(1-epsn) *psi2(+1);
    markup = 1/mc;

    %% Government consumption
    g = gstar^(1-rho_g)*g(-1)^(rho_g)*exp(geps);	

    %% Interest rate rule
    bigqp = ((inflstar/beta)^(1-delta1))*(bigq^delta1)*((infl/inflstar)^delta2)*(y/ystar)^(delta4y)*exp(qeps);
    inflp = infl(-1); 
    bigq = bigqp(+1);
    wp = w(-1);
    syp = sy(-1);

    % Exogenous shocks
	z = rho_z*z(-1) + e_z;
	geps = e_g;
    qeps = e_Q;

end;


%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
    y = Y;
    ce = CE;
    cp = CP;
    ch = CE;
    i = I;
    l = L; 
    w = W;  
    g = G;
    k = K;
    xi = XI;
    z = 0;
    q = Q;
    infl = INFL;
    inflp = INFL;
    ytil = Y;
    sy = 1;
    bigq = BIGQ;
    mc = MC;
    wa = w;
    lamb = LAMB;
    gam1 = GAM1;
    gam2 = GAM2;
    psi1 = PSI1;
    psi2 = PSI2;
    pa = PA;
    markup = 1/MC; 
    bigqp = BIGQ;
    wp = W; 
    syp = 1;

    geps = 0;
    qeps = 0;

end;

shocks;
var e_z;  stderr .0072;
var e_g;  stderr .01;
var e_Q;  stderr .02520;
end;

check;

stoch_simul(order=1,irf=20);
