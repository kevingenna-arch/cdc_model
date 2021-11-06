/*
 * Ch4rbc1.mod with DYNARE:
 * 
 * computes the model with stochastic government consumption
 * in Chapter 4.4 of Heer (2018) "Public Economics: The Macroeconomic
 * Perspective"
 *
 * Author: Burkhard Heer
 * Last Change: June 11, 2018
 *
 * To run the program, you need to install DYNARE on your PC
 * see: http:\\www.dynare.org
 */


%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------


var y, ce, cp, i, l, w, r, g, geps, k, z, xi;
varexo e_z, e_g;

parameters	beta alpha delta rhoC sigmaC gy
			phi rho_z rho_g gam gstar;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
% parameters
alpha   	= 0.36;				% share of capital in ouput
beta    	= 0.99;				% discount factor
delta  		= 0.025;			% depreciation of capital
rhoC		= 0.5;				% substitution elasticity of private/public
                                % consumption
sigmaC		= 2;				% intertemporal elasticity of substitution
gy 			= 0.2; 			 	% Public spending in GDP
phi         = 0.75;             % relative weight of private consumption in C
% shock process
rho_z   	= 0.95; 			% productivity 
rho_g   	= 0.90; 			% public spending

% steady states
L           = 0.3;
R1			= 1/beta;
R			= 1/beta-(1-delta);
KL          = (alpha/R)^(1/(1-alpha));   % caital-labor ratio
K           = KL*L;
Y			= K^alpha*L^(1-alpha);
CP			= (1-gy)*Y-delta*K;
I			= delta*K;
W			= (1-alpha)*Y/L;
G           = gy*Y;
XI           = phi*CP^(1-1/rhoC)+(1-phi)*G^(1-1/rhoC);
CE           = XI^(1/(1-1/rhoC));
temp        = (1-L)*W/CE*phi*XI^(1/(1-1/rhoC)-1)*CP^(-1/rhoC);
gam         = 1/(1+temp);
LAMB        = gam*CE^(gam*(1-sigmaC)-1)*(1-L)^((1-gam)*(1-sigmaC));
LAMB        = LAMB*phi*XI^(1/(1-1/rhoC)-1)*CP^(-1/rhoC);
gstar       = log(G);

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
model;
    %% Household
    exp(xi) = phi* exp(cp)^(1-1/rhoC)+(1-phi)*exp(g)^(1-1/rhoC);
    exp(ce) = exp(xi)^(1/(1-1/rhoC));
	% Euler
	exp(ce)^(gam*(1-sigmaC)-1)*(1-exp(l))^((1-gam)*(1-sigmaC))*exp(xi)^(1/(1-1/rhoC)-1)*exp(cp)^(-1/rhoC) = beta*exp(ce(+1))^(gam*(1-sigmaC)-1)*(1-exp(l(+1)))^((1-gam)*(1-sigmaC))*exp(xi(+1))^(1/(1-1/rhoC)-1)*exp(cp(+1))^(-1/rhoC)*(1+exp(r(+1))-delta);
	% hours supply
	exp(w)*gam*exp(ce)^(gam*(1-sigmaC)-1)*(1-exp(l))^((1-gam)*(1-sigmaC))*phi*exp(xi)^(1/(1-1/rhoC)-1)*exp(cp)^(-1/rhoC) = (1-gam)*exp(ce)^(gam*(1-sigmaC))*(1-exp(l))^((1-gam)*(1-sigmaC)-1);

    %% Firm
    exp(w) = (1-alpha) * exp(y)/exp(l);
    exp(y) = exp(z)*(exp(k(-1))^alpha)*(exp(l)^(1-alpha));
    exp(k) = exp(i)+(1-delta)*exp(k(-1));
    exp(r) = alpha * exp(y)/exp(k(-1));
    exp(y) = exp(cp) + exp(i) + exp(g);

    %% Government consumption
    exp(g) = exp(gstar)^(1-rho_g)*exp(g(-1))^(rho_g)*exp(geps);		// government consumption in period t

    % Exogenous shocks
	z = rho_z*z(-1) + e_z;
	geps = e_g;
end;


%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
    y = log(Y);
    ce = log(CE);
    cp = log(CP);
    i = log(I);
    l = log(L); 
    w = log(W); 
    r = log(R); 
    g = log(G);
    k = log(K);
    xi = log(XI);
    z = 0;
    geps = 0;
end;

shocks;
var e_z;  stderr .0072;
var e_g;  stderr .01;
end;

stoch_simul(order=1,irf=20);
