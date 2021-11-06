/*
 * Ch4rbc2.mod with DYNARE:
 * 
 * computes the model with stochastic government consumption
 * in Chapter 4.4 of Heer (2018) "Public Economics: The Macroeconomic
 * Perspective", with adaptive separable utility function
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
			phi rho_z rho_g psi0 psi1 gstar;

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
psi1        = 0.3;              % Frisch labor supply elasticity

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
LAMB        = CE^(-sigmaC)*phi*XI^(1/(1-1/rhoC)-1)*CP^(-1/rhoC);
psi0        = W*LAMB/(L^(1/psi1));
gstar       = G;

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
model;
    %% Household
    xi = phi* cp^(1-1/rhoC)+(1-phi)*g^(1-1/rhoC);
    ce = xi^(1/(1-1/rhoC));
	% Euler
	ce^(-sigmaC)*xi^(1/(1-1/rhoC)-1)*cp^(-1/rhoC) = beta*ce(+1)^(-sigmaC)*xi(+1)^(1/(1-1/rhoC)-1)*cp(+1)^(-1/rhoC)*(1+r(+1)-delta);
	% hours supply
	w*ce^(-sigmaC)*phi*xi^(1/(1-1/rhoC)-1)*cp^(-1/rhoC) = psi0*l^(1/psi1);

    %% Firm
    w = (1-alpha) * y/l;
    y = exp(z)*(k(-1)^alpha)*(l^(1-alpha));
    k = i+(1-delta)*k(-1);
    r = alpha * y/k(-1);
    y = cp + i + g;

    %% Government consumption
    g = gstar^(1-rho_g)*g(-1)^(rho_g)*exp(geps);		// government consumption in period t

    % Exogenous shocks
	z = rho_z*z(-1) + e_z;
	geps = e_g;
end;


%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
    y = Y;
    ce = CE;
    cp = CP;
    i = I;
    l = L; 
    w = W; 
    r = R; 
    g = G;
    k = K;
    xi = XI;
    z = 0;
    geps = 0;
end;

shocks;
var e_z;  stderr .0072;
var e_g;  stderr .01;
end;

stoch_simul(order=1,irf=20);
