/*
 * Ch5stochtaxlin.mod with DYNARE:
 * 
 * computes the model with stochastic labor and capital income taxes
 * in Chapter 5.7 of Heer (2018) "Public Economics: The Macroeconomic
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


var y, ce, cp, i, l, w, r, g, geps, k, z, xi, tauk, taukp, taul, taulp, taulpp, taukeps, tauleps;
varexo e_z, e_g, e_taul, e_tauk;

parameters	beta alpha delta rhoC sigmaC gy
			phi rho_z rho_g gam gstar taulstar taukstar rhok rhol1 rhol2;

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
rhok        = 0.9725;
rhol1       = 0.7841;
rhol2       = 0.2047;

% steady states
L           = 0.3;
TAUK        = 0.32;	% values from 2008.iv from Gomme et al, RED 2011
TAUL        = 0.28;
R1			= 1/beta;
R			= (1/beta-1)/(1-TAUK)+delta;
KL          = (alpha/R)^(1/(1-alpha));   % caital-labor ratio
K           = KL*L;
Y			= K^alpha*L^(1-alpha);
CP			= (1-gy)*Y-delta*K;
I			= delta*K;
W			= (1-alpha)*Y/L;
G           = gy*Y;
XI           = phi*CP^(1-1/rhoC)+(1-phi)*G^(1-1/rhoC);
CE           = XI^(1/(1-1/rhoC));
temp        = (1-L)*(1-TAUL)*W/CE*phi*XI^(1/(1-1/rhoC)-1)*CP^(-1/rhoC);
gam         = 1/(1+temp);
LAMB        = gam*CE^(gam*(1-sigmaC)-1)*(1-L)^((1-gam)*(1-sigmaC));
LAMB        = LAMB*phi*XI^(1/(1-1/rhoC)-1)*CP^(-1/rhoC);
gstar       = G;
taulstar    = TAUL;
taukstar    = TAUK;


%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
model;
    %% Household
    xi = phi* cp^(1-1/rhoC)+(1-phi)*g^(1-1/rhoC);
    ce = xi^(1/(1-1/rhoC));

	% Euler
	ce^(gam*(1-sigmaC)-1)*(1-l)^((1-gam)*(1-sigmaC))*xi^(1/(1-1/rhoC)-1)*cp^(-1/rhoC) = beta*ce(+1)^(gam*(1-sigmaC)-1)*(1-l(+1))^((1-gam)*(1-sigmaC))*xi(+1)^(1/(1-1/rhoC)-1)*cp(+1)^(-1/rhoC)*(1+(1-tauk(+1))*(r(+1)-delta));
	% hours supply
	w*(1-taul)*gam*ce^(gam*(1-sigmaC)-1)*(1-l)^((1-gam)*(1-sigmaC))*phi*xi^(1/(1-1/rhoC)-1)*cp^(-1/rhoC) = (1-gam)*ce^(gam*(1-sigmaC))*(1-l)^((1-gam)*(1-sigmaC)-1);

    %% Firm
    w = (1-alpha) * y/l;
    y = exp(z)*(k(-1)^alpha)*(l^(1-alpha));
    k = i+(1-delta)*k(-1);
    r = alpha * y/k(-1);
    y = cp + i + g;

    %% Government consumption
    g = gstar^(1-rho_g)*g(-1)^(rho_g)*exp(geps);		// government consumption in period t


	taul = taulstar^(1-rhol1-rhol2)*taulp^(rhol1)*taulpp^(rhol2)*exp(tauleps);
	tauk = taukstar^(1-rhok)*taukp^(rhok)*exp(taukeps);

    taulp = taul(-1);
    taulpp = taulp(-1);
    taukp = tauk(-1);

    % Exogenous shocks
	z = rho_z*z(-1) + e_z;
	geps = e_g;
    taukeps = e_tauk;
    tauleps = e_taul;
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
    taukeps = 0;
    tauleps = 0;
    taukp = TAUK;
    taulp = TAUL;
    taulpp = TAUL;
end;

shocks;
var e_z;  stderr .0072;
var e_g;  stderr .01;
var e_tauk;  stderr .027;
var e_taul;  stderr .028;
end;

stoch_simul(order=1,irf=20);
