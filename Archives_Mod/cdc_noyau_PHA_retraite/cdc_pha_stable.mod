%%%%% CDC clean ################################################################
% version stable qui produit un ES comme initval

% Set up vars for loops
@#define NE=1 					//%  ages of edu
@#define NT=work                   //%  working ages
@#define NTr=pens                  //%  retirement ages
@#define NLS=NT+NTr             //%  total ages
@#define IR=5                   //%  retirement wage indexation ages
@#define qualif = ["Q","NQ"]    //%  skill levels
@#include "../matrices_chocs.m" 	//%  load up shocks externally

%%%%% ENDOGENOUS VARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% endogenous vars from agents' maxing programs, by skill level
@#for ql in qualif
	% survival discount factor between periods-- CHECK
	@#for i in 2:NLS
		var beta_@{ql}_@{i};
	@#endfor
	% health stock -- CHECK
	@#for i in 2:NLS
		var h_@{ql}_@{i};
	@#endfor
	% cohort populations
	@#for i in 1:NLS
		var P_@{ql}_@{i} ;
	@#endfor
	% aggregates by skill class:
	% -	L: aggregate labour
	% - pi: share of first gen going for edu -- mirroring 
	var L_@{ql} pi_@{ql};

@#endfor

var y nbar Ptot Pret H;
var cCheck;

%%%%% EXOGENOUS VARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo A A_Q xpop shareq;

@#for ql in qualif
	% health stock or shocks?
	varexo h_@{ql}_1;

	@#for i in 2:NLS
		% shocks for migration and health stock
		varexo mig_@{ql}_@{i} med_@{ql}_@{i};
	@#endfor

@#endfor

%%%%% Parameters and values
parameters alp, beta, delta, Tr, T, LS, gam, gam1, deltah, phi, eta, rho;

@#for ql in qualif
	@#for i in 1:NT
		% a: worker productivity by age
		% e: switcher for skill level
		parameters a_@{ql}_@{i} e_@{ql}_@{i};    
	@#endfor
@#endfor

@#for i in 1:NT
	% switcher for education
	e_NQ_@{i}=0;
	@#if i <= NE
		e_Q_@{i}=1;
	@#else 
		e_Q_@{i}=0;
	@#endif
@#endfor

@#ifdef prodind
		@#if prodind == 1
			@#include "prodind_baseline.m"
		@#elseif prodind == 2
			@#include "prodind_formcont.m"
		@#elseif prodind == 3
			@#include "prodind_socact.m"
		@#endif
	@#else
	% fall back to standard
	@#include "prodind_current.m"
@#endif

gam     =   1.5;        %
gam1    =   .07;        %
alp     =   .7;         % 
eta     =   .5;         %
rho     =   .5;         % 
beta    =   0.97;       % discount factor
delta   =   0.02;       % physical capital depreciation
deltah  =   0.02;       % health depreciation
phi     =   .2;         % productivity effect for health
% phi     =   0;         % productivity effect for health
T       =   @{NT};          % working ages
Tr      =   @{NTr};          % retirement ages
LS      =   T+Tr;       % total ages

%%%%% MODEL SPECIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model;


pi_Q	=	shareq;
pi_NQ	=	1-pi_Q;

% double equations, blocs according to skill level
@#for ql in qualif

%%% Demographics
	
	% starting population by skill level
	P_@{ql}_1=pi_@{ql}*xpop;

	@#for i in 2:NLS
		% second cohort onwards, population dynamics
		P_@{ql}_@{i} =  beta_@{ql}_@{i}(-1)*P_@{ql}_@{i-1}(-1) + mig_@{ql}_@{i};
		% health stock dynamics
		h_@{ql}_@{i} =  (1-deltah)*h_@{ql}_@{i-1}(-1) + med_@{ql}_@{i};
		% health dynamics and survival probs
		h_@{ql}_@{i} = (1-beta_@{ql}_@{i})^(-1);
	@#endfor

@#endfor

% total population summing over cohorts
% including retired
Ptot= (
@#for i in 1:NLS
	@#for ql in qualif
		+ P_@{ql}_@{i}
	@#endfor
@#endfor
);

% only retired gens
Pret= (
@#for i in NT+1:NLS
	@#for ql in qualif
		+ P_@{ql}_@{i}
	@#endfor
@#endfor
);

% Health stock, aggregate
H = (
@#for i in 1:NT
	@#for ql in qualif
		+ P_@{ql}_@{i}*h_@{ql}_@{i}
	@#endfor
@#endfor
);

% Aggregate production
% kbar is tot cap
% nbar is tot labour
y = A * nbar^alp;   

% tot agg labour 
% depends on health stock
nbar = H(-1)^phi*(eta*(A_Q*L_Q)^rho + (1-eta)*L_NQ^rho)^(1/rho);

% total labour by skill
% accounts for number and prod
@#for ql in qualif
	L_@{ql} = (
	@#for i in 1:NT
		+ P_@{ql}_@{i}*a_@{ql}_@{i}                       
	@#endfor
	);
@#endfor

cCheck = Ptot - Pret - (
	@#for s in qualif
		@#for i in 1:NT
			+ P_@{s}_@{i}
		@#endfor
	@#endfor
	);

end;



%%%%%%%%%%%% STARTING VALS FOR STEADY STATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initval;
% initival provides initial guesses for solving for the actual SS
% which can differ from initial values

@#for ql in qualif

	@#for i in 1:NLS
		P_@{ql}_@{i}=1;
	@#endfor

	h_@{ql}_1=100;
	@#for i in 2:NLS
		beta_@{ql}_@{i}=.9;
		h_@{ql}_@{i}=100;

		% shocks are 0 at SS
		med_@{ql}_@{i} = 0;
		mig_@{ql}_@{i} = 0;
	@#endfor
	L_@{ql}=1.25;

@#endfor

H = 1000;
nbar = 13;
A = 1;
A_Q = 2;
xpop = 1;
pi_NQ = .7;
pi_Q = .3;
shareq = .3;
y = 6;
Ptot = 15;
cCheck = 0;


end;

% resid;

steady;
save_params_and_steady_state('./output/ss_cdcredux_pha_@{work*5+15}.txt');

%%%%% SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
perfect_foresight_setup(periods = 10);
perfect_foresight_solver(
	% linear_approximation,
    maxit = 10	
	);

%%%%% Matlab commands %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbatim;
simuls = array2table([oo_.endo_simul',oo_.exo_simul]);
simuls.Properties.VariableNames = [M_.endo_names; M_.exo_names];
match_aux = ~cellfun('isempty', regexp(simuls.Properties.VariableNames, 'AUX_', 'once'));
simuls = simuls(:, simuls.Properties.VariableNames(~match_aux));
clear AUX_* match_aux