%%%%% CDC clean ################################################################

% Set up vars for loops
@#define NE=1 					//%  ages of edu
@#define NT=9                   //%  working ages
@#define NTr=8                  //%  retirement ages
@#define NLS=NT+NTr             //%  total ages
@#define IR=5                   //%  retirement wage indexation ages
@#define qualif = ["Q","NQ"]    //%  skill levels
@#include "matrices_chocs.m" 	//

%%%%% ENDOGENOUS VARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% endogenous vars from agents' maxing programs, by skill level
@#for ql in qualif
	% health stock -- CHECK
	@#for j in 2:NLS
		var beta_@{ql}_@{j};
	@#endfor
	% aggregates by skill class:
	% -	L: aggregate labour
	% - pi: share of first gen going for edu -- mirroring 
	var L_@{ql} pi_@{ql};
	var P_@{ql}_1;


	@#for i in 2:NLS
		% shocks for migration and health stock
		var mig_@{ql}_@{i} med_@{ql}_@{i};
	@#endfor

@#endfor

var y nbar Ptot Pret H;
var cCheck;
% var junk;

%%%%% EXOGENOUS VARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo A A_Q xpop;

@#for ql in qualif
	% @#for i in 2:NLS
	% 	% shocks for migration and health stock
	% 	varexo mig_@{ql}_@{i} med_@{ql}_@{i};
	% @#endfor

	@#for i in 1:NLS
		varexo h_@{ql}_@{i};
	@#endfor
	% cohort populations
	@#for i in 2:NLS
		varexo P_@{ql}_@{i} ;
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
	% productivity by age
	a_NQ_@{i}=.7+log(@{i})/10;
	% switcher for education
	e_NQ_@{i}=0;
	@#if i <= NE
		% values for prod and switcher for skills
		a_Q_@{i}=.7+log(@{i})/10;
		e_Q_@{i}=1;
	@#else 
		a_Q_@{i}=1+log(@{i})/10;
		e_Q_@{i}=0;
	@#endif
@#endfor

gam     =   1.5;        %
gam1    =   .07;        %
alp     =   .7;         % 
eta     =   .5;         %
rho     =   .5;         % 
beta    =   0.97;       % discount factor
delta   =   0.02;       % physical capital depreciation
deltah  =   0.02;       % health depreciation
phi     =   0;          %
T       =   9;          % working ages
Tr      =   8;          % retirement ages
LS      =   T+Tr;       % total ages

%%%%% MODEL SPECIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model;
% junk=0.9*junk(+1);

pi_NQ=.7;
pi_Q=1-pi_NQ;

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
		h_@{ql}_@{i}*beta_@{ql}_@{i} = h_@{ql}_@{i} - 1;
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
% initvals and SS values from plain file
load_params_and_steady_state('./output/ss_cdcsimple_ph_STABLE.txt');


steady;
save_params_and_steady_state('./output/ss_cdcsimple_ph_EXO.txt');

%%%%%% Shocks bloc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shocks;

@#for j in 2:NLS
	@#for s in qualif
		var P_@{s}_@{j};
		periods 240:279;

		@#if old == 1
		values (s_old_P_@{s}_@{j});
		@#endif

		@#if old == 0
		values (s_P_@{s}_@{j});
		@#endif

		var h_@{s}_@{j};
		periods 240:279;
		values (100*s_h_@{s}_@{j});
	@#endfor
@#endfor

var xpop;
periods 240:279;
values (s_xpop);
end;

perfect_foresight_setup(periods = 500);
perfect_foresight_solver(
    maxit = 10	
	);

verbatim;
simuls = array2table([oo_.endo_simul',oo_.exo_simul]);
simuls.Properties.VariableNames = [M_.endo_names; M_.exo_names];
match_aux = ~cellfun('isempty', regexp(simuls.Properties.VariableNames, 'AUX_', 'once'));
simuls = simuls(:, simuls.Properties.VariableNames(~match_aux));
clear AUX_* match_aux

%% storing away shocks
% select only matching vars
matched_mig = ~cellfun('isempty', regexp(simuls.Properties.VariableNames, 'mig_', 'once'));
mig_shocks = simuls(:, simuls.Properties.VariableNames(matched_mig));
mig_flipped = rows2vars(mig_shocks);
mig_flipped.Properties.RowNames = table2array(mig_flipped(:, 1));
writetable(mig_flipped(:, (240:279)+2), './output/mig_shocks_redux.xlsx', 'WriteVariableNames',false, 'WriteRowNames', true);

matched_med = ~cellfun('isempty', regexp(simuls.Properties.VariableNames, 'med_', 'once'));
med_shocks = simuls(:, simuls.Properties.VariableNames(matched_med));
med_flipped = rows2vars(med_shocks);
med_flipped.Properties.RowNames = table2array(med_flipped(:, 1));
writetable(med_flipped(:, (240:279)+2), './output/med_shocks_redux.xlsx', 'WriteVariableNames',false, 'WriteRowNames', true);
