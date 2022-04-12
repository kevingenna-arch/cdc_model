%%%%% CDC clean ################################################################
% version avec P_ et beta_ endogenes, chocs produits par `cdcredux_cal_popbetas.mod`
% les chocs touchent les valeur de P_ et 

% Set up vars for loops
@#define NE=1 					//%  ages of edu
@#ifdef work
		@#define NT=work        //%  working ages
	@#else
		@#define NT=9
@#endif
@#ifdef pens
		@#define NTr=pens		//%  retirement ages
	@#else
		@#define NTr=8          
@#endif
@#define NLS=NT+NTr             //%  total ages
@#define IR=5                   //%  retirement wage indexation ages
@#define qualif = ["Q","NQ"]    //%  skill levels
@#include "../matrices_chocs.m" 	//%  external file collecting *all* shocks

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
	var L_@{ql};
    

@#endfor

var pi_Q;
var pi_NQ;
varexo shareq;
var y nbar Ptot Pret H;
var cCheck;

%%%%% EXOGENOUS VARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo A A_Q xpop;

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
phi     =   .2;         % productivity effect for health
T       =   9;          % working ages
Tr      =   8;          % retirement ages
LS      =   T+Tr;       % total ages

%%%%% MODEL SPECIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model;



pi_Q=shareq;
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
% initvals and SS values from plain file

@#if defined(work)
	load_params_and_steady_state('./output/ss_cdc_pha_work@{work*5+15}.txt');
@#elseif defined(prodind)
	load_params_and_steady_state('./output/ss_cdc_pha_pi@{prodind}.txt');
@#else
	% load_params_and_steady_state('./output/ss_cdc_pha.txt');
	initval;
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
		
		H = 3000;
		nbar = 1.5;
		A = 1;
		A_Q = 2;
		xpop = 1;
		pi_NQ = .7;
		pi_Q = .3;
		y = 1.5;
		Ptot = 3;
		cCheck = 0;
		shareq = .3;
		end;
@#endif

steady;

@#if defined(work)
	save_params_and_steady_state('./output/ss_cdcdyn_pha_chocs_@{work*5+15}.txt');
@#elseif defined(prodind)
	save_params_and_steady_state('./output/ss_cdcdyn_pha_chocs_pi@{prodind}.txt');
@#else
	save_params_and_steady_state('./output/ss_cdcdyn_pha_chocs.txt');
@#endif

%%%%%% Shocks bloc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shocks;

@#for j in 2:NLS
	@#for s in qualif
		var mig_@{s}_@{j};
		periods 240:279;
		% 10% shock
		@#if defined(mig10)
			values (s_ten_mig_@{s}_@{j});
		% old data
		@#elseif defined(old_pop)
			values (s_old_mig_@{s}_@{j});
		% baseline
		@#else
			values (s_mig_@{s}_@{j});
		@#endif

		var med_@{s}_@{j};
		periods 240:279;
		@#if defined(med30)
			% health shock to mid age
			@#if med30==1
			values (s_trent_med_@{s}_@{j});
			% health shock to early age
			@#elseif med30==2
			values (s_vingt_med_@{s}_@{j});
			@#endif
		@#elseif !defined(med30)
			%regular shocks
			values (s_med_@{s}_@{j});
		@#endif
	@#endfor
@#endfor

var xpop;
periods 240:279;
values (s_xpop);

@#ifdef q_shock
	var shareq;
	periods 260:276;
	values .5;
@#endif

@#if defined(TFP)
	var A;
	periods 240:280;

	@#if TFP == 0
	% pessimiste
	values (s_A_pess);
	
	@#elseif TFP == 1
	% optimiste
	values (s_A_opt);
	
	@#elseif TFP == 2
	% sans previsions(s_A_eff);
	values (s_A_eff);

	@#elseif TFP == 3
	values (s_A_des);

	@#else

	% scenario centrale
	values (s_A_ctr);

	@#endif

@#endif

end;

%%%%% SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
perfect_foresight_setup(periods = 500);
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