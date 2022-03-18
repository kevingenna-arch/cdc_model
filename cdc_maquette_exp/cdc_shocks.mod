%%%%% CDC clean ################################################################
%% VERSION NETTOYEE DE LA MAQUETTE COMPLETE CDC


% Set up vars for loops
@#define NE=1 					//%  ages of edu
@#define NT=9                   //%  working ages
@#define NTr=8                  //%  retirement ages
@#define NLS=NT+NTr             //%  total ages
@#define IR=5                   //%  retirement wage indexation ages
@#define qualif = ["Q","NQ"]    //%  skill levels
@#include "../matrices_chocs.m" 	//%  external file collecting *all* shocks


%%%%% ENDOGENOUS VARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% endogenous vars from agents' maxing programs, by skill level
@#for ql in qualif
	% capital stock by age
	@#for i in 1:NLS
		var k_@{ql}_@{i};          
	@#endfor
	% consumption by age
	@#for i in 1:NLS
		var c_@{ql}_@{i};        
	@#endfor
	% labour supply by working age
	@#for i in 1:NT
		var lab_@{ql}_@{i}; 
	@#endfor
	% labour wage by working age
	@#for i in 1:NT
		var w_@{ql}_@{i};     
	@#endfor
	% survival discount factor between periods-- CHECK
	@#for i in 2:NLS
		var beta_@{ql}_@{i};
	@#endfor
	% health stock -- CHECK
	@#for i in 2:NLS
		var h_@{ql}_@{i};
	@#endfor
	% per-period lagrange multiplier - KIM BCs are not consolidated
	@#for i in 1:NLS
		var lamda_@{ql}_@{i};   
	@#endfor
	% pension by retirement age
	@#for i in NT+1:NLS
		var pen_@{ql}_@{i};     
	@#endfor
	% cohort populations
	@#for i in 1:NLS
		var P_@{ql}_@{i} ;
	@#endfor
	% welfare function value by age
	@#for i in 1:NLS
		var welf_@{ql}_@{i};
	@#endfor
	% aggregates by skill class:
	% -	L: aggregate labour
	% - MPL: marginal produc of labour
	% -	penind: indexed pension wage
	% - pi: share of first gen going for edu -- mirroring 
	var L_@{ql} MPL_@{ql} penind_@{ql} pi_@{ql};

@#endfor

var Kmig y1 cCheck;

var y c I g r kbar nbar Ptot Pret Tw Deped rd H R beq Tc Tk penbase retire Def D lambc Defratio Dratio Gratio Penratio Edratio Tkratio Tcratio Twratio rhop tauw tauk tauc v;

%%%%% EXOGENOUS VARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo A tauf xpop lambb lambbb A_Q;

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
phi     =   0.1;          %
T       =   9;          % working ages
Tr      =   8;          % retirement ages
LS      =   T+Tr;       % total ages

%%%%% MODEL SPECIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model;

% education choice on recursive discounted welfare
0=welf_Q_1-welf_NQ_1;

pi_NQ=(lambb-lambc)/(lambb-lambbb);
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

	% FOCs in eq-by-eq form, explicit multipliers
	@#for i in 1:NLS
		% marginal utility of consumption
		% intratemporal allocation
		lamda_@{ql}_@{i} =  (1/(1+tauc))*(1/c_@{ql}_@{i});
	@#endfor

	@#for i in 1:NLS-1
		% Euler eq
		% intertemporal allocation
		lamda_@{ql}_@{i} = beta*beta_@{ql}_@{i+1}*lamda_@{ql}_@{i+1}(+1) * (1+(1-tauk(+1))*r(+1));
	@#endfor

	% start of period K stock
	% shouldn't this be in the initval bloc??
	k_@{ql}_@{1}=0;

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


% Labour supply for lower skill
% individual levels
@#for i in 1:NT
	% marginal disutility of labour equates net wage
	gam1*(1-lab_NQ_@{i})^(-gam) = lamda_NQ_@{i} * (1-tauw)*w_NQ_@{i} +
	% future discounted sum of wages due to indexation upping supply
	@#for j in 1:NTr
		+ lamda_NQ_@{NLS-j+1}(+@{NLS-i-j+1})*beta^(@{NLS-j+1-i})*rhop(+@{NLS-i-j+1})*(1-tauw)*w_NQ_@{i}*(1/5)
	@#endfor
;
@#endfor

% Labour supply for higher skill
% individual levels
@#for i in 1:NT
	% condition for education
	@#if i <= NE
		lab_Q_@{i}=0;
	@#else 
	% marginal disutility of labour equates net wage
		gam1*(1-lab_Q_@{i})^(-gam) = lamda_Q_@{i} * (1-tauw)*w_Q_@{i} +
	@#for j in 1:NTr
	% future discounted sum of wages due to indexation upping supply
		+ lamda_Q_@{NLS-j+1}(+@{NLS-i-j+1})*beta^(@{NLS-i-j+1})*rhop(+@{NLS-i-j+1})*(1-tauw)*w_Q_@{i}*(1/5)
	@#endfor
	;
	@#endif
@#endfor


% public debt and capital stock
% from individual accumulation
D + kbar = (
@#for i in 2:NLS
	@#for ql in qualif
		+ P_@{ql}_@{i-1}*k_@{ql}_@{i}
	@#endfor
@#endfor
);

% net capital from migration
Kmig=(
@#for i in 2:NLS
	@#for ql in qualif
		+ mig_@{ql}_@{i}*(1+(1-tauk)*r)*k_@{ql}_@{i}(-1)
	@#endfor
@#endfor
);

% investment
I = kbar - (1-delta)*kbar(-1);

% total consumption
c = (
@#for i in 1:NLS
	@#for ql in qualif
		+ P_@{ql}_@{i}*c_@{ql}_@{i}                            
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
y = A * kbar(-1)^(1-alp) * nbar^alp;   

% tot agg labour 
% depends on health stock
nbar = H(-1)^phi*(eta*(A_Q*L_Q)^rho + (1-eta)*L_NQ^rho)^(1/rho);

% total labour by skill
% accounts for number and prod
@#for ql in qualif
	L_@{ql} = (
	@#for i in 1:NT
		+ P_@{ql}_@{i}*lab_@{ql}_@{i}*a_@{ql}_@{i}                       
	@#endfor
	);
@#endfor

% MPL by skill from agg prod
MPL_Q =  A * alp * kbar(-1)^(1-alp) * nbar^(alp-rho) * eta * A_Q^rho * L_Q^(rho-1);       
MPL_NQ =  A * alp * kbar(-1)^(1-alp) * nbar^(alp-rho) * (1-eta) * L_NQ^(rho-1);       

% MPK 
R =  A * (1-alp) * kbar(-1)^(-alp) * nbar^alp;
% gross interest rate
R = r + delta;

% wages without labour taxes for firms
@#for i in 1:NT
	(1+tauf) * w_NQ_@{i} = MPL_NQ * a_NQ_@{i};
	@#if i <= NE
		w_Q_@{i}=0;
	@#else 
		(1+tauf) * w_Q_@{i} = MPL_Q * a_Q_@{i};
	@#endif
@#endfor


% bequests before last period (certain end of life)
Ptot*beq= (
@#for i in 1:NLS-1
	@#for ql in qualif
		+ P_@{ql}_@{i}(-1)*(1-beta_@{ql}_@{i+1}(-1))*(1+(1-tauk)*r)*k_@{ql}_@{i+1}(-1)
	@#endfor
@#endfor
);


% Retirement -- regular
% this at individual level
penbase =  ((@{NT-NE}/(@{NT-NE}+@{NT}))*(
@#for i in NE+1:NT
    + (1-tauw)*w_Q_@{i}*lab_Q_@{i}
@#endfor
))
+
((@{NT}/(@{NT-NE}+@{NT}))*(
@#for i in 1:NT
    + (1-tauw)*w_NQ_@{i}*lab_NQ_@{i}
@#endfor
));

% retirement indexation
@#for ql in qualif
	% indexing on last 5(?) periods
	penind_@{ql} =  (1/@{IR})*(
	@#for i in NT-IR+1:NT
	    + (1-tauw(@{-i}))*w_@{ql}_@{NT+1-i}(@{-i})*lab_@{ql}_@{NT+1-i}(@{-i})
	@#endfor
	);
	
	
	% weighting both pensions schemes
	@#for i in NT+1:NLS
		% rhop gives the balance, by now tunred off the base pens
		pen_@{ql}_@{i}  = 0*(1-rhop)*penbase + rhop * penind_@{ql}(@{NT-i}); 
	@#endfor

@#endfor


% total labour tax 
Tw = (
@#for i in 1:NT
	@#for ql in qualif
		+ P_@{ql}_@{i} * tauw * w_@{ql}_@{i} * lab_@{ql}_@{i} + P_@{ql}_@{i} * tauf * w_@{ql}_@{i} * lab_@{ql}_@{i}
	@#endfor
@#endfor
);

% total capital tax
Tk =  (
@#for i in 1:NLS-1
	@#for ql in qualif
		+ P_@{ql}_@{i}(-1)*tauk*r*k_@{ql}_@{i+1}(-1)
	@#endfor
@#endfor
);

% total consumption tax
Tc =  (
@#for i in 1:NLS
	@#for ql in qualif
		+ P_@{ql}_@{i}*tauc*c_@{ql}_@{i}
	@#endfor
@#endfor
);

% public spending for retirement wages
% alredy combines indexation and base
retire = (
@#for i in NT+1:NLS
	@#for ql in qualif
		+  P_@{ql}_@{i} * pen_@{ql}_@{i}
	@#endfor
@#endfor
);

% educ spending by gvt
Deped=(
@#for i in 1:NT
	@#for ql in qualif
		+ P_@{ql}_@{i}*w_NQ_@{i}*e_@{ql}_@{i}                
	@#endfor
@#endfor
);

% gvt budget constraint
Def= g + retire + v*Deped - (Tw+Tc+Tk);
% debt dynamics
D = (1+rd)*D(-1) + Def;
% debt interests
% same as any other K
rd = r;


% sequences of budget constraints
@#for ql in qualif
	@#for i in 1:NLS
		@#if i <= NT
			% workers
			(1+tauc)*c_@{ql}_@{i} = (1 + (1-tauk)*r) * k_@{ql}_@{i}(-1) + (1-tauw)*w_@{ql}_@{i}*lab_@{ql}_@{i} - w_NQ_@{i}*(1-v)*e_@{ql}_@{i} - k_@{ql}_@{i+1} + beq;
		@#else
		@#if i != NLS
			% retired before last period
	    (	1+tauc)*c_@{ql}_@{i} = (1 + (1-tauk)*r) * k_@{ql}_@{i}(-1) + pen_@{ql}_@{i} - k_@{ql}_@{i+1} + beq;
		@#else
			% budget exhaustion condition -- last period
	    	(1+tauc)*c_@{ql}_@{i} = (1 + (1-tauk)*r) * k_@{ql}_@{i}(-1) + pen_@{ql}_@{i} + beq;
		@#endif
		@#endif
	@#endfor  
@#endfor


% Welfare functiton: discounted sum of consumption/labour streams
% just NPV of U function, no health so far
@#for i in 1:NLS
	@#if i <= NT
		% WF for unskilled
		welf_NQ_@{i}=log(c_NQ_@{i}) + gam1*(1-lab_NQ_@{i})^(1-gam)/(1-gam) + beta*beta_NQ_@{i+1}*welf_NQ_@{i+1}(+1);
		% WF for skilled
		@#if i <= NE
		% accounting for educ period
		welf_Q_@{i}=log(c_Q_@{i}) + gam1*(1-lab_Q_@{i})^(1-gam)/(1-gam) - lambc + beta*beta_Q_@{i+1}*welf_Q_@{i+1}(+1);
		@#else
		% working age periods WF 
		welf_Q_@{i}=log(c_Q_@{i}) + gam1*(1-lab_Q_@{i})^(1-gam)/(1-gam) + beta*beta_Q_@{i+1}*welf_Q_@{i+1}(+1);
	@#endif
	@#else
	@#if i != NLS
		% retirement WFs
		welf_NQ_@{i}=log(c_NQ_@{i}) + gam1/(1-gam) + beta*beta_NQ_@{i+1}*welf_NQ_@{i+1}(+1);
		welf_Q_@{i}=log(c_Q_@{i}) + gam1/(1-gam) + beta*beta_Q_@{i+1}*welf_Q_@{i+1}(+1);
	@#else
		% last period WF
		welf_NQ_@{NLS}=log(c_NQ_@{NLS}) + gam1/(1-gam);
		welf_Q_@{NLS}=log(c_Q_@{NLS}) + gam1/(1-gam);
	@#endif
	@#endif
@#endfor 



% ratios:
% deficit over GDP
Defratio*y=Def;
% debt over GDP
Dratio*y=D;
% retirement spending over GDP
Penratio*y=retire;
% g is just gvt wasteful consumption, not having any effect here
% discretional spendig/consumption over GDP
Gratio*y=g;
% eduction spending over GDP
% weighted by
Edratio*y=v*Deped;
% capital income tax revenue over GDP
Tkratio*y=Tk;
% consumption tax revenue over GDP
Tcratio*y=Tc;
% labour income tax revenue over GDP
Twratio*y=Tw;

% setting some ratios !!! within model block? !!!
% v is participation to educ expenses by individual
v=.8;
% target 30% debt to GDP
Dratio=.3;
% discretional spending to 30%
Gratio=.3;
% K tax revenue to 10%
Tkratio=.1;
% C tax revenue to 10%
Tcratio=.1;
% pension exp ratio to 14%
Penratio=0.14;

% consistency checks
% y1 mops up any residual on goods' mkt
y1 + Kmig = (c + I + g + Deped);
% check ought be 0
cCheck=y + Kmig - (c + I + g + Deped);

end;



%%%%%%%%%%%% STARTING VALS FOR STEADY STATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load_params_and_steady_state('./output/ss_cdclean.txt');
steady;
%%%%%%% Save SS values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_params_and_steady_state('./output/ss_cdclean_endo.txt');

%%%%% SHOCKS BLOCK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shocks;
% in simulated series timing is shifted by one period
% if a shock hits from t=10 to t=15, its position in the 
% series is t=11 to t=16, retrieved by 11:16 in MTLB
@#for j in 2:NLS
	@#for s in qualif
		var mig_@{s}_@{j};
		% periods 240:279;
		periods 1:40;
		values (s_full_mig_@{s}_@{j});

		var med_@{s}_@{j};
		% periods 240:279;
		periods 1:40;
		values (s_full_med_@{s}_@{j});
	@#endfor
@#endfor

var xpop;
% periods 240:279;
periods 1:40;
values (s_xpop);

var A;
% periods 240:280;
periods 1:41;
values (s_A_ctr);

end;

% resid;

% check;

perfect_foresight_setup(periods = 100);
perfect_foresight_solver(
	% linear_approximation,
    maxit = 2
	);

verbatim;
simuls = array2table([oo_.endo_simul',oo_.exo_simul]);
simuls.Properties.VariableNames = [M_.endo_names; M_.exo_names];
match_aux = ~cellfun('isempty', regexp(simuls.Properties.VariableNames, 'AUX_', 'once'));
simuls = simuls(:, simuls.Properties.VariableNames(~match_aux));
clear AUX_* match_aux