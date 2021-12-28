%%%%% CDC clean ################################################################

% Set up vars for loops
@#define NE=1 					//%  ages of edu
@#define NT=9                   //%  working ages
@#define NTr=8                  //%  retirement ages
@#define NLS=NT+NTr             //%  total ages
@#define IR=5                   //%  retirement wage indexation ages
@#define qualif = ["Q","NQ"]    //%  skill levels

%%%%% Endogenous Vars
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

var Kmig y1 check;

var y c I g r kbar nbar Ptot Pret Tw Deped rd H R beq Tc Tk penbase retire Def D lambc Defratio Dratio Gratio Penratio Edratio Tkratio Tcratio Twratio rhop tauw tauk tauc v;

%%%%% Exogenous Vars
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
parameters alp, beta, delta, Tr, T, LS, gam, gam1, deltah, phi, eta, rho; //rep

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

%%%%% Model specification
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


% Labour supply for lower skill
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

// Stock de sante agrege
H = (
@#for i in 1:NT
	@#for ql in qualif
		+ P_@{ql}_@{i}*h_@{ql}_@{i}
	@#endfor
@#endfor
);

// Fonction de production
y = A * kbar(-1)^(1-alp) * nbar^alp;   

// Travail aggrege
nbar = H(-1)^phi*(eta*(A_Q*L_Q)^rho + (1-eta)*L_NQ^rho)^(1/rho);

@#for ql in qualif
	L_@{ql} = (
	@#for i in 1:NT
		+ P_@{ql}_@{i}*lab_@{ql}_@{i}*a_@{ql}_@{i}                       
	@#endfor
	);
@#endfor

// Productivite marginale du travail aggrege (CPO de la firme)
MPL_Q =  A * alp * kbar(-1)^(1-alp) * nbar^(alp-rho) * eta * A_Q^rho * L_Q^(rho-1);       
MPL_NQ =  A * alp * kbar(-1)^(1-alp) * nbar^(alp-rho) * (1-eta) * L_NQ^(rho-1);       

// Productivite marginale du capital aggrege (CPO de la firme)
R =  A * (1-alp) * kbar(-1)^(-alp) * nbar^alp;
// Taux de location du capital
R = r + delta;


// Salaires net des cotisation patronales
@#for i in 1:NT
	(1+tauf) * w_NQ_@{i} = MPL_NQ * a_NQ_@{i};
	@#if i <= NE
		w_Q_@{i}=0;
	@#else 
		(1+tauf) * w_Q_@{i} = MPL_Q * a_Q_@{i};
	@#endif
@#endfor


// Population totale
Ptot= (
@#for i in 1:NLS
	@#for ql in qualif
		+ P_@{ql}_@{i}
	@#endfor
@#endfor
);

// Population a la retraite
Pret= (
@#for i in NT+1:NLS
	@#for ql in qualif
		+ P_@{ql}_@{i}
	@#endfor
@#endfor
);


// Legs accidentels
Ptot*beq= (
@#for i in 1:NLS-1
	@#for ql in qualif
		+ P_@{ql}_@{i}(-1)*(1-beta_@{ql}_@{i+1}(-1))*(1+(1-tauk)*r)*k_@{ql}_@{i+1}(-1)
	@#endfor
@#endfor
);


// Pension de base
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


@#for ql in qualif

	// Indexation des retraites sur les salaires des 25 dernieres annees (5 periodes)
	penind_@{ql} =  (1/@{IR})*(
	@#for i in NT-IR+1:NT
	    + (1-tauw(@{-i}))*w_@{ql}_@{NT+1-i}(@{-i})*lab_@{ql}_@{NT+1-i}(@{-i})
	@#endfor
	);
	
	
	// Pensions de retraite
	@#for i in NT+1:NLS
		pen_@{ql}_@{i}  = 0*penbase + rhop * penind_@{ql}(@{NT-i}); 
	@#endfor

@#endfor


// Revenus de la taxation des revenus du travail
Tw = (
@#for i in 1:NT
	@#for ql in qualif
		+ P_@{ql}_@{i} * tauw * w_@{ql}_@{i} * lab_@{ql}_@{i} + P_@{ql}_@{i} * tauf * w_@{ql}_@{i} * lab_@{ql}_@{i}
	@#endfor
@#endfor
);

// Depenses de retraite
retire = (
@#for i in NT+1:NLS
	@#for ql in qualif
		+  P_@{ql}_@{i} * pen_@{ql}_@{i}
	@#endfor
@#endfor
);

// Revenus de la taxation du capital
Tk =  (
@#for i in 1:NLS-1
@#for ql in qualif
+ P_@{ql}_@{i}(-1)*tauk*r*k_@{ql}_@{i+1}(-1)
@#endfor
@#endfor
);

// Revenus de la taxation de la consommation
Tc =  (
@#for i in 1:NLS
@#for ql in qualif
+ P_@{ql}_@{i}*tauc*c_@{ql}_@{i}
@#endfor
@#endfor
);

// Depenses d'education
Deped=(
@#for i in 1:NT
@#for ql in qualif
+ P_@{ql}_@{i}*w_NQ_@{i}*e_@{ql}_@{i}                
@#endfor
@#endfor
);

// Contrainte budgetaire de l'Etat
Def= g + retire + v*Deped - (Tw+Tc+Tk);
D = (1+rd)*D(-1) + Def;
rd = r;


// Contraintes budgetaires

@#for ql in qualif
@#for i in 1:NLS
@#if i <= NT
// Contrainte budgaire des actifs
(1+tauc)*c_@{ql}_@{i} = (1 + (1-tauk)*r) * k_@{ql}_@{i}(-1) + (1-tauw)*w_@{ql}_@{i}*lab_@{ql}_@{i} - w_NQ_@{i}*(1-v)*e_@{ql}_@{i} - k_@{ql}_@{i+1} + beq;
@#else
@#if i != NLS
// Contrainte budgetaire des retraites (sauf derniere periode de vie)
    (1+tauc)*c_@{ql}_@{i} = (1 + (1-tauk)*r) * k_@{ql}_@{i}(-1) + pen_@{ql}_@{i} - k_@{ql}_@{i+1} + beq;
@#else
// Contrainte budgetaire de derniere periode de vie
    (1+tauc)*c_@{ql}_@{i} = (1 + (1-tauk)*r) * k_@{ql}_@{i}(-1) + pen_@{ql}_@{i} + beq;
@#endif
@#endif
@#endfor  
@#endfor



// Fonctions de bien-être

@#for i in 1:NLS
@#if i <= NT
welf_NQ_@{i}=log(c_NQ_@{i}) + gam1*(1-lab_NQ_@{i})^(1-gam)/(1-gam) + beta*beta_NQ_@{i+1}*welf_NQ_@{i+1}(+1);
@#if i <= NE
welf_Q_@{i}=log(c_Q_@{i}) + gam1*(1-lab_Q_@{i})^(1-gam)/(1-gam) - lambc + beta*beta_Q_@{i+1}*welf_Q_@{i+1}(+1);
@#else
welf_Q_@{i}=log(c_Q_@{i}) + gam1*(1-lab_Q_@{i})^(1-gam)/(1-gam) + beta*beta_Q_@{i+1}*welf_Q_@{i+1}(+1);
@#endif
@#else
@#if i != NLS
welf_NQ_@{i}=log(c_NQ_@{i}) + gam1/(1-gam) + beta*beta_NQ_@{i+1}*welf_NQ_@{i+1}(+1);
welf_Q_@{i}=log(c_Q_@{i}) + gam1/(1-gam) + beta*beta_Q_@{i+1}*welf_Q_@{i+1}(+1);
@#else
welf_NQ_@{NLS}=log(c_NQ_@{NLS}) + gam1/(1-gam);
welf_Q_@{NLS}=log(c_Q_@{NLS}) + gam1/(1-gam);
@#endif
@#endif
@#endfor 



// Ratio depenses et recettes sur PIB
Defratio*y=Def;
Dratio*y=D;
Penratio*y=retire;
Gratio*y=g;
Edratio*y=v*Deped;
Tkratio*y=Tk;
Tcratio*y=Tc;
Twratio*y=Tw;

//Defratio=.03;
//Edratio=.06;
v=.8;
Dratio=.3;
Gratio=.3;
Tkratio=.1;
Tcratio=.1;
Penratio=0.14;
//rhop=.5;


// Verification de la condition d'equilibre du marche des biens
y1 + Kmig = (c + I + g + Deped);
check=y + Kmig - (c + I + g + Deped);

end;



%%%%%%%%%%%%%%%%%%%%%%%%%%% STOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initval;

@#for ql in qualif

@#for i in 1:NLS
k_@{ql}_@{i} = 5;
lamda_@{ql}_@{i} =1;
c_@{ql}_@{i}=1;
P_@{ql}_@{i}=.9;
@#endfor

h_@{ql}_1=120;
@#for i in 2:NLS
beta_@{ql}_@{i}=.9;
h_@{ql}_@{i}=120;
@#endfor

@#for i in 1:NT
lab_@{ql}_@{i}=.7;
w_@{ql}_@{i}=.1;
@#endfor

@#for i in NT+1:NLS
pen_@{ql}_@{i}=.2;
@#endfor

penind_@{ql}=1;
k_@{ql}_1=0;

MPL_@{ql}=.6;
L_@{ql}=15;

@#endfor


mig_NQ_2	=	0.455175185711674	;
mig_NQ_3	=	0.348394685526794	;
mig_NQ_4	=	0.395840507124491	;
mig_NQ_5	=	0.353570128581158	;
mig_NQ_6	=	0.301831331248066	;
mig_NQ_7	=	0.313889015307074	;
mig_NQ_8	=	0.264746443288438	;
mig_NQ_9	=	0.237480784846447	;
mig_NQ_10	=	0.183183647039420	;
mig_NQ_11	=	0.139663033489348	;
mig_NQ_12	=	0.0818078074582580	;
mig_NQ_13	=	0.0399838025824241	;
mig_NQ_14	=	0.0123190699537774	;
mig_NQ_15	=	0.00227018225989047	;
mig_NQ_16	=	0.000240660078663681	;
mig_NQ_17	=	1.54557038047252e-05	;
			
mig_Q_2		=	0.455175185711674	;
mig_Q_3		=	0.348394685526794	;
mig_Q_4		=	0.395840507124491	;
mig_Q_5		=	0.353570128581158	;
mig_Q_6		=	0.301831331248066	;
mig_Q_7		=	0.313889015307074	;
mig_Q_8		=	0.264746443288438	;
mig_Q_9		=	0.237480784846447	;
mig_Q_10	=	0.183183647039420	;
mig_Q_11	=	0.139663033489348	;
mig_Q_12	=	0.0818078074582580	;
mig_Q_13	=	0.0399838025824241	;
mig_Q_14	=	0.0123190699537774	;
mig_Q_15	=	0.00227018225989047	;
mig_Q_16	=	0.000240660078663681	;
mig_Q_17	=	1.54557038047252e-05	;
			
			
med_NQ_2	=	-87.6010494357028	;
med_NQ_3	=	0.798487434274586	;
med_NQ_4	=	-1.60787772410298	;
med_NQ_5	=	-2.14057728001311	;
med_NQ_6	=	-2.36424723595204	;
med_NQ_7	=	-2.37002447387764	;
med_NQ_8	=	-2.95358946115637	;
med_NQ_9	=	-2.72140998967913	;
med_NQ_10	=	-2.24828983489933	;
med_NQ_11	=	-1.62167769058286	;
med_NQ_12	=	-1.14056323846574	;
med_NQ_13	=	-0.656289863733754	;
med_NQ_14	=	-0.374805069286264	;
med_NQ_15	=	-0.161039312071300	;
med_NQ_16	=	-0.0590031883432984	;
med_NQ_17	=	-0.0138986861007415	;
			
med_Q_2		=	-87.6010494357028	;
med_Q_3		=	0.798487434274586	;
med_Q_4		=	-1.60787772410298	;
med_Q_5		=	-2.14057728001311	;
med_Q_6		=	-2.36424723595204	;
med_Q_7		=	-2.37002447387764	;
med_Q_8		=	-2.95358946115637	;
med_Q_9		=	-2.72140998967913	;
med_Q_10	=	-2.24828983489933	;
med_Q_11	=	-1.62167769058286	;
med_Q_12	=	-1.14056323846574	;
med_Q_13	=	-0.656289863733754	;
med_Q_14	=	-0.374805069286264	;
med_Q_15	=	-0.161039312071300	;
med_Q_16	=	-0.0590031883432984	;
med_Q_17	=	-0.0138986861007415	;


c      		=	40;
y      		=	50;
r      		=	0.09;
kbar   		=	150;
nbar    	=  	30;
Tw 			=	12;
beq 		=	.01;
Kmig		= 	0;
Tk			=	1;
Tc			=	1;
Pret		=	5;
Ptot		=	15;
H 			= 	100;
retire 		= 	5;
penbase		= 	1;
rd 			=	.09;
D			=	1;
I			=	1;
Def 		=	1;
rhop		=	.5;
tauw		=	.2;



lambb   	=    50;
lambbb  	=    -50;
v 			=	.8;
tauc 		=	.1;
tauf 		=	.35;
tauk 		=	.2;
A 			= 	1;
A_Q 		= 	2;



xpop    	=    0.5893;


end;
steady;
resid;
check;