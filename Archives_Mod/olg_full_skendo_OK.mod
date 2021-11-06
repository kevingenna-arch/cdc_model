
@#define NE=2
@#define NT=10
@#define NTr=5
@#define NLS=NE+NT+NTr
@#define qualif = ["Q","NQ"]

k = .01*ones(2,30);
lamda = .01*ones(2,30);


//1) VARIABLES************************************************************

var c I g r MPL kbar nbar tauw Defratio Dratio Penratio Edratio Ptot Pret contr Deped rd H R beq Tc Tk penbase retire Def D lambc;

@#for ql in qualif

var penind_@{ql} pi_@{ql};

@#for i in 2:NLS
var beta_@{ql}_@{i};
@#endfor

@#for i in 2:NLS
var h_@{ql}_@{i};
@#endfor

@#for i in 1:NT
var lab_@{ql}_@{i};     
@#endfor

@#for i in 1:NT
var w_@{ql}_@{i};   
@#endfor

@#for i in 1:NLS
var k_@{ql}_@{i};         
@#endfor

@#for i in 1:NLS
var c_@{ql}_@{i};       
@#endfor

@#for i in 1:NLS
var lamda_@{ql}_@{i};     
@#endfor

@#for i in NT+1:NLS
var pen_@{ql}_@{i};  
@#endfor

@#for i in 1:NLS
var P_@{ql}_@{i} ;
@#endfor

@#for i in 1:NLS
var welf_@{ql}_@{i};
@#endfor

@#endfor

var y y1 check;

//2) EXOGENOUS VARIABLES***************************************************

varexo z tauk tauc tauf rhop xpop lambb lambbb v;

@#for ql in qualif
varexo h_@{ql}_1;

@#for i in 2:NLS
varexo mig_@{ql}_@{i} med_@{ql}_@{i};   
@#endfor

@#endfor



//3)PARAMETERS*************************************************************

parameters alp, beta, delta, Tr, T, LS, gam, gam1, deltah, phi;

@#for ql in qualif
@#for i in 1:NT
parameters a_@{ql}_@{i};    
@#endfor
@#endfor

@#for i in 1:NT
a_NQ_@{i}=.7+log(@{i})/10;
@#if i <= NE
a_Q_@{i}=.7+log(@{i})/10;
@#else 
a_Q_@{i}=1+log(@{i})/10;
@#endif
@#endfor

gam= 1.5;
gam1=.07;
alp=0.3;  
beta=0.97;  
delta=0.02;  
deltah=0.02;  
phi=.02;
T=10;       
Tr=5;  
LS=T+Tr; 


//4) MODEL EQUATION*******************************************************
model;

//Choix d'éducation
%lambc=.5;
0=welf_Q_1-welf_NQ_1;

pi_NQ=min(1,max(0,(lambb-lambc)/(lambb-lambbb)));
%pi_NQ=min(1,max(0,(lambc-lambbb)/(lambb-lambbb)));
pi_Q=1-pi_NQ;


@#for ql in qualif


// Bloc démographique
P_@{ql}_1=pi_@{ql}*xpop;

@#for i in 2:NLS
P_@{ql}_@{i} =  beta_@{ql}_@{i}(-1)*P_@{ql}_@{i-1}(-1) + mig_@{ql}_@{i}/1000;
h_@{ql}_@{i} =  (1-deltah)*h_@{ql}_@{i-1}(-1) + med_@{ql}_@{i};
beta_@{ql}_@{i} = 1 - (1/h_@{ql}_@{i});
@#endfor



// Mutliplicateur de lagrange, Utilité marginale de la consommation (CPO des ménages)
@#for i in 1:NLS
lamda_@{ql}_@{i} =  (1/(1+tauc))*(1/c_@{ql}_@{i});
@#endfor

// Equation d'Euler (CPO des ménages)
@#for i in 1:NLS-1
lamda_@{ql}_@{i} = beta*beta_@{ql}_@{i+1}*lamda_@{ql}_@{i+1}(+1) * (1+(1-tauk(+1))*r(+1));
@#endfor

// Stock de capital inital ???
k_@{ql}_@{1}=0;

@#endfor


// Offre de travail (CPO des ménages)
@#for i in 1:NT
gam1*(1-lab_NQ_@{i})^(-gam) = lamda_NQ_@{i} * (1-tauw)*w_NQ_@{i};
@#endfor

// Offre de travail (CPO des ménages)
@#for i in 1:NT
@#if i <= NE
lab_Q_@{i}=0;
@#else 
gam1*(1-lab_Q_@{i})^(-gam) = lamda_Q_@{i} * (1-tauw)*w_Q_@{i};
@#endif
@#endfor





// Equilibre du marché financier
D + kbar = (
@#for i in 2:NLS
@#for ql in qualif
+ P_@{ql}_@{i-1}*k_@{ql}_@{i}
@#endfor
@#endfor
);


// Vérification
y1= (c + I + g + Deped); //Condition d'équilibre du marché des biens
check=y - (c + I + g + Deped); //Condition d'équilibre du marché des biens

// Capital importé???


// Investissement
I = kbar - (1-delta)*kbar(-1);

// Consommation aggrégée
c = (
@#for i in 1:NLS
@#for ql in qualif
+ P_@{ql}_@{i}*c_@{ql}_@{i}                            
@#endfor
@#endfor
);

// Stock de santé aggrégé
H = (
@#for i in 1:NT
@#for ql in qualif
+ P_@{ql}_@{i}*h_@{ql}_@{i}                          
@#endfor
@#endfor
);

// Fonction de production
y = exp(z) * H(-1)^phi * kbar(-1)^alp * nbar^(1-alp);   

// Travail aggrégé
nbar = (
@#for i in 1:NT
@#for ql in qualif
+ P_@{ql}_@{i}*lab_@{ql}_@{i}*a_@{ql}_@{i}                       
@#endfor
@#endfor
);

// Productivité marginale du travail aggrégé (CPO de la firme)
MPL =  exp(z)* H(-1)^phi * (1-alp) * (kbar(-1)^alp) * (nbar^-alp);       
 
// Productivité marginale du capital aggrégé (CPO de la firme)
R =  exp(z) * H(-1)^phi * alp * kbar(-1)^(alp-1) * nbar^(1-alp);
// Taux de location du capital
R = r + delta;


// Salaires (CPO de la firme)
@#for i in 1:NT
(1+tauf) * w_NQ_@{i} = MPL * a_NQ_@{i};
@#if i <= NE
w_Q_@{i}=0;
@#else 
(1+tauf) * w_Q_@{i} = MPL * a_Q_@{i};
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

// Population à la retraite
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
penbase =  (1/2)*(1/@{NT-NE})*(
@#for ql in qualif
@#for i in NE+1:NT
    + w_@{ql}_@{i}
@#endfor
@#endfor
);


@#for ql in qualif

// Indéxation des retraites sur les salaires passés
penind_@{ql} =  (1/(@{NT-NE}))*(
@#for i in NE+1:NT
    + w_@{ql}_@{NT+1-i}(@{-i})
@#endfor
);

// Pensions de retraite
@#for i in NT+1:NLS
pen_@{ql}_@{i}  = .5 * penbase + .5 * rhop * penind_@{ql}(@{NT-i}) ;      
@#endfor

@#endfor // end for qualif


// Revenus de la taxation du travail
contr = (
@#for i in 1:NT
@#for ql in qualif
+ P_@{ql}_@{i} * tauw * w_@{ql}_@{i} * lab_@{ql}_@{i} + P_@{ql}_@{i} * tauf * w_@{ql}_@{i} * lab_@{ql}_@{i}
@#endfor
@#endfor
);

// Dépenses de retraite
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

// Dépenses d'éducation
Deped=(
@#for i in 1:NT
@#for ql in qualif
+ P_@{ql}_@{i}*w_NQ_@{i}*lambc                 
@#endfor
@#endfor
);

// Contrainte budgétaire de l'état
Def= g + retire + v*Deped - (contr+Tc+Tk);
D = (1+rd)*D(-1) + Def;
rd = r;

// Ratio
Defratio=Def/y;
Dratio=D/y;
Penratio=retire/y;
Edratio=v*Deped/y;

D=.2*y;
g=.1*y;
//Penratio=0.14*y;


// Contraintes budgétaires

@#for ql in qualif

@#for i in 1:NLS

@#if i <= NT

(1+tauc)*c_@{ql}_@{i} = (1 + (1-tauk)*r) * k_@{ql}_@{i}(-1) + (1-tauw)*w_@{ql}_@{i}*lab_@{ql}_@{i} - w_NQ_@{i}*lambc*(1-v) - k_@{ql}_@{i+1} + beq;  //worker consumption

@#else

@#if i != NLS

    (1+tauc)*c_@{ql}_@{i} = (1 + (1-tauk)*r) * k_@{ql}_@{i}(-1) + pen_@{ql}_@{i} - k_@{ql}_@{i+1} + beq;        //retiree consumption

@#else

    (1+tauc)*c_@{ql}_@{i} = (1 + (1-tauk)*r) * k_@{ql}_@{i}(-1) + pen_@{ql}_@{i} + beq;                   //retiree last period consumption

@#endif

@#endif

@#endfor  



// Fonctions de bien-être

@#for i in 1:NLS

@#if i <= NT

welf_@{ql}_@{i}=log(c_@{ql}_@{i}) + gam1*(1-lab_@{ql}_@{i})^(1-gam)/(1-gam) + beta*beta_@{ql}_@{i+1}*welf_@{ql}_@{i+1}(+1);

@#else

@#if i != NLS

welf_@{ql}_@{i}=log(c_@{ql}_@{i}) + gam1/(1-gam) + beta*beta_@{ql}_@{i+1}*welf_@{ql}_@{i+1}(+1);

@#else

welf_@{ql}_@{NLS}=log(c_@{ql}_@{NLS}) + gam1/(1-gam);

@#endif

@#endif

@#endfor  


@#endfor 



end;  




//5) INITIALISATION WITH STEADY STATE VALUES*******************************

initval;

@#for i in 2:NLS
mig_@{ql}_@{i}=0;
med_@{ql}_@{i}=0;
@#endfor

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

@#endfor

v=.8;
lambb   =    10;
lambbb  =    -5;
xpop    =    1;
c      	=	 40;
y      	=	 50;
r      	=	 0.09;
MPL     =	 1.6;
kbar   	=	 150;
nbar    =   30;
contr=12;
beq=.01;
Tk=1;
Tc=1;
Pret=5;
Ptot=15;
H=100;
retire=5;
penbase=1;
rhop=.5;
rd=.09;
D=1;
I=1;
Def=1;
tauc=.1;
tauf=.35;
tauk=.2;
tauw=.3;
z = 0;

end; 
steady;
%resid;

