

@#define NT=10
@#define NTr=5
@#define NLS=NT+NTr
@#define qualif = ["Q","NQ"]

k = .01*ones(2,30);
lamda = .01*ones(2,30);


//1) VARIABLES************************************************************

var y y1 check c I g r MPL kbar nbar tauw Defratio Dratio Ptot Pret contr rd H R beq Tc Tk penbase retire Def D;

@#for ql in qualif

var penind_@{ql};

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

@#endfor // end for qualif

//2) EXOGENOUS VARIABLES***************************************************

varexo z tauk tauc tauf rhop ;

@#for ql in qualif
varexo h_@{ql}_1 xpop_@{ql};

@#for i in 2:NLS
varexo mig_@{ql}_@{i} med_@{ql}_@{i};     
@#endfor

@#endfor // end for qualif



//3)PARAMETERS*************************************************************

parameters alp, beta, delta, Tr, T, LS, gam, gam1, deltah, phi;

@#for i in 1:NT
parameters a_@{i};     
@#endfor

@#for i in 1:NT
a_@{i}=.7+log(@{i})/10;
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

@#for ql in qualif

P_@{ql}_1=xpop_@{ql};

@#for i in 2:NLS
P_@{ql}_@{i} =  beta_@{ql}_@{i}(-1)*P_@{ql}_@{i-1}(-1) + mig_@{ql}_@{i}/1000;
h_@{ql}_@{i} =  (1-deltah)*h_@{ql}_@{i-1}(-1) + med_@{ql}_@{i};
beta_@{ql}_@{i} = 1 - (1/h_@{ql}_@{i});
@#endfor

@#for i in 1:NLS
lamda_@{ql}_@{i} =  (1/(1+tauc))*(1/c_@{ql}_@{i});
@#endfor

@#for i in 1:NLS-1
lamda_@{ql}_@{i} = beta*beta_@{ql}_@{i+1}*lamda_@{ql}_@{i+1}(+1) * (1+(1-tauk(+1))*r(+1));
@#endfor

@#for i in 1:NT
gam1*(1-lab_@{ql}_@{i})^(-gam) = lamda_@{ql}_@{i}*(1-tauw)*w_@{ql}_@{i};
@#endfor

k_@{ql}_@{1}=0;

@#endfor


D + kbar = (
@#for i in 2:NLS
@#for ql in qualif
+ P_@{ql}_@{i-1}*k_@{ql}_@{i}                        
@#endfor
@#endfor
);


%y = c + I + g;   

y1= (c + I + g);

check=y - (c + I + g);

I = kbar - (1-delta)*kbar(-1);

c = (
@#for i in 1:NLS
@#for ql in qualif
+ P_@{ql}_@{i}*c_@{ql}_@{i}    
@#endfor
@#endfor
);

H = (
@#for i in 1:NT
@#for ql in qualif
+ P_@{ql}_@{i}*h_@{ql}_@{i} 
@#endfor
@#endfor
);


y = exp(z) * H(-1)^phi * kbar(-1)^alp * nbar^(1-alp);

nbar = (
@#for i in 1:NT
@#for ql in qualif
+ P_@{ql}_@{i}*lab_@{ql}_@{i}*a_@{i}  
@#endfor
@#endfor
);

MPL =  exp(z)* H(-1)^phi * (1-alp) * (kbar(-1)^alp) * (nbar^-alp);   
R =  exp(z) * H(-1)^phi * alp * kbar(-1)^(alp-1) * nbar^(1-alp); 
R = r + delta;

@#for ql in qualif
@#for i in 1:NT
(1+tauf)*w_@{ql}_@{i}  = MPL * a_@{i};
@#endfor
@#endfor

Ptot= (
@#for i in 1:NLS
@#for ql in qualif
+ P_@{ql}_@{i}
@#endfor
@#endfor
);

Pret= (
@#for i in NT+1:NLS
@#for ql in qualif
+ P_@{ql}_@{i}
@#endfor
@#endfor
);


Ptot*beq= (
@#for i in 1:NLS-1
@#for ql in qualif
+ P_@{ql}_@{i}(-1)*(1-beta_@{ql}_@{i+1}(-1))*(1+(1-tauk)*r)*k_@{ql}_@{i+1}(-1)
@#endfor
@#endfor
);

penbase =  (1/2)*(1/T)*(
@#for ql in qualif
@#for i in 1:NT
    + w_@{ql}_@{i}
@#endfor
@#endfor
);


@#for ql in qualif

penind_@{ql} =  (1/T)*(
@#for i in 1:NT
    + w_@{ql}_@{NT+1-i}(@{-i})
@#endfor
);


@#for i in NT+1:NLS
pen_@{ql}_@{i}  = .5 * penbase + .5 * rhop * penind_@{ql}(@{NT-i}) ; 
@#endfor

@#endfor 



contr = (
@#for i in 1:NT
@#for ql in qualif
+ P_@{ql}_@{i} * tauw * w_@{ql}_@{i} * lab_@{ql}_@{i} + P_@{ql}_@{i} * tauf * w_@{ql}_@{i} * lab_@{ql}_@{i}
@#endfor
@#endfor
);

retire = (
@#for i in NT+1:NLS
@#for ql in qualif
+  P_@{ql}_@{i} * pen_@{ql}_@{i}
@#endfor
@#endfor
);

//retire=contr;

Tk =  (
@#for i in 1:NLS-1
@#for ql in qualif
+ P_@{ql}_@{i}(-1)*tauk*r*k_@{ql}_@{i+1}(-1)              
@#endfor
@#endfor
);

Tc =  (
@#for i in 1:NLS
@#for ql in qualif
+ P_@{ql}_@{i}*tauc*c_@{ql}_@{i}                        
@#endfor
@#endfor
);


// Government budget constraint

D=.2*y;
g=.1*y;
Def= g + retire - (contr+Tc+Tk);
D = (1+rd)*D(-1) + Def;
rd = r;


Defratio=Def/y;
Dratio=D/y;


@#for ql in qualif

@#for i in 1:NLS

@#if i <= NT

(1+tauc)*c_@{ql}_@{i} = (1 + (1-tauk)*r) * k_@{ql}_@{i}(-1) + (1-tauw)*w_@{ql}_@{i}*lab_@{ql}_@{i} - k_@{ql}_@{i+1} + beq;  //worker consumption

@#else

@#if i != NLS

    (1+tauc)*c_@{ql}_@{i} = (1 + (1-tauk)*r) * k_@{ql}_@{i}(-1) + pen_@{ql}_@{i} - k_@{ql}_@{i+1} + beq;        //retiree consumption

@#else

    (1+tauc)*c_@{ql}_@{i} = (1 + (1-tauk)*r) * k_@{ql}_@{i}(-1) + pen_@{ql}_@{i} + beq;                   //retiree last period consumption

@#endif

@#endif

@#endfor 




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

xpop_@{ql}=1;

k_@{ql}_1=0;

@#endfor


c      	=	 40;
y      	=	 50;
r      	=	 0.09;
MPL      	=	 1.6;
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
%check;
%resid;

%shocks;
%var tauw;
%periods 10:50;
%values 0;
%end;

endval;
%z   =   .01;
%xpop_Q=1.01;
med_NQ_12=50;
%med_50=50;
%mig_Q_5=10;
%tauk=.2;
%tauw = .15;
%tauf = .55;
%med_20=-50;
end;
steady;
%check;
%resid;

perfect_foresight_setup(periods=40);

perfect_foresight_solver;

run plotsim.m

