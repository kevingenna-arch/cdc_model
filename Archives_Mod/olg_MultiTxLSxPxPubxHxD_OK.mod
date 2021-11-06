

@#define NT=47
@#define NTr=38
@#define NLS=NT+NTr

k = .01*ones(2,30);
lamda = .01*ones(2,30);


//1) VARIABLES************************************************************

var  y c I g r MPL kbar nbar tauw Defratio Dratio Ptot Pret contr rd H R beq Tc Tk retire penind penbase Def D;  //note I have not include eps /*rhop*/

@#for i in 2:NLS
var beta_@{i};
@#endfor

@#for i in 2:NLS
var h_@{i};
@#endfor

@#for i in 1:NT
var lab_@{i};      
@#endfor

@#for i in 1:NT
var w_@{i};      
@#endfor

@#for i in 1:NLS
var ki@{i};         
@#endfor


@#for i in 1:NLS
var c_@{i};          
@#endfor

@#for i in 1:NLS
var lamda@{i};    
@#endfor

@#for i in NT+1:NLS
var pen_@{i}; 
@#endfor

@#for i in 1:NLS
var P_@{i} ;
@#endfor

@#for i in 1:NLS
var welf_@{i};
@#endfor

var check;

//2) EXOGENOUS VARIABLES***************************************************

varexo z tauk tauc tauf xpop h_1 rhop;

@#for i in 2:NLS
varexo mig_@{i} med_@{i};   
@#endfor

//predetermined_variables ki;


//3)PARAMETERS*************************************************************

parameters alp, beta, delta, Tr, T, LS, gam, gam1, deltah, phi; 

@#for i in 1:NT
parameters a_@{i};  
@#endfor

@#for i in 1:NT
a_@{i}=.7+log(@{i})/10;
@#endfor



gam= 1.5; // or value of 2
gam1=.07;
alp=0.3;    
beta=0.97; 
delta=0.02;    
deltah=0.02;   
phi=.02;
T=47;       
Tr=38;   
LS=T+Tr;     


//4) MODEL EQUATION*******************************************************
model;

//+++++++++++++++++++++++++Euler equation for worker and retiree+++++++++

P_1=xpop;

@#for i in 2:NLS
P_@{i} =  beta_@{i}(-1)*P_@{i-1}(-1) + mig_@{i}/1000;
h_@{i} =  (1-deltah)*h_@{i-1}(-1) + med_@{i};
beta_@{i} = 1 - (1/h_@{i});
@#endfor


@#for i in 1:NLS
lamda@{i} =  (1/(1+tauc))*(1/c_@{i});
@#endfor

@#for i in 1:NLS-1
lamda@{i} = beta*beta_@{i+1}*lamda@{i+1}(+1) * (1+(1-tauk(+1))*r(+1));
@#endfor


@#for i in 1:NT
gam1*(1-lab_@{i})^(-gam) = lamda@{i}*(1-tauw)*w_@{i} ;
@#endfor

/*
@#for i in 1:NT
lab_@{i} = 1;
@#endfor
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ki@{1}=0; 
D + kbar = (
@#for i in 2:NLS
+ P_@{i-1}(-1)*ki@{i}                           
@#endfor
);


%y = c + I + g;       
check=y - (c + I + g);

I = kbar - (1-delta)*kbar(-1);

c = (
@#for i in 1:NLS
+ P_@{i}*c_@{i}                   
@#endfor
);

H = (
@#for i in 1:NT
+ P_@{i}*h_@{i}                        
@#endfor
);


y = exp(z) * H(-1)^phi * kbar(-1)^alp * nbar^(1-alp);   

nbar = (
@#for i in 1:NT
+ P_@{i}*lab_@{i}*a_@{i}                           
@#endfor
);

MPL =  exp(z)* H(-1)^phi * (1-alp) * (kbar(-1)^alp) * (nbar^-alp);      
R =  exp(z) * H(-1)^phi * alp * kbar(-1)^(alp-1) * nbar^(1-alp);  
R = r + delta;

@#for i in 1:NT
(1+tauf)*w_@{i}  = MPL * a_@{i};
@#endfor

Ptot= (
@#for i in 1:NLS
+ P_@{i}
@#endfor
);

Pret= (
@#for i in NT+1:NLS
+ P_@{i}
@#endfor
);


Ptot*beq= (
@#for i in 1:NLS-1
+ P_@{i}(-1)*(1-beta_@{i+1}(-1))*(1+(1-tauk)*r)*ki@{i+1}(-1)
@#endfor
);


penind =  (1/T)*(
@#for i in 1:NT
    + w_@{NT+1-i}(@{-i})
@#endfor
);

penbase =  (1/T)*(
@#for i in 1:NT
    + w_@{i}
@#endfor
);

@#for i in NT+1:NLS
pen_@{i}  = .5 * penbase + .5 * rhop * penind(@{NT-i}) ;       
@#endfor



/*


@#for i in NT+1:NLS
pen_@{i}  =  rhop * (1-tauw) * w_40(@{NT-i}) ;   
@#endfor


// PENSIONS
// Uniform pensions
@#for i in NT+1:NLS
pen_@{i}  =  contr/Pret ;      
@#endfor


// Wage profile based pensions
pen_ind =  (
@#for i in 1:NT
    + w_@{NT+1-i}(@{-i})
@#endfor
);


@#for i in NT+1:NLS
pen_@{i}  =  (rhop/T) * pen(@{NT-i}) ;     
@#endfor
*/


contr = (
@#for i in 1:NT
+ P_@{i} * tauw * w_@{i} * lab_@{i} + P_@{i} * tauf * w_@{i} * lab_@{i}
@#endfor
);

retire = (
@#for i in NT+1:NLS
+  P_@{i} * pen_@{i}
@#endfor
);

//retire=contr;


Tk =  (
@#for i in 1:NLS-1
+ P_@{i}(-1)*tauk*r*ki@{i+1}(-1)           
@#endfor
);

Tc =  (
@#for i in 1:NLS
+ P_@{i}*tauc*c_@{i}                        
@#endfor
);


// Government budget constraint

D=.2*y;
g=.1*y;
Def= g + retire - (contr+Tc+Tk);
D = (1+rd)*D(-1) + Def;
rd = r;

//z = rho * z(-1) + eps;                              

Defratio=Def/y;
Dratio=D/y;

@#for i in 1:NLS
//Beginning loop for consumption

@#if i <= NT

(1+tauc)*c_@{i} = (1 + (1-tauk)*r) * ki@{i}(-1) + (1-tauw)*w_@{i}*lab_@{i} - ki@{i+1} + beq;  

@#else

@#if i != NLS
    //note the -1 rather than +1 WARNING
    (1+tauc)*c_@{i} = (1 + (1-tauk)*r) * ki@{i}(-1) + pen_@{i} - ki@{i+1} + beq;  

@#else

    (1+tauc)*c_@{i} = (1 + (1-tauk)*r) * ki@{i}(-1) + pen_@{i} + beq;                   

@#endif

@#endif

@#endfor  //End loop for consumption




@#for i in 1:NLS

@#if i <= NT

welf_@{i}=log(c_@{i}) + gam1*(1-lab_@{i})^(1-gam)/(1-gam) + beta*beta_@{i+1}*welf_@{i+1}(+1);

@#else

@#if i != NLS

welf_@{i}=log(c_@{i}) + gam1/(1-gam) + beta*beta_@{i+1}*welf_@{i+1}(+1);

@#else

welf_@{NLS}=log(c_@{NLS}) + gam1/(1-gam);

@#endif

@#endif

@#endfor  //End loop for welfare


end;   //----END OF MODEL BLOCK-----




//5) INITIALISATION WITH STEADY STATE VALUES*******************************
initval;

@#for i in 1:NLS
ki@{i} = 5;               
lamda@{i} =1;     
c_@{i}=1;
P_@{i}=.9;
@#endfor

h_1=120;
@#for i in 2:NLS
mig_@{i}=0;
beta_@{i}=.9;
med_@{i}=0;
h_@{i}=120;
@#endfor

@#for i in 1:NT
lab_@{i}=.7;
w_@{i}=.1;
@#endfor

@#for i in NT+1:NLS
pen_@{i}=.2;
@#endfor

ki1 =   0;
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
rhop=.5;
rd=.09;
D=1;
I=1;
Def=1;
penind=1;
penbase=1;
tauc=.1;
tauf=.35;
tauk=.1;
tauw=.3; 
xpop=1;
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
%xpop=1.01;
med_50=50;
%mig_20=10;
%tauk=.2;
%tauw = .15;
%tauf = .55;
%med_20=-50;
end;
steady;
%check;
%resid;

perfect_foresight_setup(periods=200);

perfect_foresight_solver;

run plotsim.m
