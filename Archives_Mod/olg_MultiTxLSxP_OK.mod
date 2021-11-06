

@#define NT=40
@#define NTr=20
@#define NLS=NT+NTr 


k = .01*ones(2,30);
lamda = .01*ones(2,30);


//1) VARIABLES************************************************************

var c y r w kbar nbar contr R beq Ptot Pret;  

@#for i in 1:NLS
    var ki@{i};        
@#endfor

@#for i in 1:NT
    var lab_@{i};  
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
    var P_@{i};   
@#endfor

var check;


//2) EXOGENOUS VARIABLES***************************************************

varexo z tau xpop;

@#for i in 2:NLS
   varexo mig_@{i} beta_@{i}; 
@#endfor

//predetermined_variables ki; 


//3)PARAMETERS*************************************************************

parameters alp, beta, delta, Tr, T, LS, gam, gam1;

gam= 2; 
gam1=.07;
alp=0.3;    
beta=0.97; 
delta=.02;  
%rep=0.1;    
T=40;    
Tr=20;    
LS=T+Tr; 


//4) MODEL EQUATION*******************************************************
model;

//+++++++++++++++++++++++++Euler equation for worker and retiree+++++++++

P_1=xpop;

@#for i in 2:NLS
    P_@{i} =  beta_@{i}(-1)*P_@{i-1}(-1) + mig_@{i};
@#endfor


@#for i in 1:NLS
    lamda@{i} =  1/c_@{i};
@#endfor

@#for i in 1:NLS-1
    lamda@{i} = beta*beta_@{i+1}*lamda@{i+1}(+1) * (1+r(+1));
@#endfor

@#for i in 1:NT
gam1*(1-lab_@{i})^(-gam) = lamda@{i}*(1-tau)*w ; 
@#endfor

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ki@{1}=0;
kbar = (
@#for i in 2:NLS
 
    + P_@{i-1}(-1)*ki@{i}                            
@#endfor
);


c = (                        
@#for i in 1:NLS
    + P_@{i}*c_@{i}                     
@#endfor
);

nbar = (                        
@#for i in 1:NT
    + P_@{i}*lab_@{i}                     
@#endfor
);


/*
@#for i in NT+1:NLS
    pen_@{i}  =  rep * (1-tau) * w(@{NT-i}) * lab_@{NT}(@{NT-i}) ;   
@#endfor
*/


contr = (                        
@#for i in 1:NT
    + P_@{i} * tau * w * lab_@{i}
@#endfor
);


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

@#for i in NT+1:NLS
    pen_@{i}  =  contr/Pret ;     
@#endfor


Ptot*beq= (                        
@#for i in 1:NLS-1
    + P_@{i}(-1)*(1-beta_@{i+1}(-1))*(1+r)*ki@{i+1}(-1)
@#endfor
);


/*
retire = (                        
@#for i in NT+1:NLS
    + pen_@{i}
@#endfor
);

retire=contr;
*/

//z = rho * z(-1) + eps;                                 


@#for i in 1:NLS 

    @#if i <= NT
        
        c_@{i} = (1 + r) * ki@{i}(-1) + (1-tau)*w*lab_@{i} - ki@{i+1} + beq;  

    @#else

        @#if i != NLS

            c_@{i} = (1 + r) * ki@{i}(-1) + pen_@{i} - ki@{i+1} + beq;        
        @#else
            
            c_@{i} = (1 + r) * ki@{i}(-1) + pen_@{i} + beq;                           

        @#endif

    @#endif

@#endfor 


%y = c + kbar - (1-delta)*kbar(-1);       
check=y - (c + kbar - (1-delta)*kbar(-1));
y = exp(z) * kbar(-1)^alp * nbar^(1-alp);  

w =  exp(z)*(1-alp) * (kbar(-1)^alp) * (nbar^-alp);   
R =  exp(z)*alp * kbar(-1)^(alp-1) * nbar^(1-alp) ;  

R=r+delta;

end;




//5) INITIALISATION WITH STEADY STATE VALUES*******************************
initval; 

@#for i in 1:NLS
    lamda@{i} =1;    
    c_@{i}=1;
    P_@{i}=1;
@#endfor

ki1 = 0; 
@#for i in 2:NLS
    ki@{i} = 5;           
@#endfor

@#for i in 2:NLS
    mig_@{i}=0;
    beta_@{i}=.999;
@#endfor

@#for i in 1:NT
    lab_@{i}=.7;
@#endfor

@#for i in NT+1:NLS
    pen_@{i}=.2;
@#endfor

R=0.09;
c      	=	 60;
y      	=	 70;
r      	=	 0.02;
w      	=	 1.6;
kbar   	=	 500;
nbar    =   1.75;
contr=4;
beq=.1;

Pret=5;
Ptot=15;

xpop=1;
z = 0;
tau=0.1; 
end;
steady(solve_algo = 0);
%check;
%resid;


endval;
z   =   .02;
%tau = .2;
end;
steady(solve_algo = 0);
%check;
%resid;


perfect_foresight_setup(periods=100);

perfect_foresight_solver;

run plotsim.m

