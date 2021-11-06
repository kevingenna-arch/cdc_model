

@#define NT=40
@#define NTr=20
@#define NLS=NT+NTr 


%k=[       0    0.06    0.12    0.19    0.26    0.33     0.4    0.47    0.54    0.61    0.68    0.75    0.82    0.89    0.96    1.03     1.1    1.17    1.24    1.31    1.38    1.45    1.52     1.6    1.68    1.76    1.84    1.92       2    2.07
%    2.14    2.21    2.28    2.35    2.42    2.49    2.56    2.63     2.7    2.77    2.83     2.8    2.76    2.71    2.65    2.59    2.52    2.44    2.35    2.24    2.12    1.99    1.84    1.68     1.5     1.3    1.08    0.84    0.58    0.29];

%lamda = [999  434.03  108.51   43.28   23.11   14.35    9.77    7.07    5.36    4.20    3.38    2.78    2.32    1.97    1.70    1.47    1.29    1.14    1.02    0.91    0.82    0.74    0.68    0.61    0.55    0.50    0.46    0.42    0.39    0.36
%    0.34    0.32    0.30    0.28    0.27    0.25    0.24    0.23    0.21    0.20    0.20    0.20    0.21    0.21    0.22    0.23    0.25    0.26    0.28    0.31    0.35    0.39    0.46    0.55    0.69    0.92    1.34    2.21    4.64   18.58];

k = .01*ones(2,30);
lamda = .01*ones(2,30);


//1) VARIABLES************************************************************

var  y c r kbar nbar tauw Ptot Pret contr H R MPL beq g Tc Tk retire penind penbase;  

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
a_@{i}=.5+@{i}^2/1000;  
@#endfor



gam= 1.5; // or value of 2
gam1=.07;
alp=0.3;    
beta=0.97;  
delta=0.02;    
deltah=0.02;    
phi=.02;
T=40;    
Tr=20;      
LS=T+Tr;  


//4) MODEL EQUATION*******************************************************
model;

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


ki@{1}=0;
kbar = (
@#for i in 2:NLS

+ P_@{i-1}(-1)*ki@{i}                            

@#endfor
);


%y = c + kbar - (1-delta)*kbar(-1) + g;    
check= y -( c + kbar - (1-delta)*kbar(-1) + g);       

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

retire=contr;


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


g=Tk+Tc;


@#for i in 1:NLS 

@#if i <= NT

(1+tauc)*c_@{i} = (1 + (1-tauk)*r) * ki@{i}(-1) + (1-tauw)*w_@{i}*lab_@{i} - ki@{i+1} + beq;  

@#else

@#if i != NLS

    (1+tauc)*c_@{i} = (1 + (1-tauk)*r) * ki@{i}(-1) + pen_@{i} - ki@{i+1} + beq;                  

@#else

    (1+tauc)*c_@{i} = (1 + (1-tauk)*r) * ki@{i}(-1) + pen_@{i} + beq;                   

@#endif

@#endif

@#endfor  //End loop for consumption



end;   //----END OF MODEL BLOCK-----




//5) INITIALISATION WITH STEADY STATE VALUES*******************************
initval; 


@#for i in 1:NLS
lamda@{i} =1;        
c_@{i}=1;
P_@{i}=.9;
@#endfor

ki1=0;
h_1=800;
@#for i in 2:NLS
ki@{i} = 5;           
mig_@{i}=0;
beta_@{i}=.9;
med_@{i}=0;
h_@{i}=800;
@#endfor

@#for i in 1:NT
lab_@{i}=.7;
w_@{i}=.1;
@#endfor

@#for i in NT+1:NLS
pen_@{i}=.2;
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
rhop=.5;

penind=1;
penbase=1;
tauc=.1;
tauf=.35;
tauk=.1;
tauw=.1;
xpop=1;
z = 0;

end; 
steady;
%check;
%resid;


%shocks;
%var mig_20;
%var med_50;
%periods 1;
%values 100;
%values -100;
%end;


endval;
z   =   .01;
%tauk=.2;
%tauw = .15;
%tauf = .55;
%xpop=2;
%mig_20=100;
%med_50=-100;
end;
steady;
%check;
%resid;

perfect_foresight_setup(periods=100);

perfect_foresight_solver;

run plotsim.m

