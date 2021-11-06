

periods=100;

@#define NT=6 
@#define NTr=4
@#define NLS=NT+NTr 


%k=[       0    0.06    0.12    0.19    0.26    0.33     0.4    0.47    0.54    0.61    0.68    0.75    0.82    0.89    0.96    1.03     1.1    1.17    1.24    1.31    1.38    1.45    1.52     1.6    1.68    1.76    1.84    1.92       2    2.07
%    2.14    2.21    2.28    2.35    2.42    2.49    2.56    2.63     2.7    2.77    2.83     2.8    2.76    2.71    2.65    2.59    2.52    2.44    2.35    2.24    2.12    1.99    1.84    1.68     1.5     1.3    1.08    0.84    0.58    0.29];

%lamda = [999  434.03  108.51   43.28   23.11   14.35    9.77    7.07    5.36    4.20    3.38    2.78    2.32    1.97    1.70    1.47    1.29    1.14    1.02    0.91    0.82    0.74    0.68    0.61    0.55    0.50    0.46    0.42    0.39    0.36
%    0.34    0.32    0.30    0.28    0.27    0.25    0.24    0.23    0.21    0.20    0.20    0.20    0.21    0.21    0.22    0.23    0.25    0.26    0.28    0.31    0.35    0.39    0.46    0.55    0.69    0.92    1.34    2.21    4.64   18.58];

k = .01*ones(2,30);
lamda = .01*ones(2,30);


//1) VARIABLES************************************************************

var c y r pen w kbar;  

@#for i in 1:NLS
    var ki@{i};         
@#endfor

@#for i in 1:NLS
    var c_@{i};       
@#endfor

@#for i in 1:NLS
    var lamda@{i};   
@#endfor

var check;

//2) EXOGENOUS VARIABLES***************************************************

varexo z;

//predetermined_variables ki; 


//3)PARAMETERS*************************************************************

parameters nbar, lab, alp, beta, delta, rep, tau,  Tr, T, LS;

lab = .1;            //labour supply (constant in this case)
alp=0.3;    // - production function parameter
beta=0.97;  // - discount factor
delta=.02;    // - depreciation rate
rep=0.1;     //- replacement ratio for pension
tau=0.1304; //- rep / 2+rep - income tax
T=6;       // - working period
Tr=4;      // - retirement period
LS=T+Tr;      // - Life span
nbar=T*lab;


//4) MODEL EQUATION*******************************************************
model;


//+++++++++++++++++++++++++Euler equation for worker and retiree+++++++++

@#for i in 1:NLS

    lamda@{i} =  1/c_@{i};

@#endfor

@#for i in 1:NLS-1

    lamda@{i} = beta*lamda@{i+1}(+1) * (1+r(+1));

@#endfor

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ki@{1}=0;
kbar = (
@#for i in 1:NLS
 
    + ki@{i}                             //Aggregate individual capital

@#endfor
);


%y = c + kbar - (1-delta)*kbar(-1);       //Economy wide feasibility constraint

check= y - (c + kbar - (1-delta)*kbar(-1));

y = exp(z) * kbar(-1)^alp * nbar^(1-alp);   //firm production function  

c = (                        
@#for i in 1:NLS

    + c_@{i}                             //Aggregate consumption

@#endfor
);




w =  exp(z)*(1-alp) * (kbar(-1)^alp) * (nbar^-alp);         //wage foc from firm problem
r =  exp(z)*alp * kbar(-1)^(alp-1) * nbar^(1-alp) - delta;  //interest rate foc from firm problem




//pen =  rep * (1-tau) * w(-1) * lab;       //pension equation

T * tau * w * lab = Tr * pen;                            //gov't budget balance


//z = rho * z(-1) + eps;                                    //for stochastic-use this extra expression



@#for i in 1:NLS 
//Beginning loop for consumption

    @#if i <= NT
        
        c_@{i} = (1 + r) * ki@{i}(-1) + (1-tau)*w*lab - ki@{i+1};  //worker consumption

    @#else

        @#if i != NLS
            //note the -1 rather than +1 WARNING 
            c_@{i} = (1 + r) * ki@{i}(-1) + pen - ki@{i+1};        //retiree consumption             

        @#else
            
            c_@{i} = (1 + r) * ki@{i}(-1) + pen;                   //retiree last period consumption                               

        @#endif

    @#endif

@#endfor  //End loop for consumption



end;   //----END OF MODEL BLOCK-----




//5) INITIALISATION WITH STEADY STATE VALUES*******************************
initval; 

//capital and labour supply for each cohorts


ki1    	=	 0;
ki2    	=	 0.20222;
ki3    	=	 0.138872;
ki4    	=	 0.0467899;
ki5    	=	 0.0833141;
ki6    	=	 0.259711;
ki7    	=	 0.483716;
ki8    	=	 0.571518;
ki9    	=	 0.499088;
ki10    	=	 0;
c_1    	=	 0.0151459;
c_2    	=	 0.0253418;
c_3    	=	 0.0424013;
c_4    	=	 0.0709448;
c_5    	=	 0.118703;
c_6    	=	 0.198611;
c_7    	=	 0.332312;
c_8    	=	 0.556017;
c_9    	=	 0.930314;
c_10    	=	 0.930314;

lamda1 	=	 66.0245;
lamda2 	=	 39.4605;
lamda3 	=	 23.5842;
lamda4 	=	 14.0955;
lamda5 	=	 8.42437;
lamda6 	=	 5.03496;
lamda7 	=	 3.00922;
lamda8 	=	 1.79851;
lamda9 	=	 1.07491;
lamda10 	=	 1.07491;
c      	=	 2.28979;
y      	=	 2.31508;
r      	=	 0.742893;
pen    	=	 0.0704403;
w      	=	 0.540186;
kbar   	=	 1.26456;
z = 0;
end; //----END OF INITIALISATION-----

steady;
%check;



endval;
z   =        .01;
end;
steady;
%check;
%resid;


perfect_foresight_setup;

perfect_foresight_solver;

run plotsim.m




