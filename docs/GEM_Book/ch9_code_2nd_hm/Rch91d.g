@ ------------------------------ Rch91d.g ----------------------------

    13.5.2008
    author: Burkhard Heer

    algorithm 9.1

    direct computation of the OLG model in section 9.1

-------------------------------------------------------------------------@


new;
clear all;
cls;
library pgraph;
Macheps=1e-30;
#include ch8_toolbox.src;
GraphSettings;
_plwidth=7;
bsec=hsec;


@ ---------------------------

 parameter

---------------------------- @

beta=0.99;         /* discount factor */
r=0.045;         /* initial value of the interest rate */
s=2;            /* coefficient of relative risk aversion */
alp=0.3;        /* production elasticity of capital */
rep=0.3;        /* replacement ratio */
del=0.1;          /* rate of depreciation */
tr=20;          /* retired */
t=40;           /* working time */
tau=rep/(2+rep);    /* income tax rate */
gam=2;            /* disutility from working */

kmax=10;        /* upper limit of capital grid */
kinit=0;
na=101;          /* number of grid points on assets */
a=seqa(0,kmax/(na-1),na);   /* asset grid */
psi=0.001;      /* parameter of utility function */
phi=0.8;
tol=0.001;       /* percentage deviation of final solution */
tolk=0.001;       /* percentage deviation of final solution for k_1 */
tol1=1e-10;        /* tolerance for golden section search */
neg=-1e10;        /* initial value for value function */
nq1=30;

@ ---------------------------------

initialization

--------------------------------- @

nbar=0.2;
kbar=(alp/(r+del))^(1/(1-alp))*nbar;
kold=100;
nold=2;


/* agents' policy function  */
aopt=zeros(t+tr,1);  /* optimal asset */
copt=zeros(t+tr,1);   /* optimal consumption */
nopt=0.3*ones(t,1);    /* optimal labor supply */


@ ---------------------------------------------

iteration of policy function, wealth distribution,..

---------------------------------------------- @

q=0;
do until q==30 or abs((kbar-kold)/kbar)<tol;
    krit=abs((kbar-kold)/kbar);
    krit0=abs((nbar-nold)/nbar);
    q=q+1;
    w=(1-alp)*kbar^alp*nbar^(-alp);
    r=alp*kbar^(alp-1)*nbar^(1-alp)-del;
    pen=rep*(1-tau)*w*nbar*3/2;
    kold=kbar;
    nold=nbar;
    k60q=zeros(nq1,1);
    k1q=zeros(nq1,1);

    q1=0;
    do until (q1==nq1) or ((q1>5) and (abs(aopt[1])<tolk) );
        "iteration over q1: " q1;
        q1=q1+1;
        if q1==1;
            k60=0.15;
        elseif q1==2;
            k60=0.2;
        else;
            k60=k60q[q1-1]-(k60q[q1-1]-k60q[q1-2])/(k1q[q1-1]-k1q[q1-2])*k1q[q1-1]; 
        endif;

        aopt[t+tr]=k60;
        k60q[q1]=k60;

        /* computation of the decision rules for the retired */
        i=tr;
        do until i==1;      /* all periods t=T+1,T+2,..T+TR */
            i=i-1;
            "old worker: ";
            i;
            x0=aopt[i+t+1];
            {xf,jcode}=FixVMN1(x0,&rfold);
            "jcode ~k'" jcode[1]~xf;
            "q~q1~i~K~N:";
             q~q1~i~kbar~nbar;
            aopt[i+t]=xf;
        endo;


        /* compuation of the decsion rules for the worker */
        i=t+1;
        do until i==1;      /* all periods t=1,2,..T */
            i=i-1;
            "young worker: ";
            i;
            x0=aopt[i+1]|nopt[i];
            {xf,jcode}=FixVMN1(x0,&rfyoung);
            "jcode, k, n " jcode[1]~xf';
            "q~q1~i~K~N:";
             q~q1~i~kbar~nbar;
            aopt[i]=xf[1];
            nopt[i]=xf[2];
        endo;   
    
        k1q[q1]=aopt[1];
        "q~q1: " q~q1;
        "k1: " aopt[1];
        "k1q: ";
        k1q[1:q1]; 
        

    endo;   /* q1 - iteration of k^1,..,k^60 */

    /* computation of the aggregate capital stock and employment nbar */
    knew=meanc(aopt);
    kbar=phi*kold+(1-phi)*knew;
    nnew=meanc(nopt)*2/3;
    nbar=phi*nold+(1-phi)*nnew;
    "nbar~kbar: " nbar~kbar;
endo;   /* q */

"computational time: ";
    ?etstr(hsec-bsec);  wait;
"Fig 9_1: " wait;
xlabel("Generation");
ylabel("Capital Holdings");
xy(seqa(1,1,t+tr),aopt); wait;

"Fig 9_2: " wait;
ylabel("Labor Supply");
xy(seqa(1,1,t),nopt);
wait;

"capital stock: " kbar;
"tau: " tau;
"replacement ratio: " rep;
"pensions: " pen;
"r: " r;
"aggregate employment " nbar;


@ ------------------------------------------------------

Procedures:

u           -- utility function
uc          -- marginal utility of consumption
un          -- marginal utility of leisure
rvalue(k,t) -- value function in time t at point k, linearly
               interpolated
wvalue(k,t) -- of the worker
value1      -- bellman equation for retired
value2      -- bellman equation for worker
golden      -- golden section search, Press et al., 10.1
kw,kr       -- interpolates policy function for next period's
                capital stock, worker and retired

------------------------------------------------------  @

proc u(x,y);
if s==1;
    retp(ln(x+psi)+gam*ln(y));
else;
    retp((((x+psi)*y^gam)^(1-s)-1)/(1-s));
endif;
endp;

proc uc(x,y);
    retp((x+psi)^(-s).*y^(gam*(1-s)));
endp;

proc un(x,y);
    retp(gam*(x+psi)^(1-s).*y^(gam*(1-s)-1));
endp;

proc rfold(x);
    local rf1,k0,k1,k2,c0,c1;
    rf1=0;
    
    k0=x[1];
    k1=aopt[i+t+1];
    if i==tr-1;
        k2=0;
    else;
        k2=aopt[i+t+2];
    endif;
    c0=(1+r)*k0+pen-k1;
    c1=(1+r)*k1+pen-k2;
    rf1=uc(c0,1)/beta-uc(c1,1)*(1+r);
    retp(rf1);
endp;


proc rfyoung(x);
    local rf1,rf2,k0,n0,k1,k2,n1,c0,c1;
    rf1=0;
    rf2=0;
    k0=x[1];
    n0=x[2];
    k1=aopt[i+1];
    k2=aopt[i+2];
    if i==t;
        n1=0;
        c1=(1+r)*k1+pen-k2;
    else;
        n1=nopt[i+1];
        c1=(1+r)*k1+(1-tau)*w*n1-k2;
    endif;
    c0=(1+r)*k0+(1-tau)*w*n0-k1;
    rf1=uc(c0,1-n0)/beta-uc(c1,1-n1)*(1+r);
    rf2=(1-tau)*w*uc(c0,1-n0)-un(c0,1-n0);
    retp(rf1|rf2);
endp;



/* GraphSettings: Changes Gauss' default initialization of graphic routines
**
** Usage: GraphSettings;
*/
  
proc(0)=GraphSettings;

 external matrix  _pltype, _ptitlht, _pnumht,
         _paxht, _pmcolor, _plwidth, _pcsel;
 external proc GraphSet;
 GraphSet;
   _pdate="";        @ do not plot date               @
  _pltype=6;         @ solid lines only               @
 _plwidth=4;         @ line width always 4            @
   _paxht=0.15;      @ size of axis labels            @        
  _pnumht=0.15;      @ size of axis numbering         @
 _ptitlht=0.20;      @ size of titles                 @
 _pmcolor={0,        @ color of axis: black           @
           0,        @ color of axis' numbers         @
           0,        @ color of x-axis label          @
           0,        @ color of y-axis label          @
           0,        @ color of z-axis label          @
           0,        @ color of title                 @
           0,        @ color of boxes                 @
           0,        @ color of date and author       @
           15};      @ background color:        white @
   _pcsel={0,        @ color of first line: black     @
           1,        @ color of second line: dark blue@
           2,        @ color of third line: green     @
           4,        @ color of fourth line: red      @
           3,        @ color of fifth line: light blue@
          10,        @ color of sixth line: light green@
          8};        @ color of seventh line: dark grey@
          
retp;
endp;           

