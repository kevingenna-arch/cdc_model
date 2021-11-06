@ ------------------------------ Rch83_imo.g ----------------------------

date: March 22, 2007

author: Burkhard Heer

computation of the policy functions and the distribution function

imrohoroglu, 1989, JPE

monotonocity of a' in a

concavity of v

Linear interpolation between grid points

Maximization: golden section search method

input procedures:
golden -- routine for golden section search
lininter -- linear interpolation

the following parameter needs to be specified by the user

case==1: storage economy, amin1=0

case==2: economy with intermediation

-------------------------------------------------------------------------@


new;
clear all;
cls;
library pgraph;
Graphsettings;
_plwidth=7;

/*
load path="c:\\Documents and Settings\\bheer\\My Documents\\buch\\gauss\\RCh7\\";
save path="c:\\Documents and Settings\\bheer\\My Documents\\buch\\gauss\\RCh7\\";
*/

case=2;

h0=hsec;

@ ----------------------------------------------------------------

Step 1: Parameterization and Initialization

----------------------------------------------------------------- @

/* numerical computation */
tolw=1e-4;              /* stopping criterion for value function */
tolg=1e-8; /* distribution function */
neg=-1e10;
tolgs=1e-10;   /* stopping criterium for golden section search */

phi=0.9;
eps=0.01;
beta=0.995;
sigma=1.5;
rl=0;       /* storage rate */
rb=0.01;    /* borrowing rate, a1<0 */
y=1;
th=0.25;

pp=(0.9565~0.0435|0.5~0.5);


pp=(0.9141~0.0234~0.0587~0.0038|
    0.5625~0.3750~0.0269~0.0356|
    0.0608~0.0016~0.8813~0.0563|
    0.0375~0.0250~0.4031~0.5344);


ns=rows(pp);    /* number of states */
pp1=equivec1(pp);

"ergodic distribution: ";
pp1;
wait;


if case==1;
    amin1=0;
elseif case==2;
    amin1=-8;                 /* asset grid */
else;
    "wrong case"; wait;
endif;

amax1=8;
amin1=maxc(-th/rb+0.001|amin1);

na=301;
astep=(amax1-amin1)/(na-1);
a=seqa(amin1,astep,na);

nk=3*na;            /* asset grid for distribution */
agstep=(amax1-amin1)/(nk-1);
ag=seqa(amin1,agstep,nk);

/* initialization of the value function: */
/* s=1: employed, good times
   s=2: unemployed, good times
   s=3: employed, bad times 
   s=4: unemployed, bad times */
   
y=(1|th|1|th); /* income in s=1,..,ns */

v=u((1/beta-1)*a+y[1]);
s=1;
do until s==ns;
    s=s+1;
    v=v~u((1/beta-1)*a+y[s]);
endo;

v=v/(1-beta);                /* agents consume their income */

vold=v;
copt=zeros(na,ns);           /* optimal consumption */
aopt=zeros(na,ns);           /* optimal next-period assets */

ngk=50000;
crit1=1;
kritg=1;
crita=1;
crit=1;
nit=500;

@ ----------------------------------------------------------------

Step 2: Iteration of the value function

----------------------------------------------------------------- @


j=0;
do until (j==nit);
    j=j+1; 
    "j: "; j;
    crita';
    if j/100==round(j/100); xy(a,v[.,1]~vold[.,1]);  endif;
    vold=phi*v+(1-phi)*vold;
  
    aold=aopt;
    s=0;    /* iteration over the state s=e/u, good times/bad times */
    do until s==ns;  
        s=s+1; 
        i=0;        /* iteration over asset grid a in period t */
        l0=0;
        do until i==na;
            i=i+1; 
            l=l0;
            v0=neg;
            ax=a[1]; bx=a[1]; cx=a[na];
            do until l==na; /* iteration over a' in period t*1 */
                l=l+1;
                if a[l]<0;
                    c=a[i]+y[s]-a[l]/(1+rb);
                else;
                    c=a[i]+y[s]-a[l]/(1+rl);
                endif;  
                if c>0;
                    v1=bellman(a[i],a[l]);
                    if v1>v0;
                        if l==1;
                           ax=a[1];  cx=a[2];
                        elseif l==na;
                           ax=a[na-1]; bx=a[na]; cx=a[na];
                        else;
                           ax=a[l-1]; cx=a[l+1];
                        endif;
                        v0=v1;
                        l0=maxc(l-2|0);
                    else;
                        l=na;   /* concavity of value function */
                    endif;
                else;
                    l=na;
                endif;
            endo;   /* l=1,..,na */
            if ax==bx;  /* boundary optimum, ax=bx=a[1]  */
               bx=ax+eps*(a[2]-a[1]);
               if value1(bx)<value1(ax);
                    aopt[i,s]=a[1];
               else;
                    aopt[i,s]=GSS(&value1,ax,cx);
               endif;
            elseif bx==cx;  /* boundary optimum, bx=cx=a[n] */
               bx=cx-eps*(a[na]-a[na-1]);
               if value1(bx)<value1(cx);
                  aopt[i,s]=a[na];
               else;
                  aopt[i,s]=GSS(&value1,ax,cx);
               endif;
           else;
                 aopt[i,s]=GSS(&value1,ax,cx);
           endif;
           v[i,s]=bellman(a[i],aopt[i,s]);
        endo;   /* i=1,..na */ 

    endo;   /* s=1,..ns */
    crit=meanc(abs(vold-v));
    crit=meanc(crit);
    crita=meanc(abs(aopt-aold)); 
    crita=meanc(crita);
    "crit~crita: " crit~crita;

endo;   /* j=1,..nit */

s=0;
do until s==ns;
    s=s+1;
    i=0;
    do until i==na;
        i=i+1;
        if aopt[i,s]<0; 
            r=rb; 
        else; 
            r=rl; 
        endif;
        copt[i,s]=a[i]+y[s]-aopt[i,s]/(1+r);
    endo;
endo;
 

"Net Savings";
"Case: " case; 

GraphSettings;
xtics(amin1,amax1,0.5,1);
_plctrl={0 0 0 20};
_pstype=2;
_pltype=1|3|6|6;
title("");
ylabel("Net savings");
xlabel("Individual wealth a");
xy(a,aopt-a);




if case==1;
    save imoa1=a, imoaopt1=aopt, imocopt1=copt, imov1=v;
elseif case==2;
    save imoa2=a, imoaopt2=aopt, imocopt2=copt, imov2=v;
endif;

GraphSettings;
@
s=0;
do until s==ns;
    s=s+1;
    title("value function new old");
    xy(a,v[.,s]~vold[.,s]);
    title("consumption");
    xy(a,copt[.,s]); wait;
endo;
xy(a,copt);
wait;
@

/* iteration to find invariant distribution */
q1=0;
kritg=1;
/* initialization of the distribution functions */   
gk=ones(nk,ns)/nk;
gk=gk.*pp1';

"computation of invariant distribution of wealth..";
    
do until (q1>ngk);
    q1=q1+1; 
    "q1~kritg: " q1~kritg';
    gk0=gk;
    gk=zeros(nk,ns);
    s=0;
    do until s==ns;
        s=s+1;
        i=0;
        do until i==nk;
            i=i+1;
            k0=ag[i];
            if k0<=amin1; 
                k1=aopt[1,s]; 
            elseif k0>=amax1;
                k1=aopt[na,s];
            else;    
                k1=lininter(a,aopt[.,s],k0);
            endif;
            if k1<=amin1;
                s1=0;
                do until s1==ns;
                    s1=s1+1;
                    gk[1,s1]=gk[1,s1]+gk0[i,s]*pp[s,s1];
                endo;
            elseif k1>=amax1;
                s1=0;
                do until s1==ns;
                    s1=s1+1;
                    gk[nk,s1]=gk[nk,s1]+gk0[i,s]*pp[s,s1];
                endo;
            elseif (k1>amin1) and (k1<amax1);
                j=sumc(ag.<=k1)+1;
                n0=(k1-ag[j-1])/(ag[j]-ag[j-1]);
                s1=0;
                do until s1==ns;
                    s1=s1+1;
                    gk[j,s1]=gk[j,s1]+n0*gk0[i,s]*pp[s,s1];
                    gk[j-1,s1]=gk[j-1,s1]+(1-n0)*gk0[i,s]*pp[s,s1];
                endo;
            endif;
        endo;
    endo;

    gk=gk/sumc(sumc(gk));
    kritg=sumc(abs(gk0-gk));
endo;   /* q1=1,.., invariant distribution */
 
"time elapsed: " etstr(hsec-h0);
"error value function: " crit;
"error distribution: " kritg;
"mean assets: " sumc(gk'*ag);
wait;
title("");
ylabel("Density");
xlabel("Individual wealth a");
xy(ag,gk);
wait;
xy(ag,cumsumc(gk));   
wait;


if case==1;
    save imoag1=ag, imogk1=gk;
elseif case==2;
    save imoag2=ag, imogk2=gk;
endif;

/* plotting the solution */
_plegctl=1;
_plegstr="empl, good times\000unempl, good times\000 empl, bad times\000unempl, bad times";
xtics(amin1,amax1,0.5,1);
_plctrl={0 0 0 20};
_pstype=2;
_pltype=1|3|6|6;
xlabel("Individual wealth a");
ylabel("Value function V(a,s)");
xy(a,v);
wait;

title("");
xlabel("Individual wealth a");
ylabel("Consumption c(a,s)");
xy(a,copt);
save copt;
wait;

title("change in asset level");
xy(a,aopt-a);
save aopt;
wait;

title("invariant capital distribution");
xy(ag,gk/agstep);
save ag,gk,agstep;
wait;

title("next-period assets");
xy(a,aopt);
wait;

/* entries in table 2, A. Imrohoroglu */
/* assets borrowed */
z=(ag.<=0);
nz=sumc(z); /* number of entries in ag with ag<=0 */
borrowa=-gk[1:nz,.]'*ag[1:nz];
borrowa=sumc(borrowa);
savingsa=gk[nz+1:nk,.]'*ag[nz+1:nk];
savingsa=sumc(savingsa);
"borrowing: " borrowa;
"total savings: " savingsa-borrowa;
"asset stored: " savingsa;
"credit market costs: " borrowa*rb;
"income: " pp1'*y;
"consumption: " pp1'*y-borrowa*rb;

@  ----------------------------  procedures -----------


u(x) -- utility function

value(a,e) -- returns the value of the value function for asset
                a and employment status e

value1(x) -- given a=a[i] and epsilon=e, returns the
             value of the bellman equation for a'=x

bellman -- value for the right-hand side of the Bellman equation

------------------------------------------------------- @

proc u(x);
   retp(x^(1-sigma)/(1-sigma));
endp;

proc value(x,y);
    if x>=amax1; retp(vold[na,y]); endif;
    if x==amin1; retp(vold[1,y]); endif;
    if x<amin1; retp(neg*(1+x^2)); endif;
    retp(lininter(a,vold[.,y],x));
endp;

proc value1(x);
    retp(bellman(a[i],x));
endp;


proc bellman(a0,a1);
   local c,s1,r,bell;
   if a1<0; r=rb; else; r=rl; endif;
   c=a0+y[s]-a1/(1+r);
   if c<=0;
      retp(neg);
   endif;
   if a1>a[na];
      retp(a1^2*neg);
   endif;
   bell=u(c);
   s1=0;
   do until s1==ns;
        s1=s1+1;
        bell=bell+beta*pp[s,s1]*value(a1,s1);
   endo;
   retp(bell);
endp;

proc lininter(xd,yd,x);
  local j,nj;
  j=sumc(xd.<=x');
  if j==0; retp(yd[1]); endif;
  if j==rows(xd); nj=rows(xd); retp(yd[nj]); endif;
  retp(yd[j]+(yd[j+1]-yd[j]).*(x-xd[j])./(xd[j+1]-xd[j]));
endp;

proc golden(&f,ay,by,cy,tol);
    local f:proc,x0,x1,x2,x3,xmin,r1,r2,f1,f2;
    r1=0.61803399; r2=1-r1;
    x0=ay;
    x3=cy;
    if abs(cy-by)<=abs(by-ay);
        x1=by; x2=by+r2*(cy-by);
    else;
        x2=by; x1=by-r2*(by-ay);
    endif;
    f1=-f(x1);
    f2=-f(x2);
    do until abs(x3-x0)<=tol*(abs(x1)+abs(x2));
        if f2<f1;
            x0=x1;
            x1=x2;
            x2=r1*x1+r2*x3;
            f1=f2;
            f2=-f(x2);
        else;
            x3=x2;
            x2=x1;
            x1=r1*x2+r2*x0;
            f2=f1;
            f1=-f(x1);
        endif;
    endo;
    if f1<=f2;
        xmin=x1;
        else;
        xmin=x2;
    endif;
    retp(xmin);
endp;

proc equivec1(p);
    local n,x;
    n=rows(p);
    p=diagrv(p,diag(p)-ones(n,1));
    p=p[.,1:n-1]~ones(n,1);
    x=zeros(n-1,1)|1;
    retp((x'*inv(p))');
endp;


/* GSS: Finds the maximizer of a single peaked function of a singel variable
**      in the interval [xl,xu]
**
** Usage: x=GSS(&f,xl,xu);
**
** Input: &f: pointer the the procedure that returns f(x)
**        xl: scalar, the lower boundary of the interval
**        xu: scalar, the upper boundary of the interval
**
** Output: x: scalar, the maximizer of f(x) in [xl,xu]
*/


proc(1)=GSS(&f,xl,xu);

   local tol, p, q, a, b, c, d, fb, fc, f:proc;

  /* Compute the parameters of the problem */
   tol=tolgs;
   p=(sqrt(5)-1)/2;
   q=1-p;

  /* Compute the initial interval [a,d] and the points b and c
  ** that divide it
 */
   a=xl; d=xu; b=p*a+q*d; c=q*a+p*d;

 /* Compute the function value at b and c */
   fb=f(b); fc=f(c);
   
 /* Iterate so that the inverval gets small and smaller */
   do while abs(d-a)>tol*maxc(1|(abs(a)+abs(c)));

      if fb<fc; /* choose [b,d] as next interval */
         a=b;
         b=c;
        fb=fc;
         c=p*b+q*d;
        fc=f(c);
      else; /* choose [a,c] as next interval */        
         d=c;
         c=b;
        fc=fb;
         b=p*c + q*a;         
         fb=f(b);
      endif;
   endo;

   if fb>fc; retp(b); else; retp(c); endif;

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

