@ ------------------------------ RCh8_Gues.G ----------------------

date: May 9, 2008

author: Burkhard Heer


computation of the transition dynamics: guess for time path of
r,w,tau

given initial distribution and aggregate capital stock

computation of the policy functions: value function iteration

prior to the execution of this program you need:

1) run: Rch7_disf.g and store value function, policy function
			and asset grid 

2) adjust the load path in this program so that it is identical to
	the save path of the program Rch7_disf.g

-------------------------------------------------------------------------@


new; 
clear all; 
cls; 
library pgraph; 
graphset; 


h0=hsec;

/* stationary values at the end of transition */


/*
load path="C:\\Dokumente und Einstellungen\\ba6wp1\\Eigene Dateien\\buch\\prog";
save path="C:\\Dokumente und Einstellungen\\ba6wp1\\Eigene Dateien\\buch\\prog";
*/


load gknew=gkden,ag=agden,a=aden;
load v0=vden,aopt0=aoptden;

kknew=sumc(gknew'*ag);

ng=rows(ag); 
na=rows(a); 
amin1=a[1]; 
amax1=a[na];

alpha=0.36; 
beta=0.995; 
delta=0.005; 
sigma=2; 
rep=0.25;

pp=(0.9565~0.0435|0.5~0.5); 
pp1=equivec1(pp); 
nn0=pp1[1];


/* initial distribution and capital stock */
gk=zeros(ng,2);

ng1=sumc(ag.<=300);
/* initial distribution and capital stock */
gkold=gk;
gkold[1:ng1,.]=ones(ng1,2)/(ng1);
gkold=gkold.*pp1';
kkold=sumc(gkold'*ag); /* initial aggregate capital stock */ 

nt=2000;           /* number of transition periods */

neg=-1e10; 
nq=50;

tol=0.001;  /* deviation of k */

tol1=1e-7; /* golden section search accuracy */

crit=1+tol; 

psi0=0.9; /* update of gamma */

eps=0.01;   /* test of a[1]+eps*(a[2]-a[1]) is superior to a[1] */

@ -------------------------------------------

Guess a sequence for r,w,tau,b

------------------------------------------ @


kg=exp(seqa(ln(kkold),(ln(kknew)-ln(kkold))/(nt-1),nt));


@ ---------------------------------------------

Iteration to approximate time path for K_t

--------------------------------------------- @


q=0; 
do until (crit<tol) or (q==nq);
    q=q+1;
    "iteration: " q;
    "time elapsed: " etstr(hsec-h0);
    "percentage error in K: " crit;

	/* time path for r, w, tau, b given sequence kg */
    rg=alpha.*kg^(alpha-1).*nn0^(1-alpha)-delta;
    wg=(1-alpha).*kg^(alpha).*nn0^(-alpha);
    taug=wg.*rep*pp1[2]./(wg.*nn0+rep*pp1[2].*wg+rg.*kg);
    bg=(1-taug).*rep.*wg;


@ ----------------------------------------------------------------

Backwards Iteration of the policy function/value function

----------------------------------------------------------------- @

    ve=zeros(na,nt);
    vu=zeros(na,nt);
    ve[.,nt]=v0[.,1];
    vu[.,nt]=v0[.,2];
    aopte=zeros(na,nt);
    aoptu=zeros(na,nt);
    aopte[.,nt]=aopt0[.,1];
    aoptu[.,nt]=aopt0[.,2];

    i=nt;
    do until i==1;
    i=i-1;  
	"value function iteration: q~i~crit";;
	 q~i~crit';
    e=0;    /* iteration over the employment status */
    do until e==2;  /* e=1 employed, e=2 unemployed */
        e=e+1;
        j=0;        /* iteration over asset grid a in period t */
        l0=0;
        do until j==na;
            j=j+1;
            l=l0;
            vmin=neg;
            ax=a[1]; bx=a[1]; cx=a[na];
            do until l==na; /* iteration over a' in period i */
                l=l+1; 
                if e==1;
                    c=(1+(1-taug[i])*rg[i])*a[j]+(1-taug[i])*wg[i]-a[l];
                 else;
                    c=(1+(1-taug[i])*rg[i])*a[j]+bg[i]-a[l];
                endif;
                if c>0;
                    v1=bellman(a[j],a[l],e);
                    if v1>vmin;
                        if l==1;
                           ax=a[1]; bx=a[1]; cx=a[2];
                        elseif l==na;
                           ax=a[na-1]; bx=a[na]; cx=a[na];
                        else;
                           ax=a[l-1]; bx=a[l]; cx=a[l+1];
                        endif;
                        vmin=v1;
                        l0=l-1;
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
                    if e==1;
                        aopte[j,i]=a[1];
                        ve[j,i]=bellman(a[j],a[1],e);
                    else;
                        aoptu[j,i]=a[1];
                        vu[j,i]=bellman(a[j],a[1],e);
                    endif;
               else;
                    if e==1;
                        aopte[j,i]=golden(&value1,ax,bx,cx,tol1);
                        ve[j,i]=bellman(a[j],aopte[j,i],e);
                    else;
                        aoptu[j,i]=golden(&value1,ax,bx,cx,tol1);
                        vu[j,i]=bellman(a[j],aoptu[j,i],e);
                    endif;
               endif;
            elseif bx==cx;  /* boundary optimum, bx=cx=a[n] */
               bx=cx-eps*(a[na]-a[na-1]);
               if value1(bx)<value1(cx);
                  if e==1;
                    aopte[j,i]=a[na];
                  else;
                    aoptu[j,i]=a[na];
                  endif;
               else;
                  if e==1;
                    aopte[j,i]=golden(&value1,ax,bx,cx,tol1);
                  else;
                    aoptu[j,i]=golden(&value1,ax,bx,cx,tol1);
                  endif;
               endif;
           else;
                if e==1;
                 aopte[j,i]=golden(&value1,ax,bx,cx,tol1);
                else;
                 aoptu[j,i]=golden(&value1,ax,bx,cx,tol1);
                endif;
           endif;
           if e==1;
                ve[j,i]=bellman(a[j],aopte[j,i],e);
           else;
                vu[j,i]=bellman(a[j],aoptu[j,i],e);
           endif;
        endo;   /* i=1,..na */
    endo;   /* e=1,2 */
endo;   /* i=1,..nt */


@ -------------------------------------------------------

Simulation of the distribution function

 ------------------------------------------------------ @

  knew=zeros(nt,1);
  knew[1]=kkold;
  j=1;
  kritg=1;
  do until j==nt;
	"distribution function: q~j~k";;
    q~j~knew[j];
    j=j+1; 
    if j==2;
        gk0=gkold;
    else;
        gk0=gk;
    endif;
    gk=zeros(ng,2);
    e=0;
    do until e==2;
      e=e+1;
      i=0;
      do until i==ng;
         i=i+1;
         k0=ag[i];
       	 if k0<=amin1;
            if e==1;
                k1=aopte[1,j];
            else;
                k1=aoptu[1,j];
            endif;
         elseif k0>=amax1;
            if e==1;
                k1=aopte[na,j];
            else;
                k1=aoptu[na,j];
            endif;
        else;
            if e==1;
                k1=lininter(a,aopte[.,j],k0);
            else;
                k1=lininter(a,aoptu[.,j],k0);
            endif;
        endif;

         if k1<=amin1;
            gk[1,1]=gk[1,1]+gk0[i,e]*pp[e,1];
            gk[1,2]=gk[1,2]+gk0[i,e]*pp[e,2];
         elseif k1>=amax1;
            gk[ng,1]=gk[ng,1]+gk0[i,e]*pp[e,1];
            gk[ng,2]=gk[ng,2]+gk0[i,e]*pp[e,2];
         elseif (k1>amin1) and (k1<amax1);
            q1=sumc(ag.<=k1)+1;
			q1=maxc(q1|2);
			q1=minc(q1|ng);
            n0=(k1-ag[q1-1])/(ag[q1]-ag[q1-1]);
            gk[q1,1]=gk[q1,1]+n0*gk0[i,e]*pp[e,1];
            gk[q1,2]=gk[q1,2]+n0*gk0[i,e]*pp[e,2];
            gk[q1-1,1]=gk[q1-1,1]+(1-n0)*gk0[i,e]*pp[e,1];
            gk[q1-1,2]=gk[q1-1,2]+(1-n0)*gk0[i,e]*pp[e,2];
         endif;
      endo;	/* j = ag[1],ag[2],.. */
    endo;	/* e=1,2 */
    gk=gk/sumc(sumc(gk));
    if j==10;
        gk10=gk;
    endif;
    if j==100;
        gk100=gk;
    endif;
    kk1=gk'*ag;
    if j<nt;
        knew[j]=sumc(kk1);
    else;
        knew[j]=kknew;
    endif;
    if j==nt;
        save gk2000=gk;
    endif;

  endo;   /* q1=1,.., invariant distribution */
  crit=meanc(abs(knew-kg)./kg);
  kg=psi0*kg+(1-psi0)*knew;

    save aggues=ag,gkoldgues=gkold,gk10gues=gk10,gk100gues=gk100,gkgues=gk;
endo; /* q=1,..,nq */

/* figure 8.5 */
/* transformation of discrete probabilites gk into density */
da=a[2]-a[1];
mf=ng/da;   /* number of grid points divided by average size of the interval */

    "time elapsed: " etstr(hsec-h0);
"iteration completed: " q; "deviation of capital stock: " crit; wait;
xlabel("wealth a");
ylabel("density");
title("distribution dynamics");
_plegstr="t=1\000t=10\000t=100\000t=2000";
_plegctl=1;
xy(ag,(gkold[.,1]+gkold[.,2])/mf~(gk10[.,1]+gk10[.,2])/mf~(gk100[.,1]+gk100[.,2])/mf~(gk[.,1]+gk[.,2])/mf);
wait;
title("");
_plegstr="t=2000\000new stationary density";
xy(ag,(gk[.,1]+gk[.,2])/mf~(gknew[.,1]+gknew[.,2])/mf);
wait;
save gknewgues=gknew;

@ -----------------------------------------------------

Procedures:

u(x) -- utility function

value(x,y) -- next-period value

value1(x) -- right-hand side of the bellman equation
                given a[j] and epsilon

bellman(a0,a1,e) -- right-hand side of the bellman equation

lininter -- linear interpolation

golden -- golden section search

------------------------------------------------------  @

proc u(x);
   retp(x^(1-sigma)/(1-sigma));
endp;


proc value1(x);
    retp(bellman(a[j],x,e));
endp;

proc value(x,y);
    if y==1;
        retp(lininter(a,ve[.,i+1],x));
    else;
        retp(lininter(a,vu[.,i+1],x));
    endif;
endp;

proc bellman(a0,a1,y);
   local c;
   if y==1;
      c=(1+(1-taug[i])*rg[i])*a0+(1-taug[i])*wg[i]-a1;
   else;
      c=(1+(1-taug[i])*rg[i])*a0+bg[i]-a1;
   endif;
   if c<0;
      retp(neg);
   endif;
   if a1>a[na];
      retp(a1^2*neg);
   endif;
   retp(u(c)+beta*(pp[y,1]*value(a1,1)+pp[y,2]*value(a1,2)));
endp;

proc lininter(xd,yd,x);
  local j;
  j=sumc(xd.<=x');
  j=minc(j|(rows(xd)-1));
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

