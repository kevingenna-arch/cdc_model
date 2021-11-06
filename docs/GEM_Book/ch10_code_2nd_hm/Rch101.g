@ ------------------------------ Rch101.g ----------------------------

    22.10.2008, Burkhard Heer

    finite value function approximation
    
    exogenous labor supply: solving non-

    linear interpolation between grid points

    Golden Section Search

    Adjust load/save PATH in line 30-31!!

-------------------------------------------------------------------------@


new;
clear all;
cls;
library pgraph;
library pgraph,user;
GraphSettings;
_plwidth=7;


/* ADJUST LOAD/SAVE PATH HERE !! */

/*
load path="d:\\buch\\prog\\"; 
save path="d:\\buch\\prog\\";
*/


@ ---------------------------

 parameter

---------------------------- @


/* survival probabilities of the 1-60 year old */
load sp1,ef;
/* sp1=ones(59,1)|0; */

ef=ef/meanc(ef);
mass=ones(60,1);       /* measure of living persons */

/* ef=ones(40,1); */

j=1;
do until j==60;
    j=j+1;
    mass[j]=mass[j-1]*sp1[j-1];
endo;

title("mean effective labor supply");
xlabel("age");
xy(seqa(20,1,rows(ef)),ef);
wait;
title("survival probability");
xy(seqa(20,1,rows(sp1)-1),sp1[1:59]);
wait;

b=1.011;         /* discount factor */
r=0.045;         /* initial value of the interest rate */
s=1.5;            /* coefficient of relative risk aversion */
alp=0.36;        /* production elasticity of capital */
rep=0.3;        /* replacement ratio */
del=0.06;          /* rate of depreciation */
tr=20;          /* retired */
t=40;           /* working time */
hbar=0.3;       /* fixed shift length */

kmin=0;
kmax=40;        /* upper limit of capital grid */
kinit=0;
na=601;          /* number of grid points on assets */
a=seqa(kmin,(kmax-kmin)/(na-1),na);   /* asset grid */
nag=2*na;          /* number of grid points for distribution  */
ag=seqa(kmin,(kmax-kmin)/(nag-1),nag);   /* asset grid for distribution*/

@  -------------------------------------------------------------

productivity process

markov transition matrix py

productivity distribution of 20-year-old: mu[1,.]

------------------------------------------------------     @

lamb0=.96;  /* autoregressive parameter */
sigmay1=0.38;   /* variance for 20-year old, log earnings */
sigmae=0.045;   /* earnings disturbance term variance */

sy=sqrt(sigmae/(1-lamb0^2));    /* earnings variance */
x0=1;   
nsim=10;
ny=9;   /* number of productivities */
m=2;    
ye=seqa(-m*sy,(2*m*sy)/(ny-1),ny);  /* grid over productivities */
{x,py} = armc(nsim,x0,lamb0,sigmae,ny,m);

muy=zeros(t,ny); 
w=ye[2]-ye[1];
muy[1,1]=cdfn((ye[1]+w/2)/sqrt(sigmay1));
muy[1,ny]=1-cdfn((ye[ny]-w/2)/sqrt(sigmay1));
i=1;
do until i==ny-1;
    i=i+1;
    muy[1,i]=cdfn((ye[i]+w/2)/sqrt(sigmay1))-cdfn((ye[i]-w/2)/sqrt(sigmay1));
endo;



@ ------------------------------------------

 computation of the labor supply 

 and the Gini coefficient of earnings (workers)
  
------------------------------------------  @

i=1;
earn=muy[1,.]'~exp(ye)*ef[1]*hbar;

do until i==t;
    i=i+1;
    muy[i,.]=muy[i-1,.]*py*sp1[i];
    earn=earn|muy[i,.]'~exp(ye)*ef[i]*hbar;
endo;
nbar=sumc(sumc(muy.*exp(ye').*ef))*hbar/sumc(mass[1:t]);
"nbar: " nbar;
earn=sortc(earn,2);
earn[.,1]=earn[.,1]/sumc(earn[.,1]);
"sumc(earn[.,1])=1 " sumc(earn[.,1]);
earnm=earn[.,1]'*earn[.,2];

ea=(0~0|
    0.05~0|
    0.1~0|
    0.2~0|
    0.4~0.032|
    0.6~0.125|
    0.8~0.233|
    0.9~0.6139-0.1238-0.1637-0.1476|
    0.95~0.1238|
    0.99~0.1637|
    1~1);
    
    
ea[2:rows(ea)-1,2]=cumsumc(ea[2:rows(ea)-1,2]);
ea=ea|ones(rows(earn)-rows(ea),2);
    

xlabel("proportion of workers");
ylabel("proportion of earnings");
title("Lorenz curve for earnings");
plegctl=1;
_plegctl={ 40 0.02 0.1 0.7};
_pltype=6|1|3;
_pgrid={1,1};
_plctrl={0 0 0 1};
_pstype=1;
xtics(0,1,0.1,0.05);
ytics(0,1,0.1,0.05);
 _plegstr="model\000US Earnings\000equal distribution";
xy(cumsumc(earn[.,1])~ea[.,1]~ea[.,1],cumsumc(earn[.,2])/sumc(earn[.,2])~ea[.,2]~ea[.,1]);
wait;
graphsettings;
save earn,ea;

j=1;
fj=earn[.,1];
hj=zeros(rows(earn),1);
hj[1]=fj[1]*earn[1,2]/earnm;
gini=1-hj[1]*earn[1,1];
do until j==rows(earn);
  j=j+1;
  hj[j]=hj[j-1]+fj[j]*earn[j,2]/earnm;
  gini=gini-(hj[j]+hj[j-1])*fj[j];
endo;
"gini earnings: " gini;


wait;


@  ----------------------------

computational parameters

--------------------------------- @

psi=0.001;      /* parameter of utility function */
eps=0.05;
phi=0.8;
tol=0.001;       /* percentage deviation of final solution */
tol1=1e-10;        /* tolerance for golden section search */
neg=-1e10;        /* initial value for value function */

@ ---------------------------------

initialization

--------------------------------- @


kbar=(alp/(r+del))^(1/(1-alp))*nbar;
kold=100;
tau=0.2;    /* income tax rate, initial guess */

@ ---------------------------------------------

iteration of policy function, wealth distribution,..

---------------------------------------------- @

q=0;
do until q==30 or abs((kbar-kold)/kbar)<tol;
krit=abs((kbar-kold)/kbar);

q=q+1;
w=(1-alp)*kbar^alp*nbar^(-alp);
r=alp*kbar^(alp-1)*nbar^(1-alp)-del;
pen=rep*(1-tau)*w*nbar*sumc(mass)/sumc(mass[1:t]);
tau=pen*sumc(mass[t+1:t+tr])/sumc(mass)/(w*nbar);
kold=kbar;

/* retired agents' value function  */
vr=zeros(na,tr); /* value function */
aropt=ones(na,tr);  /* optimal asset */
cropt=zeros(na,tr); /* optimal consumption */
l=0;
do until l==na;
l=l+1;
vr[l,tr]=u(a[l]*(1+r)+pen);
cropt[l,tr]=a[l]*(1+r)+pen;
endo;

/* workers' value function */
vw=zeros(na*ny,t);
awopt=ones(na*ny,t);
cwopt=zeros(na*ny,t);

/* computation of the decision rules for the retired */
i=tr;
do until i==1;      /* all periods t=T+1,T+2,..T+TR */
i=i-1;
l=0;
do until l==na;     /* asset holding in period t */
l=l+1;
y=varput(a[l],"k0");

/* triple ax,bx,cx */
ax=0; bx=-1; cx=-2;
v0=neg;
m=0;
do until ax<=bx and bx<=cx;
m=m+1;
v1=value1(a[m]);
if v1>v0;
    if m==1; ax=a[m]; bx=a[m];
        else;
        bx=a[m]; ax=a[m-1];
     endif;
    v0=v1;
else;
    cx=a[m];
endif;
if m==na; ax=a[m-1]; bx=a[m]; cx=a[m]; endif;
endo;

if ax==bx;  /* randminimum */
        aropt[l,i]=kmin;
    elseif bx==cx;
        aropt[l,i]=a[na];
    else;
        aropt[l,i]=golden(&value1,ax,bx,cx,tol1);
endif;

vr[l,i]=value1(aropt[l,i]);
cropt[l,i]=(1+r)*a[l]+pen-aropt[l,i];

"q~i~l~n:" q~i~l~nbar;
endo;


if q==1 and i==10;
xlabel("capital stock");
title("next period's capital stock");
xy(a,aropt[.,i]);
title("consumption");
xy(a,cropt[.,i]);
title("value function");
xy(a,vr[.,i]);
endif;
endo;

@ -----------------------------------------------------------------------

 compuation of the decsion rules for the worker 

-------------------------------------------------------------------------  @

i=t+1;
do until i==1;      /* all periods t=1,2,..T */
i=i-1;
i;
j=0;
do until j==ny; /* productivity at age t */
    j=j+1;
    l=0;
        m0=0;
    do until l==na;     /* asset holding at age t */
        l=l+1;
        /* triple ax,bx,cx */
        ax=0; bx=-1; cx=-2;
        v0=neg;
        m=m0;
        do until ax<=bx and bx<=cx;
            m=m+1;
            v1=value2(a[m]);
            if v1>v0;
                m0=maxc(0|m-2);
                if m==1; 
                   ax=a[m]; bx=a[m];
                       else;
                        bx=a[m]; ax=a[m-1];
                   endif;
                   v0=v1;
            else;
                   cx=a[m];
            endif;
            if m==na; 
                ax=a[m-1]; bx=a[m]; cx=a[m]; 
            endif;
        endo;

        if ax==bx;  /* boundary optimum, ax=bx=a[1]  */
               bx=ax+eps*(a[2]-a[1]);
               if value2(bx)<value2(ax);
                    awopt[(l-1)*ny+j,i]=a[1];
               else;
                    awopt[(l-1)*ny+j,i]=golden(&value2,ax,bx,cx,tol1);
               endif;
        elseif bx==cx;
               bx=cx-eps*(a[na]-a[na-1]);
               if value2(bx)<value2(cx);
                  awopt[(l-1)*ny+j,i]=a[na];
               else;
                  awopt[(l-1)*ny+j,i]=golden(&value2,ax,bx,cx,tol1);
               endif;
        else;
                awopt[(l-1)*ny+j,i]=golden(&value2,ax,bx,cx,tol1);
        endif;

        vw[(l-1)*ny+j,i]=value2(awopt[(l-1)*ny+j,i]);
                y=0;
        cwopt[(l-1)*ny+j,i]=(1+r)*a[l]+w*exp(ye[j])*hbar*(1-tau)*ef[i]-awopt[(l-1)*ny+j,i];
        "q~i~l~j~k:" q~i~l~j~kbar;
    endo;   /* l */
endo;   /* j */


if (q==1 or q==5) and (i==20);
j=0;
do until j==ny;
j=j+1;
a5=awopt[j,i];
c5=cwopt[j,i];
v5=vw[j,i];

l=1;
do until l==na;
    l=l+1;
    a5=a5|awopt[(l-1)*ny+j,i];
    c5=c5|cwopt[(l-1)*ny+j,i];
    v5=v5|vw[(l-1)*ny+j,i];
endo;
xlabel("capital stock");
title("savings");
xy(a,a5-a);
title("consumption");
xy(a,c5);
title("value function");
xy(a,v5);
endo;
endif;

endo;   /* i */

cls;

/* computation of the aggregate capital stock and employment nbar */
gkw=zeros(nag*ny,t);
gkr=zeros(nag,tr);
kgen=zeros(t+tr,1);
gkw[1:ny,1]=muy[1,1:ny]';


gk=zeros(nag,t+tr);
gk[1,1]=mass[1];
i=0;
do until i==t-1;
i=i+1; i~kgen[i];
 j=0;
    do until j==ny;
      j=j+1;
      l=0;
      do until l==nag;
        l=l+1;
        if ag[l]<=kmin; 
            k1=awopt[j,i];
        elseif ag[l]>=kmax;
            k1=awopt[(na-1)*ny+j,i];
        else;
            j0=sumc(a.<ag[l])+1;
            j0=minc(j0|nag);
            n0=(ag[l]-a[j0-1])/(a[j0]-a[j0-1]);
            k1=(1-n0)*awopt[(j0-2)*ny+j,i]+n0*awopt[(j0-1)*ny+j,i];
        endif;
        
        if k1<=kmin;
            m=0;
            do until m==ny;
                m=m+1;
                gkw[m,i+1]=gkw[m,i+1]+py[j,m]*sp1[i]*gkw[(l-1)*ny+j,i];
            endo;
        elseif k1>=kmax;
            m=0;
            do until m==ny;
                m=m+1;
                gkw[(nag-1)*ny+m,i+1]=gkw[(nag-1)*ny+m,i+1]+py[j,m]*sp1[i]*gkw[(l-1)*ny+j,i];
            endo;
        elseif (k1>kmin) and (k1<kmax);
            j0=sumc(ag.<k1)+1;
            j0=minc(j0|nag);
            n0=(k1-ag[j0-1])/(ag[j0]-ag[j0-1]);
            m=0;
            do until m==ny;
                m=m+1;
                gkw[(j0-2)*ny+m,i+1]=gkw[(j0-2)*ny+m,i+1]+(1-n0)*py[j,m]*sp1[i]*gkw[(l-1)*ny+j,i];
                gkw[(j0-1)*ny+m,i+1]=gkw[(j0-1)*ny+m,i+1]+n0*py[j,m]*sp1[i]*gkw[(l-1)*ny+j,i];
            endo;            
        endif;
      endo; /* l */
    endo; /* j */
    kg0=reshape(gkw[.,i+1],nag,ny);
    kg0=sumc(kg0');
    gk[.,i+1]=kg0;
    kgen[i+1]=kg0'*ag;
endo;

i=t;
j=0;
    do until j==ny;
      j=j+1;
      l=0;
      do until l==nag;
        l=l+1;
        if ag[l]<=kmin; 
            k1=awopt[j,i];
        elseif ag[l]>=kmax;
            k1=awopt[(na-1)*ny+j,i];
        else;
            j0=sumc(a.<ag[l])+1;
            n0=(ag[l]-a[j0-1])/(a[j0]-a[j0-1]);
            k1=(1-n0)*awopt[(j0-2)*ny+j,i]+n0*awopt[(j0-1)*ny+j,i];
        endif;
       
        if k1<=kmin;
                gkr[1,1]=gkr[1,1]+sp1[i]*gkw[(l-1)*ny+j,i];
        elseif k1>=kmax;
                gkr[nag,1]=gkr[nag,1]+sp1[i]*gkw[(l-1)*ny+j,i];
        elseif (k1>kmin) and (k1<kmax);
            j0=sumc(ag.<k1)+1;
            n0=(k1-ag[j0-1])/(ag[j0]-ag[j0-1]);
                gkr[j0-1,1]=gkr[j0-1,1]+(1-n0)*sp1[i]*gkw[(l-1)*ny+j,i];
                gkr[j0,1]=gkr[j0,1]+n0*sp1[i]*gkw[(l-1)*ny+j,i];     
        endif;
      endo; /* l */
    endo; /* j */    
    kgen[t+1]=gkr[.,1]'*ag; 
    t+1~kgen[t+1];


i=0;
do until i==tr-1;
    i=i+1; i+t;
      l=0;
      do until l==nag;
        l=l+1;
        if ag[l]<=kmin; 
            k1=aropt[1,i];
        elseif ag[l]>=kmax;
            k1=awopt[na,i];
        else;
            j0=sumc(a.<ag[l])+1;
            n0=(ag[l]-a[j0-1])/(a[j0]-a[j0-1]);
            k1=(1-n0)*aropt[j0-1,i]+n0*aropt[j0,i];
        endif;
        if k1<=kmin;
            gkr[1,i+1]=gkr[1,i+1]+sp1[i+t]*gkr[l,i];
        elseif k1>=kmax;
            gkr[nag,i+1]=gkr[nag,i+1]+sp1[i+t]*gkr[l,i];
        elseif (k1>kmin) and (k1<kmax);
            j0=sumc(ag.<k1)+1;
            n0=(k1-ag[j0-1])/(ag[j0]-ag[j0-1]);
            gkr[j0-1,i+1]=gkr[j0-1,i+1]+(1-n0)*sp1[i+t]*gkr[l,i];
            gkr[j0,i+1]=gkr[j0,i+1]+n0*sp1[i+t]*gkr[l,i];
        endif;
      endo; /* l */
      kgen[i+t+1]=gkr[.,i+1]'*ag;
    i+t+1~kgen[i+t+1];
      save kgen;
endo; /* i */

gk=gk~gkr;
gk=gk/sumc(sumc(gk));
/* computation of the gini coefficient of capital distribution */
gk1=sumc(gk');
gkbar=gk1'*ag;
j=1;
fj=zeros(nag,1);
fj[1]=gk1[1];
hj=zeros(nag,1);
hj[1]=fj[1]*ag[1]/gkbar;
gini=1-hj[1]*gk1[1];
do until j==nag;
  j=j+1;
  fj[j]=gk1[j];
  hj[j]=hj[j-1]+fj[j]*ag[j]/gkbar;
  gini=gini-(hj[j]+hj[j-1])*fj[j];
endo;
"gini: " gini;
save gini;
    
if q==1 or q==5 or q>=10;
xlabel("a");
title("distribution of capital stock");
xy(ag,gk1);

wm=gk1.*ag;
wm=gk1~wm;

we=(0~0|
    0.05~0|
    0.1~0|
    0.2~0|
    0.4~0.0174|
    0.6~0.0572|
    0.8~0.1343|
    0.9~0.7949-0.1262-0.2395-0.2955|
    0.95~0.1262|
    0.99~0.2395|
    1~1);

we[2:rows(we)-1,2]=cumsumc(we[2:rows(we)-1,2]);
  we=we|ones(rows(wm)-rows(we),2);
          
title("Lorenz curve of wealth");
xlabel("proportion of households");
ylabel("proportion of wealth");
_plegctl=1;
_plegctl={ 40 0.02 0.1 0.7};
_pltype=6|1|3;
_pgrid={1,1};
_plctrl={0 0 0 1};
_pstype=1;
xtics(0,1,0.1,0.05);
ytics(0,1,0.1,0.05);
 _plegstr="model\000US wealth\000equal wealth";
xy(cumsumc(wm[.,1])~we[.,1]~we[.,1],cumsumc(wm[.,2])/sumc(wm[.,2])~we[.,2]~we[.,1]);
graphsettings;
save wm,we;


title("");
xlabel("generation");
ylabel("average capital holdings");
xy(seqa(20,1,t+tr+1),((kgen./mass)|0));
endif;


"sumc(mass[1:t])=sumc(sumc(gkw))?";
sumc(mass[1:t]); sumc(sumc(gkw));
"sumc(mass[t+1:tr])=sumc(sumc(gkr))?";
sumc(mass[t+1:tr]); sumc(sumc(gkr));

knew=sumc(kgen)/sumc(mass);
kbar=phi*kold+(1-phi)*knew;
endo;   /* q */

wm=gk1.*ag;
wm=gk1~wm;

we=(0~0|
    0.05~0|
    0.1~0|
    0.2~0|
    0.4~0.0174|
    0.6~0.0572|
    0.8~0.1343|
    0.9~0.7949-0.1262-0.2395-0.2955|
    0.95~0.1262|
    0.99~0.2395|
    1~1);

we[2:rows(we)-1,2]=cumsumc(we[2:rows(we)-1,2]);
  we=we|ones(rows(wm)-rows(we),2);
          
title("Lorenz curve of wealth");
xlabel("proportion of households");
ylabel("proportion of wealth");
_plegctl=1;
_plegctl={ 40 0.02 0.1 0.7};
_pltype=6|1|3;
_pgrid={1,1};
_plctrl={0 0 0 1};
_pstype=1;
xtics(0,1,0.1,0.05);
ytics(0,1,0.1,0.05);
 _plegstr="model\000US wealth\000equal wealth";
xy(cumsumc(wm[.,1])~we[.,1]~we[.,1],cumsumc(wm[.,2])/sumc(wm[.,2])~we[.,2]~we[.,1]);
graphsettings;
wait;

title("");
xlabel("generation");
title("average capital holdings");
ylabel("");
xy(seqa(20,1,t+tr+1),((kgen./mass)|0));
wait;

"capital stock: " kbar;
"tau: " tau;
"replacement ratio: " rep;
"pensions: " pen;
"r: " r;



@ ------------------------------------------------------

Procedures:

u           -- utility function
uc          -- marginal utility of consumption

rvalue(k,t) -- value function in time t at point k, linearly
               interpolated
wvalue(k,t) -- of the worker
value1      -- bellman equation for retired
value2      -- bellman equation for worker
golden      -- golden section search, Press et al., 10.1


------------------------------------------------------  @

proc u(x);
if s==1;
    retp(ln(x+psi));
else;
    retp(((x+psi)^(1-s)-1)/(1-s));
endif;
endp;

proc uc(x);
    retp((x+psi)^(-s));
endp;

proc rvalue(k,t);
    local k0,n1,n2;
    k0=k/kmax*(na-1)+1;
    n2=floor(k0);
    n1=k0-n2;
    if k<=0; retp(vr[1,t]-(1-k0)*(vr[2,t]-vr[1,t])); endif;
    if k>=kmax; retp(vr[na,t]); endif;
    retp((1-n1)*vr[n2,t]+n1*vr[n2+1,t]);
endp;

proc wvalue(k,t);
    local k0,n1,n2,r,ww;
    k0=k/kmax*(na-1)+1;
    n2=floor(k0);
    n1=k0-n2;
    ww=0;
    if k<=0; 
        r=0;
        do until r==ny;
            r=r+1;
            ww=ww+py[j,r]*(vw[r,t]-(1-k0)*(vw[na+r,t]-vw[r,t])); 
        endo;
    elseif k>=kmax; 
        r=0;
        do until r==ny;
            r=r+1;
            ww=ww+py[j,r]*vw[(na-1)*ny+r,t]; 
        endo;
    else;
        r=0;
        do until r==ny;
            r=r+1;
            ww=ww+py[j,r]*((1-n1)*vw[(n2-1)*ny+r,t]+n1*vw[n2*ny+r,t]); 
        endo;
    endif;
    retp(ww);
endp;


proc value1(x);
    local c;
    c=(1+r)*a[l]+pen-x;
    if c<=0; retp(-1e10); endif;
    retp(u(c)+sp1[t+i]*b*rvalue(x,i+1));
endp;

proc value2(x);
    local y,xf,fvp,jk,jcode,c,k0;
    k0=a[l];
    c=(1+r)*k0+(1-tau)*w*ef[i]*exp(ye[j])*hbar-x;
    if c<=0; retp(-1e10); endif;
    if i==t;
        retp(u(c)+sp1[i]*b*rvalue(x,1));
    else;
        retp(u(c)+sp1[i]*b*wvalue(x,i+1));
    endif;
endp;


/* Procedure zur Berechnung vom Minimum - Golden Section Search
   Inputs:  &f  Funktion, fï¿½r die das Minimum bestimmt werden soll
            0,1,2 stuetzstellen mit x0<x1<x2
   Output:      minimum
   Remark:      Siehe auch press et al., chapter 10.1 */

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


/* procedue for the discrete approximation of an AR(1) process */
/* input: n - scalar, number of values for the simulaton
          x0 - scalar, initial value, 1<x0<NN
          a - autoregression coefficient
          s2 - scalar, variance of the disturbance epsilon
          NN - scalar, number of grid points
          m - scalar, spread parameter, e.g. m=3
   output: (n,1) vector
*/

proc (2) = armc(n,x0,a,s2,NN,m);
    local f1,f2,i,ix,p,s,sy,w,y;
    sy=sqrt(s2/(1-a^2));
    y=seqa(-m*sy,(2*m*sy)/(NN-1),NN);
    s=sqrt(s2);
    w=y[2]-y[1];
    f1=zeros(NN,1)~cdfn((y[2:NN]'-a*y-w/2)/s);
    f2=cdfn((y[1:NN-1]'-a*y+w/2)/s)~ones(NN,1);
    p=f2-f1;
    ix=zeros(n,1);
    ix[1]=x0;
    i=1;
    do until i==n;
        i=i+1;
        ix[i]=rndintw(1,p[ix[i-1],.]');
    endo;
    retp(y[ix],p);
endp;

proc rndintw(n,p);
    local k,u;
    k=rows(p);
    u=rndu(n,1);
    retp(subscat(rndu(n,1),cumsumc(p),seqa(1,1,k)));
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

