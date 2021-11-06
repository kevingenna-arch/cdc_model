@ ------------------------------ Rch91v.g ----------------------------

    April 16, 2008, Burkhard Heer

    finite value function approximation

    LINEAR or CUBIC SPLINE interpolation between grid points:
    please choose _IV_method=1 or 2 in line 25!!

    Golden Section Search

-------------------------------------------------------------------------@


new;
clear all;
cls;
library pgraph;
Macheps=1e-30;
#include ch8_toolbox.src;
#include function.src;
graphset;


_VI_method=1;       /* 1 --- linear interpolation, 2 --- cubic spline */

bsec=hsec;

@ ---------------------------

 parameter

---------------------------- @

b=0.99;         /* discount factor */
r=0.045;         /* initial value of the interest rate */
s=2;            /* coefficient of relative risk aversion */
alp=0.3;        /* production elasticity of capital */
rep=0.3;        /* replacement ratio */
del=0.1;          /* rate of depreciation */
tr=20;          /* retired */
t=40;           /* working time */
tau=rep/(2+rep);    /* income tax rate */
gam=2;            /* disutility from working */

kmax=5;        /* upper limit of capital grid */
kinit=0;
na=51;          /* number of grid points on assets */
a=seqa(0,kmax/(na-1),na);   /* asset grid */
psi=0.001;      /* parameter of utility function */
phi=0.8;
tol=0.001;       /* percentage deviation of final solution */
tol1=1e-10;        /* tolerance for golden section search */
neg=-1e10;        /* initial value for value function */

@ ---------------------------------

initialization

--------------------------------- @

nbar=0.2;
kbar=(alp/(r+del))^(1/(1-alp))*nbar;
kold=100;
nold=2;

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

/* retired agents' value function  */
vr=zeros(na,tr); /* value function */
aropt=ones(na,tr);  /* optimal asset */
cropt=zeros(na,tr); /* optimal consumption */
l=0;
do until l==na;
l=l+1;
vr[l,tr]=u(a[l]*(1+r)+pen,1);
cropt[l,tr]=a[l]*(1+r)+pen;
endo;

/* workers' value function */
vw=zeros(na,t);
awopt=ones(na,t);
cwopt=zeros(na,t);
nwopt=zeros(na,t);

/* computation of the decision rules for the retired */
i=tr;
do until i==1;      /* all periods t=T+1,T+2,..T+TR */

    m0=0;

    if _VI_method==2;   /* cubic spline interpolation */
        vr0=vr[.,i];
        y2=CSpline(a,vr0,2,0|0);
    endif;


    i=i-1;
    l=0;
    do until l==na;     /* asset holding in period t */
        l=l+1;
        y=varput(a[l],"k0");


        /* triple ax,bx,cx */
        ax=0; bx=-1; cx=-2;
        v0=neg;
        m=maxc(0|m0-2);
        do until ax<=bx and bx<=cx;
            m=m+1;
            v1=value1(a[m]);
            if v1>v0;
                if m==1; 
                    ax=a[m]; bx=a[m];
                else;
                    bx=a[m]; ax=a[m-1];
                endif;
                v0=v1;
                m0=m;   /* monotonocity of the value function */
            else;
                cx=a[m];
            endif;
            if m==na; 
                ax=a[m-1]; bx=a[m]; cx=a[m]; 
            endif;
    endo;

    if ax==bx;  /* randminimum */
        aropt[l,i]=0;
    elseif bx==cx;
        aropt[l,i]=a[na];
    else;
        aropt[l,i]=golden(&value1,ax,bx,cx,tol1);
    endif;

    vr[l,i]=value1(aropt[l,i]);
    cropt[l,i]=(1+r)*a[l]+pen-aropt[l,i];

    "q~i~l~n:" q~i~l~nbar;
endo;
endo;

/* compuation of the decsion rules for the worker */
i=t+1;
do until i==1;      /* all periods t=1,2,..T */

    m0=0;

    if _VI_method==2;   /* cubic spline interpolation */
        if i==t+1;
            vw0=vr[.,1];
        else;
            vw0=vw[.,i];
        endif;
        y2=CSpline(a,vw0,2,0|0);
    endif;

    i=i-1;
    i;
    l=0;
    do until l==na;     /* asset holding in period t */
        l=l+1;
        y=varput(a[l],"k0");

        /* triple ax,bx,cx */
        ax=0; bx=-1; cx=-2;
        v0=neg;
        m=maxc(0|m0-2);
        do until ax<=bx and bx<=cx;
            m=m+1;
            v1=value2(a[m]);
            if v1>v0;
                if m==1; 
                    ax=a[m]; bx=a[m];
                else;
                    bx=a[m]; ax=a[m-1];
                endif;
                v0=v1;
                m0=m;
            else;
                cx=a[m];
            endif;
            if m==na; 
                ax=a[m-1]; bx=a[m]; cx=a[m]; 
            endif;
        endo;

        if ax==bx;  /* randminimum */
            awopt[l,i]=0;
        elseif bx==cx;
            awopt[l,i]=a[na];
        else;
            awopt[l,i]=golden(&value2,ax,bx,cx,tol1);
        endif;

        vw[l,i]=value2(awopt[l,i]);
        y=varput(awopt[l,i],"kx");
        x0=0.3|0.3;
        {xf,jcode}=FixVMN1(x0,&rf);
        cwopt[l,i]=xf[1];
        nwopt[l,i]=xf[2];
    "q~i~l~k:" q~i~l~kbar;
    endo;   /* l */
endo;   /* i */

/* computation of the aggregate capital stock and employment nbar */
kgen=zeros(t+tr,1);
ngen=zeros(t,1);
cgen=zeros(t+tr,1);
kgen[1]=kinit;
vexact=0;
j=0;
do until j==t+tr-1;
j=j+1;
if j<t+1;
kgen[j+1]=kw(kgen[j],j);
    y=varput(kgen[j+1],"kx");
    y=varput(kgen[j],"k0");  
    x0=0.3|0.3;
    {xf,jcode}=FixVMN1(x0,&rf);  
    cgen[j]=xf[1];
    ngen[j]=xf[2];
    vexact=vexact+b^(j-1)*u(cgen[j],1-ngen[j]);
else;
    kgen[j+1]=kr(kgen[j],j-t);
    cgen[j]=(1+r)*kgen[j]+pen-kgen[j+1];
    vexact=vexact+b^(j-1)*u(cgen[j],1);
endif;
endo;
cgen[t+tr]=(1+r)*kgen[t+tr]+pen;
vexact=vexact+b^(t+tr-1)*u(cgen[t+tr],1);

knew=meanc(kgen);
kbar=phi*kold+(1-phi)*knew;
nnew=meanc(ngen)*2/3;
nbar=phi*nold+(1-phi)*nnew;
"nbar~kbar: " nbar~kbar;
vexact1=zeros(na,1);
j=0;
do until j==na;
j=j+1;
vexact1[j]=wvalue(a[j],1);
endo;
save vexact1,a;
save cv=cgen;
save kv=kgen;

endo;   /* q */


    "time elapsed: " etstr(hsec-bsec);
wait;
xlabel("generation");
ylabel("capital holdings");
xy(seqa(1,1,t+tr),kgen);
ylabel("consumption");
xy(seqa(1,1,t+tr),cgen);
ylabel("labor supply");
xy(seqa(1,1,t),ngen[1:t]);


"capital stock: " kbar;
"exact value for value function: " vexact;
"value function: " wvalue(kinit,1);
"tau: " tau;
"replacement ratio: " rep;
"pensions: " pen;
"r: " r;
"aggregate employment " nbar;
"mean consumption: " meanc(cgen);

vexact1=zeros(na,1);
j=0;
do until j==na;
j=j+1;
vexact1[j]=wvalue(a[j],1);
endo;
save vexact1,a;
save cv=cgen;
save kv=kgen;

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

proc rvalue(k,t0);
    local k0,n1,n2,rv;
    /* at the border of the interval */
        k0=k/kmax*(na-1)+1;
        n2=floor(k0);
        n1=k0-n2;
        if k<=0; retp(vr[1,t0]-(1-k0)*(vr[2,t0]-vr[1,t0])); endif;
        if k>=kmax; retp(vr[na,t0]); endif;

    if _VI_method==1;   /* linear interpolation */
        retp((1-n1)*vr[n2,t0]+n1*vr[n2+1,t0]);
    elseif _VI_method==2;   /* cubic spline */ 
        rv=Splint(a,vr0,y2,1,k);
        retp(rv);
    endif;
endp;

proc wvalue(k,t0);
    local k0,n1,n2,wv;
    k0=k/kmax*(na-1)+1;
    n2=floor(k0);
    n1=k0-n2;
    if k<=0; 
        if t0<=t;
            retp(vw[1,t0]-(1-k0)*(vw[2,t0]-vw[1,t0])); 
        else;
            retp(vr[1,1]-(1-k0)*(vr[2,1]-vr[1,1])) ; 
        endif;
    endif;
    if k>=kmax; 
        if t0<=t;
            retp(vw[na,t0]); 
        else;
            retp(vr[na,1]);
        endif;
    endif;
    if _VI_method==1;   /* linear interpolation */
        retp((1-n1)*vw[n2,t0]+n1*vw[n2+1,t0]);
    elseif _VI_method==2;   /* cubic spline */ 
        wv=Splint(a,vw0,y2,1,k);
        retp(wv);
    endif;
endp;


proc value1(x);
    local c;
    c=(1+r)*a[l]+pen-x;
    if c<=0; retp(-1e10); endif;
    retp(u(c,1)+b*rvalue(x,i+1));
endp;

proc value2(x);
    local y,xf,fvp,jk,jcode,c,empl;
    y=varput(x,"kx");   
    
    {xf,jcode}=FixVMN1(0.3|0.3,&rf);
    if jcode[1]>1; "jcode>1"; wait; endif;
    c=xf[1];
    empl=xf[2];
    if c<=0; retp(-1e10); endif;
    if _VI_method==1;
        if i==t;
            retp(u(c,1-empl)+b*rvalue(x,1));
        else;
            retp(u(c,1-empl)+b*wvalue(x,i+1));
        endif;
    elseif _VI_method==2;
        retp(u(c,1-empl)+b*wvalue(x,i+1));
    endif;
endp;


proc rf(x);
    local n,c,rf1,rf2,kx,k0;
    c=x[1];
    n=x[2];
    kx=varget("kx");
    k0=varget("k0");
    rf1=(1+r)*k0+(1-tau)*w*n-kx-c;
    rf2=un(c,1-n)./uc(c,1-n)-w*(1-tau);
    retp(rf1|rf2);
endp;

/* Procedure zur Berechnung vom Minimum - Golden Section Search
   Inputs:  &f  Funktion, fr die das Minimum bestimmt werden soll
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

proc kw(k,y);
local k0,n2,n1;
    k0=k/kmax*(na-1)+1;
    n2=floor(k0);
    n1=k0-n2;
    retp((1-n1)*awopt[n2,y]+n1*awopt[n2+1,y]);
endp;


proc kr(k,y);
local k0,n2,n1;
if y==tr;
    retp(0);
else;
    k0=k/kmax*(na-1)+1;
    n2=floor(k0);
    n1=k0-n2;
    retp((1-n1)*aropt[n2,y]+n1*aropt[n2+1,y]);
endif;
endp;


