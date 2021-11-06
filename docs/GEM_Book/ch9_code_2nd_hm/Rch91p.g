@ ------------------------------ Rch91p.g ----------------------------

    18.4.2008

    author: Burkhard Heer

    60-period OLG computed with PROJECTION METHOD for
    the consumption and labor supply function, but not the value
    function

    COLLOCATION

    t years worker (employed or unemployed)
    tr years retired


-------------------------------------------------------------------------@


new;
clear all;
cls;
library pgraph;
Macheps=1e-30;
#include ch8_toolbox.src;
graphset;

bsec=hsec;


@ ---------------------------

 parameter

---------------------------- @

b=0.99;         /* discount factor */
r=0.01;         /* initial value of the interest rate */
s=2;            /* coefficient of relative risk aversion */
alp=0.3;        /* production elasticity of capital */
rep=0.3;        /* replacement ratio */
del=0.1;          /* rate of depreciation */
j=60;           /* total time periods */
tr=20;          /* retired */
t=40;           /* working time */
tau=rep/(2+rep);    /* income tax rate */
gam=2;            /* disutility from working */
kinit=0;

ncheb=3;        /* number of chebychev polynomials for value function */
kmax=5;        /* upper limit of capital grid */
kmin=0;
knode=(cos((2*seqa(1,1,ncheb+1)-1)/(2*(ncheb+1))*pi)+1)*(kmax-kmin)/2+kmin;

psi=0.001;      /* parameter of utility function */
phi=0.5;
tol=0.0001;       /* percentage deviation of final solution */

@ ---------------------------------

initialization

--------------------------------- @

nbar=0.3;
kbar=(alp/(r+del))^(1/(1-alp))*nbar;
kold=100;
nold=2;

@ ---------------------------------------------

iteration of policy function, wealth distribution,..

---------------------------------------------- @

q=0;
do until q==50 or abs((kbar-kold)/kbar)<tol;
krit=abs((kbar-kold)/kbar);
krit1=abs((nbar-nold)/nbar);
"q= " q;
q=q+1;
w=(1-alp)*kbar^alp*nbar^(-alp);
r=alp*kbar^(alp-1)*nbar^(1-alp)-del;
pen=rep*(1-tau)*w*nbar*(t+tr)/t;
"r= " r;
"pensions: " pen;
"k= " kbar;
"krit: " krit~krit1;
"nbar: " nbar;
kold=kbar;
nold=nbar;

/* chebyshev coefficients */
ac=zeros(ncheb+1,t+tr);         /* conumption function */
an=zeros(ncheb+1,t);         /* labor supply */

knode=(cos((2*seqa(1,1,ncheb+1)-1)/(2*(ncheb+1))*pi)+1)*(kmax-kmin)/2+kmin;

/* value and policy function in period t+tr */
ac[.,t+tr]=chebcoef(&c0,ncheb,2*ncheb,kmax|kmin);
an[.,t]=chebcoef(&n0,ncheb,2*ncheb,kmax|kmin);

i=t+tr;
do until i==1;
i=i-1;
/* solve for the first-order conditions */
if i>t;
    x0=ac[.,i+1];
    {xf,jcode}=FixVMN1(x0,&rf);
    ac[.,i]=xf[1:ncheb+1];
else;
    x0=ac[.,i+1]|an[.,t];
    {xf,jcode}=FixVMN1(x0,&rf0);
    ac[.,i]=xf[1:ncheb+1];
    an[.,i]=xf[ncheb+2:2*ncheb+2];
endif;
endo;


/* computation of the aggregate capital stock and employmen */
kgen=zeros(t+tr,1);
cgen=zeros(t+tr,1);
ngen=zeros(t,1);
kgen[1]=kinit;
vexact=0;
i=1;
do until i==t+tr;
cgen[i]=chebeval(kgen[i],ac[.,i],kmax|kmin);
if i<41;
    ngen[i]=chebeval(kgen[i],an[.,i],kmax|kmin);
    y=(1-tau)*ngen[i]*w;
    vexact=vexact+b^(i-1)*u(cgen[i],1-ngen[i]);
else;
    y=pen;
    vexact=vexact+b^(i-1)*u(cgen[i],1);
endif;

kgen[i+1]=(1+r)*kgen[i]+y-chebeval(kgen[i],ac[.,i],kmax|kmin);
i=i+1;
endo;
cgen[t+tr]=chebeval(kgen[t+tr],ac[.,t+tr],kmax|kmin);
    vexact=vexact+b^(t+tr-1)*u(cgen[t+tr],1);


knew=meanc(kgen);
nnew=meanc(ngen)*t/(t+tr);
kbar=phi*kold+(1-phi)*knew;
nbar=phi*nold+(1-phi)*nnew;
/*
if q==10;
xlabel("generation");
ylabel("capital holdings");
xy(seqa(1,1,t+tr),kgen);
ylabel("employment");
xy(seqa(1,1,t),ngen);
endif;
*/
endo;

"computational time: ";
    ?etstr(hsec-bsec);  wait;

save cp=cgen;
save kp=kgen;
save ac,an;
save rp=r;
save wp=w;
save np=nbar;
save penp=pen;

xlabel("generation");
ylabel("capital holdings");
xy(seqa(1,1,t+tr),kgen);
ylabel("consumption");
xy(seqa(1,1,t+tr),cgen);
ylabel("employment");
xy(seqa(1,1,t),ngen);

"capital stock: " kbar;
"exact value for value function: " vexact;
"tau: " tau;
"replacement ratio: " rep;
"pensions: " pen;
"r: " r;
"aggregate employment " nbar;
"mean consumption: " meanc(cgen);

@ ------------------------------------------------------

Procedures:

u  -- utility function
c0 -- consumption in period t+tr
rf -- residual function for retired
rf1 -- residual function for employed
uc -- marginal utility
un -- marginal utility of employment
n0 -- labor supply
u0 -- old-age utility initialization

------------------------------------------------------  @

proc u(x,y);
if s==1;
    retp(ln(x+psi)+gam*ln(y));
else;
    retp((((x+psi).*y^gam)^(1-s)-1)/(1-s));
endif;
endp;

proc u0(x);
    retp(u(c0(x),1));
endp;

proc c0(x);
retp((1+r)*x+pen);
endp;

proc n0(x);
retp(0.3*(1+x/kmax));
endp;

proc uc(x,y);
    retp((x+psi)^(-s).*y^(gam*(1-s)));
endp;

proc un(x,y);
    retp(gam*(x+psi)^(1-s).*y^(gam*(1-s)-1));
endp;


proc rf(x);
    local ac0,c,k1,rf1,c1,y;
    ac0=x;
    c=chebeval(knode,ac0,kmax|kmin);
    /* next period's capital stock */
    k1=(1+r)*knode+pen-c;
    c1=chebeval(k1,ac[.,i+1],kmax|kmin);
    /* first order condition with respect to k_t+1 */
    rf1=uc(c,1)./uc(c1,1)-b*(1+r);
    retp(rf1);
endp;


proc rf0(x);
    local n,ac0,an0,c,k1,rf1,rf2,empl,empl1,c1,y;
    n=rows(x)/2;
    ac0=x[1:n];
    an0=x[n+1:2*n];
    c=chebeval(knode,ac0,kmax|kmin);
    empl=chebeval(knode,an0,kmax|kmin);
    /* next period's capital stock, employment, consumption */
    k1=(1+r)*knode+(1-tau)*w*empl-c;
    if i==t;
    empl1=zeros(n,1);
    else;
    empl1=chebeval(k1,an[.,i+1],kmax|kmin);
    endif;
    c1=chebeval(k1,ac[.,i+1],kmax|kmin);
    /* first order condition with respect to k_t+1 */
    rf1=uc(c,1-empl)./uc(c1,1-empl1)-b*(1+r);
    rf2=un(c,1-empl)./uc(c,1-empl)-w*(1-tau);
    retp(rf1|rf2);
endp;


proc chebcoef(&f,n,m,d);
  local f:proc,a,i,t0,t1,t,x,y,z;
  z=-cos((2.*seqa(1,1,m)-1)./(2*m).*pi);
  x=(z+1).*(d[1]-d[2])/2+d[2];
  y=f(x);
  t0=ones(m,1);
  t1=z;
  a=zeros(n+1,1);
  a[1]=sumc(y)/m;
  a[2]=(y'*t1)./(t1'*t1);
  i=2;
  do until i==n+1;
    i=i+1;
    t=2.*z.*t1-t0;
    a[i]=(y'*t)./(t'*t);
    t0=t1;
    t1=t;
  endo;
  retp(a);
endp;


proc chebeval(x,a,d);
  local i,n,p,t0,t1,t,y,z;
  n=rows(a);
  p=rows(x);
  z=2*(x-d[2])./(d[1]-d[2])-1;
  t0=ones(p,1);
  t1=z;
  y=a[1]+a[2].*z;
  i=2;
  do until i==n;
    i=i+1;
    t=2.*z.*t1-t0;
    y=y+a[i].*t;
    t0=t1;
    t1=t;
  endo;
  retp(y);
endp;

