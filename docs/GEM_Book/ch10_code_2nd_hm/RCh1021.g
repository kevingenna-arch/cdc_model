@ ------------------------------ Rch1021.g ----------------------------

    25.5.2008
    Burkhard Heer

    direct computation of the OLG model in section 10.2

-------------------------------------------------------------------------@


new;
clear all;
cls;
library pgraph;
Macheps=1e-30;
#include ch8_toolbox.src;
graphset;


@ ---------------------------

 parameter

---------------------------- @

beta=0.99;         /* discount factor */
s=2;            /* coefficient of relative risk aversion */
alpha=0.3;        /* production elasticity of capital */
rep=0.3;        /* replacement ratio */
delta=0.04;          /* rate of depreciation */
gam=2;            /* disutility from working */

tau=rep/(2+rep);    /* income tax rate */
t=40;
tr=20;
kmax=10;        /* upper limit of capital grid */
kinit=0;
na=101;          /* number of grid points on assets */
a=seqa(0,kmax/(na-1),na);   /* asset grid */
psi=0.0;      /* parameter of utility function */
phi=0.8;
tol=0.001;       /* percentage deviation of final solution */
tolk=0.001;       /* percentage deviation of final solution for k_1 */
tol1=1e-10;        /* tolerance for golden section search */
neg=-1e10;        /* initial value for value function */
nq1=30;

/* AR(1) process of technology */
rho0=0.95;  /* quarterly process */
sigma0=0.00763;
rho=0.95^4;
sigma=sqrt(1+rho0^2+rho0^4+rho0^6)*sigma0;

/* number of times for the simulation of the moments */
nsim=100;

@ ---------------------------------

initialization

--------------------------------- @

r=0.04;
nbar=0.2;
kbar=(alpha/(r+delta))^(1/(1-alpha))*nbar;
kold=100;
nold=2;


/* agents' policy function  */
aopt=zeros(t+tr,1);  /* optimal asset */
copt=zeros(t+tr,1);   /* optimal consumption */
nopt=0.3*ones(t,1);    /* optimal labor supply */
lambs=zeros(t+tr,1);	/* lambda */

@ ---------------------------------------------

iteration of policy function, wealth distribution,..

---------------------------------------------- @

q=0;
do until q==30 or abs((kbar-kold)/kbar)<tol;
    krit=abs((kbar-kold)/kbar);
    krit0=abs((nbar-nold)/nbar);
    q=q+1;
    w=(1-alpha)*kbar^alpha*nbar^(-alpha);
    r=alpha*kbar^(alpha-1)*nbar^(1-alpha)-delta;
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
            k60=0.4;
        elseif q1==2;
            k60=0.5;
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
            i; x0=aopt[i+1]|nopt[i];
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

xlabel("generation");
title("capital holdings");
xy(seqa(1,1,t+tr),aopt); wait;
title("labor supply");
xy(seqa(1,1,t),nopt);


i=0;
do until i==40;
	i=i+1;
	copt[i]=(1+r)*aopt[i]+(1-tau)*w*nopt[i]-aopt[i+1];
	lambs[i]=uc(copt[i],1-nopt[i]);
endo;
i=40;
do until i==59;
	i=i+1;
	copt[i]=(1+r)*aopt[i]+pen-aopt[i+1];
	lambs[i]=uc(copt[i],1);
endo;
copt[60]=(1+r)*aopt[60]+pen;
lambs[60]=uc(copt[60],1);
rbar=r;

"capital stock: " kbar;
"tau: " tau;
"replacement ratio: " rep;
"pensions: " pen;
"r: " rbar;
"aggregate employment " nbar;

cs=copt;
ks=aopt;
nst=nopt;
wbar=w;
cbar=meanc(copt);
ybar=kbar^alpha*nbar^(1-alpha);
ibar=delta*kbar;


@ ----- dynamics -------------------------------------------------- @

nc=103;   /* number of control variables */
nsc=59;  /* number of constraints = costate variable lambda */
ns=59;   /* number of state variables */
nex=1;  /* number of exogenous variables */



@ --------------  Definition of the M matrices  ------------------ @
/* states: k */
/* costates: lambda */
/* controls: c,n,r,w,tau */
/* exogenous: z*/


mcc=zeros(103,103);
mcs=zeros(103,118);
mce=zeros(103,1);

/* first-order condition labor supply */
i=0;
do until i==40;
	i=i+1;	
	mcc[i,i]=1-s;
	mcc[i,i+60]=(1-gam*(1-s))*nst[i]/(1-nst[i]);
	mcc[i,102]=-1;
	mcc[i,103]=tau/(1-tau);
	mcs[i,59+i]=1;
endo;

/* lambda=u_c */
j=0;
do until j==40;
	j=j+1;
	mcc[j+40,j]=-s;
	mcc[j+40,j+60]=-gam*(1-s)*nst[j]/(1-nst[j]);
	mcs[j+40,j+59]=1;
endo;

j=40;
do until j==59;
	j=j+1;
	mcc[j+40,j]=-s;
	mcs[j+40,j+59]=1;
endo;

/* consumption at s=60 */
mcc[100,60]=cs[60];
mcc[100,101]=-rbar*ks[60];
mcs[100,59]=(1+rbar)*ks[60];

/* w,r,tau */
j=0;
do until j==40;
	j=j+1;
	mcc[101,60+j]=alpha*nst[j]/nbar*1/60;
	mcc[102,60+j]=-(1-alpha)*nst[j]/nbar*1/60;
	mcc[103,60+j]=nst[j]/nbar*1/60;
endo;
mcc[101,102]=1;
mcc[102,101]=1;
mcc[103,102]=1;
mcc[103,103]=1;

j=0;
do until j==59;
	j=j+1;
	mcs[101,j]=alpha*ks[j+1]/kbar*1/60;
	mcs[101,j]=-(1-alpha)*ks[j+1]/kbar*1/60;
endo;

mce[101]=1;
mce[102]=1;


@ -----------------------------------------

dynamic equation in E_t

------------------------------------------ @

mss0=zeros(118,118);
mss1=zeros(118,118);
msc0=zeros(118,103);
msc1=zeros(118,103);
mse0=zeros(118,1);
mse1=zeros(118,1);


/* first order condition with respect to k^s_t+1 */
i=0;
do until i==58;
	i=i+1;
	mss0[i,59+1+i]=-1;	/* E_t lambda^s+1_t+1 */
	mss1[i,59+i]=1;		/* lambda^s_t */
	msc0[i,101]=rbar/(1+rbar);
endo;

/* foc for s=59 */
mss1[59,118]=1;
msc0[59,60]=-s;
msc0[59,101]=rbar/(1+rbar);

/* budget constraint */
j=0;
do until j==40;
	j=j+1;
	mss0[59+j,j]=ks[j+1];
	mss1[60+j,j]=-(1+rbar)*ks[j+1];
	msc1[59+j,j]=-cs[j];
	msc1[59+j,60+j]=(1-tau)*wbar*nst[j];
	msc1[59+j,101]=rbar*ks[j];
	msc1[59+j,102]=(1-tau)*wbar*nst[j];
	msc1[59+j,103]=-tau*wbar*nst[j];		
endo;


/* budget constraint */
j=40;
do until j==58;
	j=j+1;
	mss0[59+j,j]=ks[j+1];
	mss1[60+j,j]=-(1+rbar)*ks[j+1];
	msc1[59+j,j]=-cs[j];	
	msc1[59+j,101]=rbar*ks[j];
endo;

	j=59;
	mss0[59+j,j]=ks[j+1];
	msc1[59+j,j]=-cs[j];	
	msc1[59+j,101]=rbar*ks[j];


nc=103;   /* number of control variables */
nl=59;  /* number of constraints = costate variable  */
ns=59;   /* number of state variables */
nn=1;  /* number of exogenous variables */
nfx=4; 
/*  K, N, C, Y*/
fvv=zeros(nfx,ns+nn);
fvv[1,1:59]=ks[2:60]'/kbar*1/60;
fvv[4,1:59]=alpha*ks[2:60]'/kbar*1/60;
fvv[4,60]=1;


fvu=zeros(nfx,nc);
fvu[2,61:100]=nst'/nbar*1/60;
fvu[3,1:60]=cs'/cbar*1/60;
fvu[4,61:100]=(1-alpha)*nst'/nbar*1/60;

fvl=zeros(nfx,nl);


w=-pinv(mss0-msc0*inv(mcc)*mcs)*(mss1-msc1*inv(mcc)*mcs);
rr=pinv(mss0-msc0*inv(mcc)*mcs)*(mse0+msc0*inv(mcc)*mce);
q=pinv(mss0-msc0*inv(mcc)*mcs)*(mse1+msc1*inv(mcc)*mce);

{lamb0,p0}=eigv(w);

"number of eigenvalues<1";
z=abs(lamb0);
sumc(z.<1.0);
wait;

{hmat,mvmat}=kpr(mcc,mcs,mce,msc0,msc1,mss0,mss1,mse0,mse1,fvu,fvv,fvl,rho);
periods=25;
shocknr=1;
{states,exogs,controls,costates}=kprirfs(hmat,mvmat,periods,shocknr,ns,nn,nc,nl,nfx);
"maximum coefficient of imaginary part: ";
"hmat: " maxc(maxc(imag(hmat)));
"mvmat: " maxc(maxc(imag(mvmat)));
wait;
hmat=real(hmat);
mvmat=real(mvmat);

ct=controls[.,1];
i=1;
do until i==60;
    i=i+1;
    ct=ct+controls[.,i];
endo;
ct=ct/60;


nt=controls[.,61];
i=1;
do until i==40;
    i=i+1;
    nt=nt+controls[.,60+i];
endo;
nt=nt/40;

rt=controls[.,101];
wt=controls[.,102];
taut=controls[.,103];
kt=states[.,1];

i=1;
do until i==58;
    i=i+1;
    kt=kt+states[.,i+1];
endo;
kt=kt/60;

yt=exogs+alpha*kt+(1-alpha)*nt;
invest=yt*ybar/(delta*kbar)-ct*cbar/(delta*kbar);

select1=real(exogs~yt~invest~kt~ct~nt);


varnames="Z"|"Y"|"I"|"K"|"C"|"N";
graphset;
/*
_pmcolor=4;_pcolor=4|1|2;
*/
_pdate="";
           begwind;
             window(2,3,0);
             setwind(1);
             i=1;
             do while i<=cols(select1);
                xlabel(varnames[i]);
            @    ylabel("percent deviation");@
                _paxht=0.35;
                _pnumht=0.35;
            /*
				xtics(0,18,2,0);
            */
                @
				if i==6;
						_plegctl=1|5|1|0.63;
                      	_plegstr="gamma=0.1\00gamma=1\00gamma=2";
				endif;
                @
                xy(seqa(0,1,rows(select1)),select1[.,i]~zeros(rows(select1),1));
                i=i+1;
                nextwind;
             endo;
            endwind;

wait;

@ ------------------------------------------------------------

time series simulation 

------------------------------------------------------------ @

yresults=zeros(nsim,3);
iresults=zeros(nsim,3);
cresults=zeros(nsim,3);
nresults=zeros(nsim,3);
wresults=zeros(nsim,3);

isim=0;
do until isim==nsim;
    isim=isim+1;

tx=zeros(1001,60);  /* state variable */
tu=zeros(1000,nc+nl+nfx); /* controls */
ty=zeros(1000,1);   /* production */
tk=zeros(1001,1);   /* capital stock */
tks=zeros(1001,60);
tr=zeros(1001,1);   /* real interest rate */
tc=zeros(1000,1);  /* consumption in t */
tcs=zeros(1000,60);
tn=zeros(1000,1);   /* employment */
tns=zeros(1000,40);
tinvest=zeros(1000,1);   /* investment */
tw=zeros(1001,1);       /* wage */
shocks=rndn(1001,1);
tz=zeros(1001,1);       /* technology level */


/* initial values for pre-determined state variables */
tk[1]=kbar;
tks[1,.]=ks';
tw[1]=wbar;
tz[1]=0; /* z[0]==0 */
tr[1]=rbar;
tn[1]=nbar;
tns[1,.]=nst';
tc[1]=cbar;
tcs[1,.]=cs';
tinvest[1]=delta*kbar;
tx[1,.]=zeros(1,60);
xbar=ks|0;

i=0;
do until i==1000;
	i=i+1;
    /* percentage deviations from steady state */
	x=mvmat*tx[i,.]';
    tx[i+1,.]=x';
	tx[i+1,60]=tx[i+1,60]+shocks[i+1]*sigma;
	u=hmat*tx[i,.]';
    tu[i,.]=u';
    ty[i]=ybar*(1+tu[i,nc+nl+4]);
    tn[i]=nbar*(1+tu[i,nc+nl+2]);
    tc[i]=cbar*(1+tu[i,nc+nl+3]);
	tk[i]=kbar*(1+tu[i,nc+nn+1]);
    tr[i]=rbar*(1+tu[i,101]);
	tinvest[i]=ty[i]-tc[i];
	tw[i]=wbar*(1+tu[i,102]);
	i~ty[i]~isim;
endo;

@
graphset;
t0=seqa(1,1,1000);
begwind;
window(2,2,0);
setwind(1);
title("capital k");
/* _protate=1; */
_ptitlht=0.28;
_pnumht=0.22;
xy(t0,tk[1:1000]);
nextwind;
title("employment");
xy(t0,tn[1:1000]);
nextwind;
title("output");
xy(t0,ty[1:1000]);
nextwind;
title("technology");
xy(t0,tx[1:1000,60]);
endwind;
wait;
graphset;

"standard deviation technology: ";
stdc(tx[1:1000,60]);
wait;

title("technology");
xy(t0,tx[1:1000,60]);
wait;
@
if round(isim/20)==isim/20; cls; endif;
"Simulation: " i;
"HP-Filtered Momente: ";
{x0,yt}=hptrend(ln(ty[200:1000]),100);
sy=stdc(yt);
"s_y= " sy; 

{x0,ct}=hptrend(ln(tc[200:1000]),100);
sc=stdc(ct);
{x0,invest}=hptrend(ln(tinvest[200:1000]),100);
sinvest=stdc(invest);
{x0,nt}=hptrend(ln(tn[200:1000]),100);
sn=stdc(nt);
{x0,rt}=hptrend(ln(tr[200:1000]),100);
sr=stdc(rt);
{x0,wt}=hptrend(ln(tw[200:1000]),100);
sw=stdc(wt);

yresults[isim,1]=sy;
cresults[isim,1]=sc;
iresults[isim,1]=sinvest;
nresults[isim,1]=sn;
wresults[isim,1]=sw;


"std korrelationen std(x)/std(y) ";
"c: " sc~sc/sy; ccor=corrx(yt~ct); ccor;
"inv: " sinvest~sinvest/sy; icor=corrx(invest~yt); icor;
"n: " sn~sn/sy; ncor=corrx(nt~yt); ncor;
"w: " sw~sw/sy; wcor=corrx(wt~yt); wcor;
"r: " sr~sr/sy; corrx(rt~yt); 

yresults[isim,2]=1;
cresults[isim,2]=ccor[1,2];
iresults[isim,2]=icor[1,2];
nresults[isim,2]=ncor[1,2];
wresults[isim,2]=wcor[1,2];


"first-order autocorrelation: ";
"y: "; yar=corrx(yt[1:799]~yt[2:800]); yar;
"invest: "; iar=corrx(invest[1:799]~invest[2:800]); iar;
"c: "; car=corrx(ct[1:799]~ct[2:800]); car;
"n: "; nar=corrx(nt[1:799]~nt[2:800]); nar;
"w: "; war=corrx(wt[1:799]~wt[2:800]); war;


yresults[isim,3]=yar[1,2];
cresults[isim,3]=car[1,2];
iresults[isim,3]=iar[1,2];
nresults[isim,3]=nar[1,2];
wresults[isim,3]=war[1,2];

endo;

"Table 10.2: ";
meanc(yresults);
meanc(iresults);
meanc(cresults);
meanc(nresults);
meanc(wresults);


@ -------------------------------------------------------------------------------

procedures 

----------------------------------------------------------------------------------- @

proc(2)=kpr(muu,mus,mue,mu0,mu1,ms0,ms1,me0,me1,fvu,fvv,fvl,rho);
local nl,nt,ns,qus,que,qus1,qus2,msss0,msss1,msse0,msse1,w,a,lv,mvmat,uv,fv,hmat;

nl=cols(fvl);
nt=cols(mus);
ns=nt-nl;
qus=inv(muu)*mus;
que=inv(muu)*mue;
qus1=qus[.,1:ns];
qus2=qus[.,ns+1:nt];
msss0=ms0-mu0*qus;
msss1=ms1-mu1*qus;
msse0=me0+mu0*que;
msse1=me1+mu1*que;
w=-inv(msss0)*msss1;
a=inv(msss0)*(msse0*rho+msse1);
{lv,mvmat}=bk(w,a,rho,ns);
uv=(qus1~que)+qus2*lv; 
fv=fvu*uv+fvv+fvl*lv;
hmat=uv|fv|lv;
retp(hmat,mvmat);
endp;

proc(2)=bk(w,a,rho,ns);
local nt,nl,nn,w11,w12,evec,eval,mu1,ind,mu,p,mu2,ps,b,d1,d2,lee,ls,lemat,lv,mvmat;

nt=rows(w);
nl=nt-ns;
nn=rows(rho);
w11=w[1:ns,1:ns];
w12=w[1:ns,ns+1:nt];
{eval,evec}=eigv(w);
mu1=sortc(abs(eval),1);
ind=sortind(abs(eval));
mu=eval[ind];

print;
"Eigenvalues and absolute eigenvalues:";
mu~mu1;
print;
"Note: number of endogenous states is";;ns;
"Draw your own conclusions, then press key to continue...";
wait;

p=evec[.,ind];
mu2=diagrv(zeros(nt-ns,nt-ns),mu[ns+1:nt]);
ps=inv(p)*eye(rows(p));
b=ps*a;
d1=inv(mu2)*eye(rows(mu2));
d2=d1*b[ns+1:nt,.];
lee=-inv(eye(nl*nn)-rho'.*.d1)*vec(d2);
lee=reshape(lee,nl,nn);
ls=-inv(ps[ns+1:nt,ns+1:nt])*ps[ns+1:nt,1:ns];
lemat=inv(ps[ns+1:nt,ns+1:nt])*lee;
lv=ls~lemat;
mvmat=(w11+w12*ls)~(w12*lemat+a[1:ns,.])|(zeros(nn,ns)~rho);
retp(lv,mvmat);
endp;

proc(4)=kprirfs(hmat,mvmat,periods,shocknr,ns,nn,nc,nl,nfx);
local t,vt,shock,ut,yt,thisvt1,thisyt,controls,exogs,states,costates;
t=1;
vt=zeros(periods,ns+nn);
shock=zeros(periods,nn);
shock[2,shocknr]=1;
ut=zeros(periods,ns)~shock;
yt=zeros(periods,nc+nfx+nl);
do while t<=periods-1;
   thisvt1=mvmat*vt[t,.]'+ut[t,.]';
   vt[t+1,.]=thisvt1';
   thisyt=hmat*vt[t,.]';
   yt[t+1,.]=thisyt';
   t=t+1;
endo;
vt=real(vt[2:periods-1,.]);
yt=real(yt[2:periods-1,.]);
/*
in vt: endog states | exog states
in yt: controls | others (minimum one line) | costates
*/
states=vt[.,1:ns];
exogs=vt[.,ns+1:cols(vt)];
controls=yt[.,1:nc];
costates=yt[.,nc+nfx+1:cols(yt)];
retp(states,exogs,controls,costates);
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



proc (2)=HPTREND(x,mu);
   local m, t, i, j, a;
   /* Anlegen der TxT-Matrix m */
   t=rows(x);
   m=zeros(t,t);

   /* Belegen der Matrix */

   m[1,1]=1+mu;
   m[1,2]=-2*mu;
   m[1,3]=mu;
   m[2,1]=-2*mu;
   m[2,2]=1+5*mu;
   m[2,3]=-4*mu;
   m[2,4]=mu;
   m[3,1]=mu;
   m[3,2]=-4*mu;
   m[3,3]=1+6*mu;
   m[3,4]=-4*mu;
   m[3,5]=mu;

   i=4;
   do while i<=t-2;
    j=1;
    do while j<=5;
     m[i,i-3+j]=m[i-1,i-4+j];
     j=j+1;
    endo;
    i=i+1;
   endo;

   m[t-1,t-3]=mu;
   m[t-1,t-2]=-4*mu;
   m[t-1,t-1]=1+5*mu;
   m[t-1,t]=-2*mu;
   m[t,t-2]=mu;
   m[t,t-1]=-2*mu;
   m[t,t]=1+mu;
   m=invpd(m);
   a=m*x;
  retp(a,x-a);
endp;


