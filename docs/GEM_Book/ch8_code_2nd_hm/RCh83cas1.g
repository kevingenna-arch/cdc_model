@ ------------------------------ RCh83cas1.g ----------------------------

Burkard Heer, March 22, 2007

computation of the dynamics 
in the heterogenous-agent neoclassical growth model
with value function iteration

algorithm: see krusell/smith, 1998, JPE

-   interpolation between kbar (aggregate capital stock):
    linear 
    lininter.g
    
-   interpolation between k (individual capital stock):
   linear
   lininter.g

-   interpolation between a and k (optimal next period capital stock)
    bilinear
    bilina.g
   
-   Maximization: golden section search method
    golden.g
    

input procedures:
golden -- routine for golden section search
lininter -- linear interpolation
bilina -- bilinear interpolation of policy function
bellman -- right-hand side of bellman equation
matroot (x,n) - nth root of a matrix x

-------------------------------------------------------------------------@


new;
clear all;
cls;
library pgraph; 
graphset; 
Macheps=1e-20;
#include ch8_toolbox.src;


@
load path="d:\\buch\\prog";
save path="d:\\buch\\prog";
@


h0=hsec;


@ ----------------------------------------------------------------

Step 1: Parameterization and Initialization

----------------------------------------------------------------- @



amin1=0;                 /* asset grid individual agents */
amax1=800;
na=50;
a=seqa(amin1,(amax1-amin1)/(na-1),na);

kmin=80;
kmax=400;       /* asset grid aggregate capital stock */
nk=3;
k=seqa(kmin,(kmax-kmin)/(nk-1),nk);

/* simulation parameters */
nt=2000; /* number of transition periods */
nt1=100; /* stationary region nt1<nt */
nh=5000;    /* number of households for simulation */
kind00=180;   /* initial capital stock of all households */

/* numerical computation */
tol=0.01;              /* stopping criterion for value function */
tol1=1e-7;            /* stopping criterion for golden section search */
eps=0.05;
phi1=0.9;   /* updating of gam0,gam1 */
nit=30;           /* number of iterations over value function */
nq=10;                  /* number of iterations over aggrgate capital stock dynamics */
kritgam=1;
critfoc=100;
tolgam=0.01;
alpha=0.3750;
beta=0.96^(1/8);
delta=1-0.9^(1/8);
sigma=1.5;
r=0.03;
w=0.5;
w0=0.05;
home=0.25;

/* calibration of z */
zg=1;
zb=0.9130;

/* working hours in good and bad times */
hz=0.32|0.3;

zg=zg/hz[1]^(1-alpha);
zb=zb/hz[2]^(1-alpha);
z0=zg;      /* initial technology level */


/* calibration of p in zg,zb */
ppz=(0.9722~0.0278|
    0.0278~0.9722);

ni=5;	/* number of productivities */

ef=(0.509|0.787|1.000|1.290|2.081);


@ ----------------------------------------------------------------------------

calibration as in Castaneda et al 

----------------------------------------------------------------  @

nz=(0.8612~0.8232|
	0.9246~0.8854|
	0.9376~0.9024|
	0.9399~0.9081|
	0.9375~0.9125);


un=1-nz;	/* unemployment in good
				and bad times of productivity type i */

/* calibration of transition matrix of employment */
/* from good to good times */
/* z=zg, z'=zg */

udur=0.4;
nrate=nz[1,1];
x0=0.95|0.05|0.95|0.05;


{xf,crit1}=FixVMN1(x0,&transp);
ppgg1=xf[1]~xf[2]|xf[3]~xf[4];

"ppgg1: ";
ppgg1;
"ergodic ppgg1: ";
equivec1(ppgg1);
wait;

nrate=nz[2,1];
x0=0.95|0.05|0.95|0.05;
{xf,crit1}=FixVMN1(x0,&transp);
ppgg2=xf[1]~xf[2]|xf[3]~xf[4];

nrate=nz[3,1];
x0=0.95|0.05|0.95|0.05;
{xf,crit1}=FixVMN1(x0,&transp);
ppgg3=xf[1]~xf[2]|xf[3]~xf[4];

nrate=nz[4,1];
x0=0.95|0.05|0.95|0.05;
{xf,crit1}=FixVMN1(x0,&transp);
ppgg4=xf[1]~xf[2]|xf[3]~xf[4];

nrate=nz[5,1];
x0=0.95|0.05|0.95|0.05;
{xf,crit1}=FixVMN1(x0,&transp);
ppgg5=xf[1]~xf[2]|xf[3]~xf[4];

save ppgg1,ppgg2,ppgg3,ppgg4,ppgg5;

"ppgg1,ppgg5:";
ppgg1;
ppgg5;
wait;



/* z=zb, z'=zb */

udur=1-6/14;
udur=0.57143;
nrate=nz[1,2];
x0=0.95|0.05|0.95|0.05;
{xf,crit1}=FixVMN1(x0,&transp);
ppbb1=xf[1]~xf[2]|xf[3]~xf[4];

nrate=nz[2,2];
x0=0.95|0.05|0.95|0.05;
{xf,crit1}=FixVMN1(x0,&transp);
ppbb2=xf[1]~xf[2]|xf[3]~xf[4];

nrate=nz[3,2];
x0=0.95|0.05|0.95|0.05;
{xf,crit1}=FixVMN1(x0,&transp);
ppbb3=xf[1]~xf[2]|xf[3]~xf[4];

nrate=nz[4,2];
x0=0.95|0.05|0.95|0.05;
{xf,crit1}=FixVMN1(x0,&transp);
ppbb4=xf[1]~xf[2]|xf[3]~xf[4];

nrate=nz[5,2];
x0=0.95|0.05|0.95|0.05;
{xf,crit1}=FixVMN1(x0,&transp);
ppbb5=xf[1]~xf[2]|xf[3]~xf[4];


"ppbb1,ppbb5:";
ppbb1;
ppbb5;

save ppbb1,ppbb2,ppbb3,ppbb4,ppbb5;

/* z=zg,z'=zb */
ppgb1=(nz[1,2]/nz[1,1]~1-nz[1,2]/nz[1,1]|0~1);
ppgb2=(nz[2,2]/nz[2,1]~1-nz[2,2]/nz[2,1]|0~1);
ppgb3=(nz[3,2]/nz[3,1]~1-nz[3,2]/nz[3,1]|0~1);
ppgb4=(nz[4,2]/nz[4,1]~1-nz[4,2]/nz[4,1]|0~1);
ppgb5=(nz[5,2]/nz[5,1]~1-nz[5,2]/nz[5,1]|0~1);

/* z=zb,z'=zg */
ppbg1=(1~0|1-(1-nz[1,1])/(1-nz[1,2])~(1-nz[1,1])/(1-nz[1,2]));
ppbg2=(1~0|1-(1-nz[2,1])/(1-nz[2,2])~(1-nz[2,1])/(1-nz[2,2]));
ppbg3=(1~0|1-(1-nz[3,1])/(1-nz[3,2])~(1-nz[3,1])/(1-nz[3,2]));
ppbg4=(1~0|1-(1-nz[4,1])/(1-nz[4,2])~(1-nz[4,1])/(1-nz[4,2]));
ppbg5=(1~0|1-(1-nz[5,1])/(1-nz[5,2])~(1-nz[5,1])/(1-nz[5,2]));


"Press any key to continue.. ";
wait;

/* transition matrix of efficiency type i */
/* as in castaneda et al. 1998 */
ppef=(	1.0~0~0~0~0|
		0~1.0~0~0~0|
		0~0~1.0~0~0|
		0~0~0~1.0~0|
		0~0~0~0~1.0);

/* yearly transition matrix of earnings as in 
diaz-giminez et al, 1997 */
@
ppef=( 	0.858~0.116~0.014~0.006~0.005|
		0.186~0.409~0.300~0.071~0.034|
		0.071~0.120~0.470~0.262~0.076|
		0.075~0.068~0.175~0.465~0.217|
		0.058~0.041~0.055~0.183~0.663);

@

ppef=(    0.58~0.28~0.09~0.03~0.02|
		0.22~0.44~0.22~0.08~0.03|
		0.10~0.15~0.43~0.23~0.09|
		0.06~0.09~0.18~0.46~0.21|
		0.06~0.02~0.06~0.21~0.65);
	
	ppef=ppef./sumc(ppef');

"yearly transtion matrix: ";
ppef;

ppef=matroot(ppef,20);
"6-weekly transition matrix: ";
ppef;
"press any key... ";

wait;


/* initialization of value functions */
/* ve - employed agent */
/* vu - unemployed agent */
ve1=u(r*a+ef[1]*w)/(1-beta);
vu1=u(r*a+w0)/(1-beta);
ve=ve1;
vu=vu1;

i=1;
do until i==nk;
    i=i+1;
    ve=ve~u(r*a+ef[1]*w)/(1-beta);
    vu=vu~vu1;
endo;

/* initialization of the value function: */
v1eg=ve;    /* productivity one, employed, good state */
v1eb=ve;    /* productivity one, employed, bad state */
v1ug=vu;    /* productivity one, unemployed, good state */
v1ub=vu;    /* productivity one, unemployed, bad state */
v2eg=ve;    /* productivity two, employed, good state */
v2eb=ve;    /* productivity two, employed, bad state */
v2ug=vu;    /* productivity two, unemployed, good state */
v2ub=vu;    /* productivity two, unemployed, bad state */
v3eg=ve;    
v3eb=ve;  
v3ug=vu;   
v3ub=vu;
v4eg=ve;    
v4eb=ve;  
v4ug=vu;   
v4ub=vu;
v5eg=ve;    
v5eb=ve;  
v5ug=vu;   
v5ub=vu;


copt1eg=zeros(na,nk);           /* optimal consumption, employed */
copt1eb=zeros(na,nk);
copt2eg=zeros(na,nk);
copt2eb=zeros(na,nk);
copt3eg=zeros(na,nk);
copt3eb=zeros(na,nk);
copt4eg=zeros(na,nk);
copt4eb=zeros(na,nk);
copt5eg=zeros(na,nk);
copt5eb=zeros(na,nk);
copt1ug=zeros(na,nk);           /* optimal consumption, unemployed*/
copt1ub=zeros(na,nk);
copt2ug=zeros(na,nk);
copt2ub=zeros(na,nk);
copt3ug=zeros(na,nk);
copt3ub=zeros(na,nk);
copt4ug=zeros(na,nk);
copt4ub=zeros(na,nk);
copt5ug=zeros(na,nk);
copt5ub=zeros(na,nk);

aopt1eg=zeros(na,nk);           /* optimal next-period assets, employed */
aopt1eb=zeros(na,nk);
aopt2eg=zeros(na,nk);       
aopt2eb=zeros(na,nk);
aopt3eg=zeros(na,nk);       
aopt3eb=zeros(na,nk);
aopt4eg=zeros(na,nk);       
aopt4eb=zeros(na,nk);
aopt5eg=zeros(na,nk);       
aopt5eb=zeros(na,nk);
aopt1ug=zeros(na,nk);           /* optimal next-period assets, unemployed */
aopt1ub=zeros(na,nk);
aopt2ug=zeros(na,nk);       
aopt2ub=zeros(na,nk);
aopt3ug=zeros(na,nk);       
aopt3ub=zeros(na,nk);
aopt4ug=zeros(na,nk);       
aopt4ub=zeros(na,nk);
aopt5ug=zeros(na,nk);       
aopt5ub=zeros(na,nk);

/* law of motion for capital: initial guess, 
ln(kbar')=gam0+gam1*ln(kbar) */
gam0g=0.03;
gam1g=0.99;
gam0b=0.01;
gam1b=0.99;

kmean=kind00;


gam0q=zeros(nq,2);
gam1q=zeros(nq,2);
kq=zeros(nq,1); 

/* gini coefficients of wealth and income */
giniw=zeros(nt,1);
giniy=zeros(nt,1);
ginie=zeros(nt,1);


@ ---------------------------------------------------------------------

Step 1: computation of the employment in good times and in bad times 

----------------------------------------------------------------------- @

/*employment in good and bad times */
nn0=ef'*nz[.,1]|ef'*nz[.,2];
nn0=nn0/5;
"nn0= "; nn0;
wait;

q=0;	/* iteration over capital dynamics */
do until q==nq or (q>5 and kritgam<tolgam);
    q=q+1;

    @ ----------------------------------------------------------------

    Step 2: Iteration of the value function over kbar and a

    ----------------------------------------------------------------- @

    crit=tol+1;
    neg=-1e10;
    j=0; 
    do until (crit<tol and j>10) or (j==nit);  
    j=j+1; aopt=0;
        "iteration: " j;
        "time elapsed: " etstr(hsec-h0);
        "error: " crit;

        vold1eg=v1eg;
        vold1eb=v1eb;
        vold1ug=v1ug;
        vold1ub=v1ub;
        vold2eg=v2eg;
        vold2eb=v2eb;
        vold2ug=v2ug;
        vold2ub=v2ub;
        vold3eg=v3eg;
        vold3eb=v3eb;
        vold3ug=v3ug;
        vold3ub=v3ub;
        vold4eg=v4eg;
        vold4eb=v4eb;
        vold4ug=v4ug;
        vold4ub=v4ub;
        vold5eg=v5eg;
        vold5eb=v5eb;
        vold5ug=v5ug;
        vold5ub=v5ub;

        /* home production = 25% of average wage income */
        w0=home*(1-alpha)*meanc(zg|zb)*kmean^alpha*meanc(hz.*nn0)^(-alpha)*meanc(hz);

        /* technology state of the economy */
        s=0;
        do until s==2;
            s=s+1;

            if s==1; 
                zz=zg; 
                gam0=gam0g; 
                gam1=gam1g;  
            else; 
                zz=zb; 
                gam0=gam0b; 
                gam1=gam1g;
            endif;

            /* iteration over k */
            m=0;
            do until m==nk;
                m=m+1;
                r=rrate(zz,k[m],nn0[s]*hz[s]); 
                w=wage(zz,k[m],nn0[s]*hz[s]); 
                k1=exp(gam0+gam1*ln(k[m]));  /* next period aggregate capital stock */
                v1eg1=zeros(na,1);    /* value function of employed with k'=k1, z'=zg, eps'=1 in t+1 */
                v1eb1=zeros(na,1);
                v1ug1=zeros(na,1);    /* value function of unemployed with k'=k1, z'=zg, eps'=1 */
                v1ub1=zeros(na,1);
                v2eg1=zeros(na,1);  
                v2eb1=zeros(na,1);
                v2ug1=zeros(na,1);   
                v2ub1=zeros(na,1);
                v3eg1=zeros(na,1);  
                v3eb1=zeros(na,1);
                v3ug1=zeros(na,1);   
                v3ub1=zeros(na,1);
                v4eg1=zeros(na,1);  
                v4eb1=zeros(na,1);
                v4ug1=zeros(na,1);   
                v4ub1=zeros(na,1);
                v5eg1=zeros(na,1);  
                v5eb1=zeros(na,1);
                v5ug1=zeros(na,1);   
                v5ub1=zeros(na,1);

                i=0; 
                do until i==na;
                    i=i+1;
                    if k1<=kmin; 
                        v1eg1[i]=v1eg[i,1];
                        v1eb1[i]=v1eb[i,1];
                        v1ug1[i]=v1ug[i,1];
                        v1ub1[i]=v1ub[i,1];
                        v2eg1[i]=v2eg[i,1];
                        v2eb1[i]=v2eb[i,1];
                        v2ug1[i]=v2ug[i,1];
                        v2ub1[i]=v2ub[i,1];
                        v3eg1[i]=v3eg[i,1];
                        v3eb1[i]=v3eb[i,1];
                        v3ug1[i]=v3ug[i,1];
                        v3ub1[i]=v3ub[i,1];
                        v4eg1[i]=v4eg[i,1];
                        v4eb1[i]=v4eb[i,1];
                        v4ug1[i]=v4ug[i,1];
                        v4ub1[i]=v4ub[i,1];
                        v5eg1[i]=v5eg[i,1];
                        v5eb1[i]=v5eb[i,1];
                        v5ug1[i]=v5ug[i,1];
                        v5ub1[i]=v5ub[i,1];
                    elseif k1>=kmax;
                        v1eg1[i]=v1eg[i,nk];
                        v1eb1[i]=v1eb[i,nk];
                        v1ug1[i]=v1ug[i,nk];
                        v1ub1[i]=v1ub[i,nk];
                        v2eg1[i]=v2eg[i,nk];
                        v2eb1[i]=v2eb[i,nk];
                        v2ug1[i]=v2ug[i,nk];
                        v2ub1[i]=v2ub[i,nk];
                        v3eg1[i]=v3eg[i,nk];
                        v3eb1[i]=v3eb[i,nk];
                        v3ug1[i]=v3ug[i,nk];
                        v3ub1[i]=v3ub[i,nk];
                        v4eg1[i]=v4eg[i,nk];
                        v4eb1[i]=v4eb[i,nk];
                        v4ug1[i]=v4ug[i,nk];
                        v4ub1[i]=v4ub[i,nk];
                        v5eg1[i]=v5eg[i,nk];
                        v5eb1[i]=v5eb[i,nk];
                        v5ug1[i]=v5ug[i,nk];
                        v5ub1[i]=v5ub[i,nk];
                    else;        
                        v1eg1[i]=lininter(k,v1eg[i,.]',k1);
                        v1eb1[i]=lininter(k,v1eb[i,.]',k1);
                        v1ug1[i]=lininter(k,v1ug[i,.]',k1);
                        v1ub1[i]=lininter(k,v1ub[i,.]',k1);
                        v2eg1[i]=lininter(k,v2eg[i,.]',k1);
                        v2eb1[i]=lininter(k,v2eb[i,.]',k1);
                        v2ug1[i]=lininter(k,v2ug[i,.]',k1);
                        v2ub1[i]=lininter(k,v2ub[i,.]',k1);
                        v3eg1[i]=lininter(k,v3eg[i,.]',k1);
                        v3eb1[i]=lininter(k,v3eb[i,.]',k1);
                        v3ug1[i]=lininter(k,v3ug[i,.]',k1);
                        v3ub1[i]=lininter(k,v3ub[i,.]',k1);
                        v4eg1[i]=lininter(k,v4eg[i,.]',k1);
                        v4eb1[i]=lininter(k,v4eb[i,.]',k1);
                        v4ug1[i]=lininter(k,v4ug[i,.]',k1);
                        v4ub1[i]=lininter(k,v4ub[i,.]',k1);
                        v5eg1[i]=lininter(k,v5eg[i,.]',k1);
                        v5eb1[i]=lininter(k,v5eb[i,.]',k1);
                        v5ug1[i]=lininter(k,v5ug[i,.]',k1);
                        v5ub1[i]=lininter(k,v5ub[i,.]',k1);
                    endif;
                endo;

                h=0;  /* iteration over productivity types eps=1,., ni */
                do until h==ni;	
                    h=h+1;
                    e=0;    /* iteration over the employment status */
                    do until e==2;  /* e=1 employed, e=2 unemployed */
                        e=e+1;
                        i=0;        /* iteration over asset grid a in period t */
                        l0=amin1;
                        do until i==na;
                            i=i+1;
                            v0=neg;
                            aopt=GSS(&value1,l0,amax1);
                            l0=aopt;

                            if e==1;
			                    if s==1;
       			                    if h==1;
				                        aopt1eg[i,m]=aopt;
				                        v1eg[i,m]=bellman(a[i],aopt,h);
			                        elseif h==2;
				                        aopt2eg[i,m]=aopt;
				                        v2eg[i,m]=bellman(a[i],aopt,h);
			                        elseif h==3;				
				                        aopt3eg[i,m]=aopt;
				                        v3eg[i,m]=bellman(a[i],aopt,h);
			                        elseif h==4;				
				                        aopt4eg[i,m]=aopt;
				                        v4eg[i,m]=bellman(a[i],aopt,h);
			                        elseif h==5;
				                        aopt5eg[i,m]=aopt;
				                        v5eg[i,m]=bellman(a[i],aopt,h);
			                        endif;
			                    else; /* s==2*/
       		       	                if h==1;
				                        aopt1eb[i,m]=aopt;
				                        v1eb[i,m]=bellman(a[i],aopt,h);
			                        elseif h==2;
				                        aopt2eb[i,m]=aopt;
				                        v2eb[i,m]=bellman(a[i],aopt,h);
			                        elseif h==3;				
				                        aopt3eb[i,m]=aopt;
				                        v3eb[i,m]=bellman(a[i],aopt,h);
			                        elseif h==4;				
				                        aopt4eb[i,m]=aopt;
				                        v4eb[i,m]=bellman(a[i],aopt,h);
			                        elseif h==5;
				                        aopt5eb[i,m]=aopt;
				                        v5eb[i,m]=bellman(a[i],aopt,h);
			                        endif; /* h */
 			                    endif; /* s */
                            else;	/* e==2 */
			                    if s==1;
        		                    if h==1;
				                        aopt1ug[i,m]=aopt;
				                        v1ug[i,m]=bellman(a[i],aopt,h);
			                        elseif h==2;
				                        aopt2ug[i,m]=aopt;
				                        v2ug[i,m]=bellman(a[i],aopt,h);
			                        elseif h==3;				
				                        aopt3ug[i,m]=aopt;
				                        v3ug[i,m]=bellman(a[i],aopt,h);
			                        elseif h==4;				
				                        aopt4ug[i,m]=aopt;
				                        v4ug[i,m]=bellman(a[i],aopt,h);
			                        elseif h==5;
				                        aopt5ug[i,m]=aopt;
				                        v5ug[i,m]=bellman(a[i],aopt,h);
			                        endif;
			                    else;
       		       	                if h==1;
				                        aopt1ub[i,m]=aopt;
				                        v1ub[i,m]=bellman(a[i],aopt,h);
			                        elseif h==2;
				                        aopt2ub[i,m]=aopt;
				                        v2ub[i,m]=bellman(a[i],aopt,h);
			                        elseif h==3;				
				                        aopt3ub[i,m]=aopt;
				                        v3ub[i,m]=bellman(a[i],aopt,h);
			                        elseif h==4;				
				                        aopt4ub[i,m]=aopt;
				                        v4ub[i,m]=bellman(a[i],aopt,h);
			                        elseif h==5;
				                        aopt5ub[i,m]=aopt;
				                        v5ub[i,m]=bellman(a[i],aopt,h);
			                        endif; /* h */
			                    endif;  /* s */
                            endif; /* e */
                        endo;   /* i=1,..na */

                            "run time: " etstr(hsec-h0);
                            "q~j~s~meanc(crit)~critfoc~kmean";    
                            q~j~s~meanc(crit)~critfoc~kmean;
                            "r~m~e~h: " r~m~e~h;
	                endo;  /* h=1,..ni */
                endo;   /* e=1,2 */
            endo;   /* m=1..nk */

        endo;   /* s=1,2 -- zz=zg,zb */


        if (q==nq) and (j==nit);
            title("aopt1eg-a");
            xy(a,aopt1eg-a);
            wait;
            title("aop2teg-a");
            xy(a,aopt2eg-a);
            wait;
            title("aopt3eg-a");
            xy(a,aopt3eg-a);
            wait;
            title("aopt4eg-a");
            xy(a,aopt4eg-a);
            wait;
            title("aopt5eg-a");
            xy(a,aopt5eg-a);
            wait;
            title("aopt1ug-a");
            xy(a,aopt1ug-a);
            title("aopt2ug-a");
            xy(a,aopt2ug-a);
            title("aopt3ug-a");
            xy(a,aopt3ug-a);
            title("aopt4ug-a");
            xy(a,aopt4ug-a);
            title("aopt5ug-a");
            xy(a,aopt5ug-a);
            title("aopt1eb-a");
            xy(a,aopt1eb-a);
            title("aopt2eb-a");
            xy(a,aopt2eb-a);
            title("aopt3eb-a");
            xy(a,aopt3eb-a);
            title("aopt4eb-a");
            xy(a,aopt4eb-a);
            title("aopt5eb-a");
            xy(a,aopt5eb-a);
            title("aopt1ub-a");
            xy(a,aopt1ub-a);
            title("aopt2ub-a");
            xy(a,aopt2ub-a);
            title("aopt3ub-a");
            xy(a,aopt3ub-a);
            title("aopt4ub-a");
            xy(a,aopt4ub-a);
            title("aopt5ub-a");
            xy(a,aopt5ub-a);
            graphset;
        endif;

        crit1=meanc(abs(vold1eg-v1eg))|meanc(abs(vold1eb-v1eb))|meanc(abs(vold1ug-v1ug))|meanc(abs(vold1ub-v1ub));
        crit2=meanc(abs(vold2eg-v2eg))|meanc(abs(vold2eb-v2eb))|meanc(abs(vold2ug-v2ug))|meanc(abs(vold2ub-v2ub));
        crit3=meanc(abs(vold3eg-v3eg))|meanc(abs(vold3eb-v3eb))|meanc(abs(vold3ug-v3ug))|meanc(abs(vold3ub-v3ub));
        crit4=meanc(abs(vold4eg-v4eg))|meanc(abs(vold4eb-v4eb))|meanc(abs(vold4ug-v4ug))|meanc(abs(vold4ub-v4ub));
        crit5=meanc(abs(vold5eg-v5eg))|meanc(abs(vold5eb-v5eb))|meanc(abs(vold5ug-v5ug))|meanc(abs(vold5ub-v5ub));
        crit=meanc(crit1|crit2|crit3|crit4|crit5);

        "iteration q~j " q~j;
        "time elapsed: " etstr(hsec-h0);
        "error value function: " crit;


    endo;   /* j=1,..nit */


                "test: violation of the FOC: ";
                critfoc=0; typ1=seqa(1,1,10);
                s=0;
                do until s==2;
                    s=s+1;
                    if s==1;
                        zz=zg;
                    else;
                        zz=zb;
                    endif;
                    e=0;
                    do until e==2;
                        e=e+1;
                        h=0;
                        do until h==5;
                            h=h+1;
                            i=(e-1)*5+h;
                            resid=zeros(na,nk);
                            ia=0;
                            do until ia==na;
                                ia=ia+1;
                                a0=a[ia];
                                ik=0;
                                do until ik==nk;
                                    ik=ik+1;
                                    k0=k[ik];
                                    resid[ia,ik]=rf2(a0,k0);   
                                endo;
                            endo;
                            critfoc=critfoc+meanc(meanc(resid));
                        endo;
                    endo;
                endo;
                critfoc=critfoc/20;


    save v1eg,v2eg,v3eg,v4eg,v5eg,aopt1eg,aopt2eg,aopt3eg,aopt4eg,aopt5eg;
    save v1ug,v2ug,v3ug,v4ug,v5ug,aopt1ug,aopt2ug,aopt3ug,aopt4ug,aopt5ug;
    save v1eb,v2eb,v3eb,v4eb,v5eb,aopt1eb,aopt2eb,aopt3eb,aopt4eb,aopt5eb;
    save v1ub,v2ub,v3ub,v4ub,v5ub,aopt1ub,aopt2ub,aopt3ub,aopt4ub,aopt5ub;

    @   ------------------------------------------------------------

    Step 5:     simulation of the dynamics

    -------------------------------------------------------------- @


    kt=zeros(nt,2); /* aggregate capital stock K and technology z over time */
    yt=zeros(nt-1,1);	/* output */
    incomet=zeros(nt-1,6);	/* income quintiles, 80-95%, top 5% */
    incom=zeros(nh,1);		/* income of individual in period t-1 */
    earnt=zeros(nt-1,6);
    earn=zeros(nh,1);
    kind0=zeros(nh,2);  /* individual capital stock k and employment e/efficiency h in t-1 */
    kind1=zeros(nh,2);  /* k,(e,i) in t */

    nnh=nh/5;


    /* typ[1]: number of employed agents, efficiency=1 */
    /* typ[2]: number of employed agents, efficiency=2 */
    /* typ[3]m,typ[4],typ[5]... */
    /* typ[6]: number of unemployed agents, efficiency = 1 */


    kind00=maxc(kmean|kmin);
    kind00=minc(kmax|kmean);
    kind0[1:nh,1]=kind00*ones(nh,1);

    if rndu(1,1)<0.5; z0=zg; else; z0=zb; endif;
    kt1=meanc(kind0[.,1]);
    cls;
    kt[1,.]=kt1~z0;

    t=1;
    do until t==nt; /* transition of distribution function */
 
        t=t+1; 
        q~t~kt[t-1,.]; 
        
        "kmean: " kmean;
        xx=maxc(1|q-10);
        "K: " kq[xx:q];   
        gam0g~gam0b~gam1g~gam1b;

        if kt[t-1,2]==zg; 
            typ=(1-un[.,1])*nnh|un[.,1]*nnh; 
            s0=1;
        else; 
	        s0=2; 
            typ=(1-un[.,2])*nnh|un[.,2]*nnh;
        endif;

	    typ1=cumsumc(typ);
	    typ1=round(typ1);
        kind0[.,2]=zeros(nh,1);
	    kind0[1:typ1[5],2]=ones(typ1[5],1);


        x0=rndu(1,1);  
        if x0<ppz[s0,1];
            kt[t,2]=zg;
            s=1;
        else;
            kt[t,2]=zb;
            s=2;
        endif;

 
        xs=rndu(nh,1);  /* employed/unemployed in t+1 */
        xs1=rndu(nh,1);  /* efficiency type in t+1 */
        xs2=rndu(nh,1);  	/* random order of households */


        k0=kt[t-1,1];   /* capital stock in t */
        k0=minc(k0|kmax);
        k0=maxc(k0|kmin);
        z0=kt[t-1,2];
        r=rrate(z0,k0,nn0[s0]*hz[s0]); /* interest rate in t */
        w=wage(z0,k0,nn0[s0]*hz[s0]);  /* wage in t */
        "r~w~w0: " r~w~w0;
        "critfoc: " critfoc;
        i=0;
        do until i==nh;
            i=i+1;
            kind1[i,1]=bilina(kind0[i,1],kt[t-1,1],i,kt[t-1,2]);
            if i<=typ1[1];
                ppgg=ppgg1[1,1];
                ppgb=ppgb1[1,1];
                ppbg=ppbg1[1,1];
                ppbb=ppbb1[1,1];                                        
                j0=sumc( cumsumc(ppef[1,.]').<=xs1[i] )+1; 
	            incom[i]=r*kind0[i,1]+ef[1]*w*hz[s0];
	            earn[i]=ef[1]*w*hz[s0];
            endif;  
            if typ1[1]<i and i<=typ1[2];
                ppgg=ppgg2[1,1];
                ppgb=ppgb2[1,1];
                ppbg=ppbg2[1,1];
                ppbb=ppbb2[1,1];
                j0=sumc(cumsumc(ppef[2,.]').<=xs1[i])+1;
	            incom[i]=r*kind0[i,1]+ef[2]*w*hz[s0];
	            earn[i]=ef[2]*w*hz[s0];
	        endif;
            if typ1[2]<i and i<=typ1[3];
                ppgg=ppgg3[1,1];
                ppgb=ppgb3[1,1];
                ppbg=ppbg3[1,1];
                ppbb=ppbb3[1,1];
                j0=sumc(cumsumc(ppef[3,.]').<=xs1[i])+1;      
	            incom[i]=r*kind0[i,1]+ef[3]*w*hz[s0];
  	            earn[i]=ef[3]*w*hz[s0];
	        endif;
            if typ1[3]<i and i<=typ1[4];
                ppgg=ppgg4[1,1];
                ppgb=ppgb4[1,1];
                ppbg=ppbg4[1,1];
                ppbb=ppbb4[1,1];
                j0=sumc(cumsumc(ppef[4,.]').<=xs1[i])+1;      
	            incom[i]=r*kind0[i,1]+ef[4]*w*hz[s0];
	            earn[i]=ef[4]*w*hz[s0];
	        endif;
            if typ1[4]<i and i<=typ1[5];
                ppgg=ppgg5[1,1];
                ppgb=ppgb5[1,1];
                ppbg=ppbg5[1,1];
                ppbb=ppbb5[1,1];
                j0=sumc(cumsumc(ppef[5,.]').<=xs1[i])+1;      
	            incom[i]=r*kind0[i,1]+ef[5]*w*hz[s0];
	            earn[i]=ef[5]*w*hz[s0];
	        endif;
            if typ1[5]<i and i<=typ1[6];
                ppgg=ppgg1[2,1];
                ppgb=ppgb1[2,1];
                ppbg=ppbg1[2,1];
                ppbb=ppbb1[2,1];
                j0=sumc(cumsumc(ppef[1,.]').<=xs1[i])+1;          
	            incom[i]=r*kind0[i,1]+w0;
	            earn[i]=w0;
	        endif;
            if typ1[6]<i and i<=typ1[7];
                ppgg=ppgg2[2,1];
                ppgb=ppgb2[2,1];
                ppbg=ppbg2[2,1];
                ppbb=ppbb2[2,1];
                j0=sumc(cumsumc(ppef[2,.]').<=xs1[i])+1;    
	            incom[i]=r*kind0[i,1]+w0;      
	            earn[i]=w0;
	        endif;
            if typ1[7]<i and i<=typ1[8];
                ppgg=ppgg3[2,1];
                ppgb=ppgb3[2,1];
                ppbg=ppbg3[2,1];
                ppbb=ppbb3[2,1];
                j0=sumc(cumsumc(ppef[3,.]').<=xs1[i])+1;          
	            incom[i]=r*kind0[i,1]+w0;
	            earn[i]=w0;
	        endif;
            if typ1[8]<i and i<=typ1[9];
                ppgg=ppgg4[2,1];
                ppgb=ppgb4[2,1];
                ppbg=ppbg4[2,1];
                ppbb=ppbb4[2,1];
                j0=sumc(cumsumc(ppef[4,.]').<=xs1[i])+1;          
	            incom[i]=r*kind0[i,1]+w0;
	            earn[i]=w0;
            endif;
	        if i>=typ1[9];
                ppgg=ppgg5[2,1];
                ppgb=ppgb5[2,1];
                ppbg=ppbg5[2,1];
                ppbb=ppbb5[2,1];            
                j0=sumc(cumsumc(ppef[5,.]').<=xs1[i])+1;
	            incom[i]=r*kind0[i,1]+w0;
	            earn[i]=w0;
            endif;
  
            if kt[t-1,2]==zg and kt[t,2]==zg;
                if xs[i]<ppgg;  /* employed in period t+1 */
                    kind1[i,2]=j0;
                else;
                    kind1[i,2]=5+j0;
                endif;
            elseif kt[t-1,2]==zg and kt[t,2]==zb;
                if xs[i]<ppgb;
                    kind1[i,2]=j0;
                else;
                    kind1[i,2]=5+j0;
                endif;
             elseif kt[t-1,2]==zb and kt[t,2]==zg;
                if xs[i]<ppbg;
                    kind1[i,2]=j0;
                else;
                    kind1[i,2]=5+j0;
                endif;   
             else;
                if xs[i]<ppbb;
                    kind1[i,2]=j0;
                else;
                    kind1[i,2]=5+j0;
                endif;
             endif;   
        endo; /* i=1,..,nh */

        /* random order of households */
	    kind2=kind1~xs2;
	    kind3=sortc(kind2,3);	
	    kind4=kind3[.,1:2];
        kind1=sortc(kind4,2);
	    kind0=kind1;
	    kt[t,1]=meanc(kind1[.,1]);
	    incom1=sortc(incom,1);
	    income2=cumsumc(incom1);
	    income2=income2/income2[nh];
	    n1=round(nh/5);
	    incomet[t-1,1]=income2[n1];
	    n2=round(2*nh/5);
	    incomet[t-1,2]=income2[n2]-income2[n1];
	    n3=round(3*nh/5);
	    incomet[t-1,3]=income2[n3]-income2[n2];
	    n4=round(4*nh/5);
	    incomet[t-1,4]=income2[n4]-income2[n3];
	    n5=round(0.95*nh);
	    incomet[t-1,5]=income2[n5]-income2[n4];
	    incomet[t-1,6]=1-income2[n5];

        earn1=sortc(earn,1);
	    earn2=cumsumc(earn1);
	    earn2=earn2/earn2[nh];
	    n1=round(nh/5);
	    earnt[t-1,1]=earn2[n1];
	    n2=round(2*nh/5);
	    earnt[t-1,2]=earn2[n2]-earn2[n1];
	    n3=round(3*nh/5);
	    earnt[t-1,3]=earn2[n3]-earn2[n2];
	    n4=round(4*nh/5);
	    earnt[t-1,4]=earn2[n4]-earn2[n3];
	    n5=round(0.95*nh);
	    earnt[t-1,5]=earn2[n5]-earn2[n4];
	    earnt[t-1,6]=1-earn2[n5];
	    yt[t-1]=kt[t-1,2]*kt[t-1,1]^(alpha)*(nn0[s0]*hz[s0])^(1-alpha);

        gk=sortc(kind1,1);
        gk=gk[.,1]~seqa(1/nh,0,nh);


        /* computation of the gini coefficient of capital distribution */
        gkbar=meanc(gk[.,1]);
        j=1;
        fj=zeros(nh,1);
        fj[1]=gk[1,2];
        hj=zeros(nh,1);
        hj[1]=fj[1]*gk[1,1]/gkbar;
        gini=1-hj[1]*gk[1,1];
        do until j==nh;
            j=j+1;
            fj[j]=gk[j,2];
            hj[j]=hj[j-1]+fj[j]*gk[j,1]/gkbar;
            gini=gini-(hj[j]+hj[j-1])*fj[j];
        endo;
        giniw0=(2*(gk[.,1]'*seqa(1,1,nh))/sumc(kind1[.,1])-nh-1)/nh;
        "wealth gini: " gini~giniw0;
        giniw[t]=giniw0;

        /* computation of the gini coefficient of income distribution */
        ybar=meanc(incom1);
        j=1;
        fj=zeros(nh,1);
        fj[1]=1/nh;
        hj=zeros(nh,1);
        hj[1]=fj[1]*incom1[1]/ybar;
        gini=1-hj[1]*incom1[1];
        do until j==nh;
            j=j+1;
            fj[j]=1/nh;
            hj[j]=hj[j-1]+fj[j]*incom1[j]/ybar;
            gini=gini-(hj[j]+hj[j-1])*fj[j];
        endo;
        incom3=incom1/sumc(incom1);
        giniy0=(2*(incom3'*seqa(1,1,nh))-nh-1)/nh;
        "income gini: " gini~giniy0;
        giniy[t]=giniy0;
    
        /* computation of the gini coefficient of earnings distribution */
        ebar=meanc(earn1);
        j=1;
        fj=zeros(nh,1);
        fj[1]=1/nh;
        hj=zeros(nh,1);
        hj[1]=fj[1]*earn1[1]/ebar;
        gini=1-hj[1]*earn1[1];
        do until j==nh;
            j=j+1;
            fj[j]=1/nh;
            hj[j]=hj[j-1]+fj[j]*earn1[j]/ebar;
            gini=gini-(hj[j]+hj[j-1])*fj[j];
        endo;
        earn3=earn1/sumc(earn1);
        ginie0=(2*(earn3'*seqa(1,1,nh))-nh-1)/nh;
        "earnings gini: " gini~ginie0;
        ginie[t]=ginie0;
   endo;  /* t=1,..,nt */






    if q==nq;
        graphset;
	    title("capital over time");
        xlabel("");
        ylabel("");
	    xy(seqa(1,1,rows(kt[.,1])),kt[.,1]);
        wait;
	    title("income shares over time");
	    xy(seqa(1,1,rows(incomet)),incomet);
        wait;
	    title("production y ");
	    xy(seqa(1,1,rows(yt)),yt);
        wait;
	    title("distribution of wealth in t");
	    xy(gk[.,1],cumsumc(gk[.,2]));
        wait;
	    title("distribution of income in t");   
	    xy(incom1,seqa(1/nh,1/nh,nh));
        wait;
        title("distribution of earnings in t");
	    xy(earn1,seqa(1/nh,1/nh,nh));
        wait;
	    title("Lorenz curve earnings");
	    xy(seqa(1/nh,1/nh,nh),cumsumc(earn1)/sumc(earn1)~seqa(1/nh,1/nh,nh));
        wait;
	    title("time series of wealth gini");
	    xy(seqa(nt1,1,nt-nt1+1),giniw[nt1:nt]);
        wait;
	    title("time series of income gini");
	    xy(seqa(nt1,1,nt-nt1+1),giniy[nt1:nt]);
        wait;
	    title("time series of earnings gini");
	    xy(seqa(nt1,1,nt-nt1+1),ginie[nt1:nt]);
        wait;
	    title("Lorenz curve wealth");
        _plegctl=1;
        
        _plegctl={1 4 0.1 0.7};
        
        _plegstr="equal distribution\000Model\000US economy";
        /* table 6 in castaneda */
        x0=zeros(nh-6,2)|-0.0060~0.2|0.0180~0.2|0.057~0.2|0.1340~0.2|0.2580~0.15|0.539~0.05;
        x02=cumsumc(x0[.,1])/sumc(x0[.,1]);
        x01=cumsumc(x0[.,2])/sumc(x0[.,2]);
        z01=cumsumc(gk[.,2])/sumc(gk[.,2]);
        z02=cumsumc(gk[.,1])/sumc(gk[.,1]);
	    xy(z01~z01~x01,z01~z02~x02);
        save wz01=z01,wx01=x01,wz02=z02,wx02=x02;
        wait;
	    title("Lorenz curve income");
        x0=zeros(nh-6,2)|0.0505~0.2|0.1195~0.2|0.1756~0.2|0.2391~0.2|0.2556~0.15|0.1597~0.05;
        x02=cumsumc(x0[.,1])/sumc(x0[.,1]);
        x01=cumsumc(x0[.,2])/sumc(x0[.,2]);
        z01=seqa(1/nh,1/nh,nh);
        z01=z01/z01[nh];
        z02=cumsumc(incom1)/sumc(incom1);
	    xy(z01~z01~x01,z01~z02~x02);
        save iz01=z01,ix01=x01,ix02=x02,iz02=z02;
        wait;
        graphset;

    /* correlation of income shares and output, table 4 in castaneda,
    yearly observations */ 

        "Correlation of income shares with logged,detrended output: ";

        yt1=yt[nt1:nt-1]; yt1=reshape(yt1,(nt-nt1-1)/8,8); yt1=yt1';
        yt1=meanc(yt1);

        incomet1=incomet[nt1:nt-1,1];
        incomet1=reshape(incomet1,(nt-nt1-1)/8,8); incomet1=incomet1';
        incomet1=meanc(incomet1);

        incomet2=incomet[nt1:nt-1,2];
        incomet2=reshape(incomet2,(nt-nt1-1)/8,8); incomet2=incomet2';
        incomet2=meanc(incomet2);

        incomet3=incomet[nt1:nt-1,3];
        incomet3=reshape(incomet3,(nt-nt1-1)/8,8); incomet3=incomet3';
        incomet3=meanc(incomet3);

        incomet4=incomet[nt1:nt-1,4];
        incomet4=reshape(incomet4,(nt-nt1-1)/8,8); incomet4=incomet4';
        incomet4=meanc(incomet4);

        incomet5=incomet[nt1:nt-1,5];
        incomet5=reshape(incomet5,(nt-nt1-1)/8,8); incomet5=incomet5';
        incomet5=meanc(incomet5);

        incomet6=incomet[nt1:nt-1,6];
        incomet6=reshape(incomet6,(nt-nt1-1)/8,8); incomet6=incomet6';
        incomet6=meanc(incomet6);


        et1=earnt[nt1:nt-1,1]; 
        et1=reshape(et1,(nt-nt1-1)/8,8); et1=et1';
        et1=meanc(et1);

        et2=earnt[nt1:nt-1,2]; et2=reshape(et2,(nt-nt1-1)/8,8); et2=et2';
        et2=meanc(et2);

        et3=earnt[nt1:nt-1,3]; et3=reshape(et3,(nt-nt1-1)/8,8); et3=et3';
        et3=meanc(et3);

        et4=earnt[nt1:nt-1,4]; et4=reshape(et4,(nt-nt1-1)/8,8); et4=et4';
        et4=meanc(et4);

        et5=earnt[nt1:nt-1,5]; et5=reshape(et5,(nt-nt1-1)/8,8); et5=et5';
        et5=meanc(et5);

        et6=earnt[nt1:nt-1,6]; et6=reshape(et6,(nt-nt1-1)/8,8); et6=et6';
        et6=meanc(et6);


        {x0,yt1}=hptrend(ln(yt1),100);
        corin=corrx(incomet1~incomet2~incomet3~incomet4~incomet5~incomet6~yt1);
        corin=corin[7,1:6];

        corearn=corrx(et1~et2~et3~et4~et5~et6~yt1);
        corearn=corearn[7,1:6];


        save corin,corearn; 
        "correlation of income and logged, detrended yearly output"; 
        corin; 
        "correlation of earnings and logged,detrended yearly output:"; 
        corearn;

        "standard deviation of yearly output: "
        stdc(yt1); wait;
        "standard deviation of yearly income shares 1-6: ";
        stdc(incomet1);
        stdc(incomet2);
        stdc(incomet3);
        stdc(incomet4);
        stdc(incomet5);
        stdc(incomet6);
        "standard deviation of yearly earning shares 1-6:";
        stdc(et1);
        stdc(et2);
        stdc(et3);
        stdc(et4);
        stdc(et5);
        stdc(et6);
        "Press any key...";
        wait;

    endif;




    "Gini coefficient of wealth: ";
    meanc(giniw[nt1:nt]);

    "Gini coefficient of income: ";
    meanc(giniy[nt1:nt]);

    "Gini coefficient of earnings: ";
    meanc(ginie[nt1:nt]);

    /* transition dynamics for z=zg and z=zb */
    ktx=kt[nt1:nt-1,.]~kt[nt1+1:nt,1];
    kmean0=meanc(kt[nt1:nt,1]);
    kq[q]=kmean0;
    /* update of mean capital stock */
    kmean=phi1*kmean+(1-phi1)*kmean0;
    ktx1=selif(ktx,ktx[.,2].==zg);
    ktx2=selif(ktx,ktx[.,2].==zb);

    /* ols-estimate gam0,gam1 */
    xi=ln(ktx1[.,1]); 
    xi=ones(rows(xi),1)~xi;
    yi=ln(ktx1[.,3]);
    betai=inv(xi'*xi)*xi'*yi;
    gam01g=betai[1];
    gam11g=betai[2];
    xi=ln(ktx2[.,1]); 
    xi=ones(rows(xi),1)~xi;
    yi=ln(ktx2[.,3]);
    betai=inv(xi'*xi)*xi'*yi;
    gam01b=betai[1];
    gam11b=betai[2];
    "gam0g~gam1g: " gam01g~gam11g;
    "gam0b~gam1b: " gam01b~gam11b;
    
    if q==nq;
        title("fit of capital stock law of motion for z==zg");
        xlabel("k");
        ylabel("k'");
        xy(ktx1[.,1],ktx1[.,3]~exp(gam01g+gam11g*ln(ktx1[.,1])));
        wait;
        title("fit of capital stock law of motion for z==zb");
        xlabel("k");
        ylabel("k'");
        xy(ktx2[.,1],ktx2[.,3]~exp(gam01b+gam11b*ln(ktx2[.,1])));
        wait;
    endif;

    save kkt=kt;
    kritgam=abs(ln(gam0g)-ln(gam01g)|ln(gam1g)-ln(gam11g)|ln(gam0b)-ln(gam01b)|ln(gam1b)-ln(gam11b));
    "kritgam: " kritgam';
    gam0g=phi1*gam0g+(1-phi1)*gam01g;
    gam1g=phi1*gam1g+(1-phi1)*gam11g;
    gam0b=phi1*gam0b+(1-phi1)*gam01b;
    gam1b=phi1*gam1b+(1-phi1)*gam11b;
    save kgam0g=gam0g,kgam1g=gam1g; 
    save kgam0b=gam0b,kgam1b=gam1b;
    gam0q[q,.]=gam0g~gam0b;
    gam1q[q,.]=gam1g~gam1b;
    kq[q]=kmean0;
endo;   /* q=1,..,nq */

"convergence after iteration: " q;

"press any key to continue..";
wait;	


        graphset;
	    title("capital over time");
        xlabel("");
        ylabel("");
	    xy(seqa(1,1,rows(kt[.,1])),kt[.,1]);
        wait;
	    title("income shares over time");
	    xy(seqa(1,1,rows(incomet)),incomet);
        wait;
	    title("production y ");
	    xy(seqa(1,1,rows(yt)),yt);
        wait;
	    title("distribution of wealth in t");
	    xy(gk[.,1],cumsumc(gk[.,2]));
        wait;
	    title("distribution of income in t");   
	    xy(incom1,seqa(1/nh,1/nh,nh));
        wait;
        title("distribution of earnings in t");
	    xy(earn1,seqa(1/nh,1/nh,nh));
        wait;
	    title("Lorenz curve earnings");
	    xy(seqa(1/nh,1/nh,nh),cumsumc(earn1)/sumc(earn1)~seqa(1/nh,1/nh,nh));
        wait;
	    title("time series of wealth gini");
	    xy(seqa(nt1,1,nt-nt1+1),giniw[nt1:nt]);
        wait;
	    title("time series of income gini");
	    xy(seqa(nt1,1,nt-nt1+1),giniy[nt1:nt]);
        wait;
	    title("time series of earnings gini");
	    xy(seqa(nt1,1,nt-nt1+1),ginie[nt1:nt]);
        wait;
	    title("Lorenz curve wealth");
        _plegctl=1;
        
        _plegctl={1 4 0.1 0.7};
        
        _plegstr="equal distribution\000Model\000US economy";
        /* table 6 in castaneda */
        x0=zeros(nh-6,2)|-0.0060~0.2|0.0180~0.2|0.057~0.2|0.1340~0.2|0.2580~0.15|0.539~0.05;
        x02=cumsumc(x0[.,1])/sumc(x0[.,1]);
        x01=cumsumc(x0[.,2])/sumc(x0[.,2]);
        z01=cumsumc(gk[.,2])/sumc(gk[.,2]);
        z02=cumsumc(gk[.,1])/sumc(gk[.,1]);
	    xy(z01~z01~x01,z01~z02~x02);
        save wz01=z01,wx01=x01,wz02=z02,wx02=x02;
        wait;
	    title("Lorenz curve income");
        x0=zeros(nh-6,2)|0.0505~0.2|0.1195~0.2|0.1756~0.2|0.2391~0.2|0.2556~0.15|0.1597~0.05;
        x02=cumsumc(x0[.,1])/sumc(x0[.,1]);
        x01=cumsumc(x0[.,2])/sumc(x0[.,2]);
        z01=seqa(1/nh,1/nh,nh);
        z01=z01/z01[nh];
        z02=cumsumc(incom1)/sumc(incom1);
	    xy(z01~z01~x01,z01~z02~x02);
        save iz01=z01,ix01=x01,iz02=z02,ix02=x02;
        wait;
        graphset;

    /* correlation of income shares and output, table 4 in castaneda,
    yearly observations */ 

        "Correlation of income shares with logged,detrended output: ";

        yt1=yt[nt1:nt-1]; yt1=reshape(yt1,(nt-nt1-1)/8,8); yt1=yt1';
        yt1=meanc(yt1);

        incomet1=incomet[nt1:nt-1,1];
        incomet1=reshape(incomet1,(nt-nt1-1)/8,8); incomet1=incomet1';
        incomet1=meanc(incomet1);

        incomet2=incomet[nt1:nt-1,2];
        incomet2=reshape(incomet2,(nt-nt1-1)/8,8); incomet2=incomet2';
        incomet2=meanc(incomet2);

        incomet3=incomet[nt1:nt-1,3];
        incomet3=reshape(incomet3,(nt-nt1-1)/8,8); incomet3=incomet3';
        incomet3=meanc(incomet3);

        incomet4=incomet[nt1:nt-1,4];
        incomet4=reshape(incomet4,(nt-nt1-1)/8,8); incomet4=incomet4';
        incomet4=meanc(incomet4);

        incomet5=incomet[nt1:nt-1,5];
        incomet5=reshape(incomet5,(nt-nt1-1)/8,8); incomet5=incomet5';
        incomet5=meanc(incomet5);

        incomet6=incomet[nt1:nt-1,6];
        incomet6=reshape(incomet6,(nt-nt1-1)/8,8); incomet6=incomet6';
        incomet6=meanc(incomet6);


        et1=earnt[nt1:nt-1,1]; 
        et1=reshape(et1,(nt-nt1-1)/8,8); et1=et1';
        et1=meanc(et1);

        et2=earnt[nt1:nt-1,2]; et2=reshape(et2,(nt-nt1-1)/8,8); et2=et2';
        et2=meanc(et2);

        et3=earnt[nt1:nt-1,3]; et3=reshape(et3,(nt-nt1-1)/8,8); et3=et3';
        et3=meanc(et3);

        et4=earnt[nt1:nt-1,4]; et4=reshape(et4,(nt-nt1-1)/8,8); et4=et4';
        et4=meanc(et4);

        et5=earnt[nt1:nt-1,5]; et5=reshape(et5,(nt-nt1-1)/8,8); et5=et5';
        et5=meanc(et5);

        et6=earnt[nt1:nt-1,6]; et6=reshape(et6,(nt-nt1-1)/8,8); et6=et6';
        et6=meanc(et6);


        {x0,yt1}=hptrend(ln(yt1),100);
        corin=corrx(incomet1~incomet2~incomet3~incomet4~incomet5~incomet6~yt1);
        corin=corin[7,1:6];

        corearn=corrx(et1~et2~et3~et4~et5~et6~yt1);
        corearn=corearn[7,1:6];


        save corin,corearn; 
        "correlation of income and logged, detrended yearly output"; 
        corin; 
        "correlation of earnings and logged,detrended yearly output:"; 
        corearn;

        "standard deviation of yearly output: "
        stdc(yt1); wait;
        "standard deviation of yearly income shares 1-6: ";
        stdc(incomet1);
        stdc(incomet2);
        stdc(incomet3);
        stdc(incomet4);
        stdc(incomet5);
        stdc(incomet6);
        "standard deviation of yearly earning shares 1-6:";
        stdc(et1);
        stdc(et2);
        stdc(et3);
        stdc(et4);
        stdc(et5);
        stdc(et6);
        "Press any key...";
        wait;


"Gini coefficient of wealth: ";
meanc(giniw[nt1:nt]);

"Gini coefficient of income: ";
meanc(giniy[nt1:nt]);

"Gini coefficient of earnings: ";
meanc(ginie[nt1:nt]);

"gam0g~gam1g: " gam01g~gam11g;
"gam0b~gam1b: " gam01b~gam11b;
    
"Deviation of value function: ";
crit;

"press any key to continue..";
wait;


graphset;
title("aopt1eg-a");
xlabel("a");
ylabel("");
xy(a,aopt1eg-a);
wait;
title("aop2teg-a");
xy(a,aopt2eg-a);
wait;
title("aopt3eg-a");
xy(a,aopt3eg-a);
wait;
title("aopt4eg-a");
xy(a,aopt4eg-a);
wait;
title("aopt5eg-a");
xy(a,aopt5eg-a);
wait;
title("aopt1ug-a");
xy(a,aopt1ug-a);
wait;
title("aopt2ug-a");
xy(a,aopt2ug-a);
wait;
title("aopt3ug-a");
xy(a,aopt3ug-a);
wait;
title("aopt4ug-a");
xy(a,aopt4ug-a);
wait;
title("aopt5ug-a");
xy(a,aopt5ug-a);
wait;
title("aopt1eb-a");
xy(a,aopt1eb-a);
wait;
title("aopt2eb-a");
xy(a,aopt2eb-a);
wait;
title("aopt3eb-a");
xy(a,aopt3eb-a);
wait;
title("aopt4eb-a");
xy(a,aopt4eb-a);
wait;
title("aopt5eb-a");
xy(a,aopt5eb-a);
wait;
title("aopt1ub-a");
xy(a,aopt1ub-a);
wait;
title("aopt2ub-a");
xy(a,aopt2ub-a);
wait;
title("aopt3ub-a");
xy(a,aopt3ub-a);
wait;
title("aopt4ub-a");
xy(a,aopt4ub-a);
wait;
title("aopt5ub-a");
xy(a,aopt5ub-a);
wait;


@  ----------------------------  procedures -----------


u(x) -- utility function

value(a,e) -- returns the value of the value function for asset
                a and employment status e

value1(x) -- given a=a[i] and epsilon=e, returns the
             value of the bellman equation for a'=x

bellman -- value for the right-hand side of the Bellman equation

transp - computes the transition matrix for given equilibrium 
		employment

------------------------------------------------------- @


proc transp(x);
	local a11,a12,a21,a22,z,a,transp;
	a11=x[1];
	a12=x[2];
	a21=x[3];
	a22=x[4];
	transp=x;
	transp[1]=a11+a12-1;
	transp[2]=a21+a22-1;
	a=(a11~a12|a21~a22);
	z=equivec1(a);
	transp[3]=a22-udur;
	transp[4]=z[1]-nrate;
	retp(transp);
endp;

proc u(x);
   retp(x^(1-sigma)/(1-sigma));
endp;


proc value1(x);
    retp(bellman(a[i],x,h));
endp;


proc bellman(a0,a1,y);
    local c,bell,bell1,bell2,ppgg,ppgb,ppbg,ppbb;
   
    if a1<a[1]; 
	    retp((1+a1^2)*neg); 
    endif; 

    if e==1;
        c=(1+r)*a0+ef[y]*hz[s]*w-a1;
    else;
        c=(1+r)*a0+w0-a1;
    endif;

    if c<0;
        retp((c^2+1)*neg);
    endif;
  
    bell=u(c);
    bell1=0;
    bell2=0;
    if y==1;
        ppgg=ppgg1;
	    ppgb=ppgb1;
	    ppbg=ppbg1;
	    ppbb=ppbb1;
    elseif y==2;
	    ppgg=ppgg2;
	    ppgb=ppgb2;
	    ppbg=ppbg2;
	    ppbb=ppbb2;
    elseif y==3;
	    ppgg=ppgg3;
	    ppgb=ppgb3;
	    ppbg=ppbg3;
	    ppbb=ppbb3;
    elseif y==4;
	    ppgg=ppgg4;
	    ppgb=ppgb4;
	    ppbg=ppbg4;
	    ppbb=ppbb4;
    elseif y==5;
	    ppgg=ppgg5;
	    ppgb=ppgb5;
	    ppbg=ppbg5;
	    ppbb=ppbb5;
    endif;

    if a1>=a[na];    
        if zz==zg;
            bell1=bell1+beta*ppgg[e,1]*ppef[y,1]*v1eg1[na]+beta*ppgg[e,2]*ppef[y,1]*v1ug1[na];
            bell1=bell1+beta*ppgg[e,1]*ppef[y,2]*v2eg1[na]+beta*ppgg[e,2]*ppef[y,2]*v2ug1[na];
            bell1=bell1+beta*ppgg[e,1]*ppef[y,3]*v3eg1[na]+beta*ppgg[e,2]*ppef[y,3]*v3ug1[na];
            bell1=bell1+beta*ppgg[e,1]*ppef[y,4]*v4eg1[na]+beta*ppgg[e,2]*ppef[y,4]*v4ug1[na];
            bell1=bell1+beta*ppgg[e,1]*ppef[y,5]*v5eg1[na]+beta*ppgg[e,2]*ppef[y,5]*v5ug1[na];
 	        bell=bell+ppz[1,1]*bell1;
	
            bell2=bell2+beta*ppgb[e,1]*ppef[y,1]*v1eb1[na]+beta*ppgb[e,2]*ppef[y,1]*v1ub1[na];
            bell2=bell2+beta*ppgb[e,1]*ppef[y,2]*v2eb1[na]+beta*ppgb[e,2]*ppef[y,2]*v2ub1[na];
            bell2=bell2+beta*ppgb[e,1]*ppef[y,3]*v3eb1[na]+beta*ppgb[e,2]*ppef[y,3]*v3ub1[na];
            bell2=bell2+beta*ppgb[e,1]*ppef[y,4]*v4eb1[na]+beta*ppgb[e,2]*ppef[y,4]*v4ub1[na];
            bell2=bell2+beta*ppgb[e,1]*ppef[y,5]*v5eb1[na]+beta*ppgb[e,2]*ppef[y,5]*v5ub1[na];
	        bell=bell+ppz[1,2]*bell2;
         
        else;

            bell1=bell1+beta*ppbg[e,1]*ppef[y,1]*v1eg1[na]+beta*ppbg[e,2]*ppef[y,1]*v1ug1[na];
            bell1=bell1+beta*ppbg[e,1]*ppef[y,2]*v2eg1[na]+beta*ppbg[e,2]*ppef[y,2]*v2ug1[na];
            bell1=bell1+beta*ppbg[e,1]*ppef[y,3]*v3eg1[na]+beta*ppbg[e,2]*ppef[y,3]*v3ug1[na];
            bell1=bell1+beta*ppbg[e,1]*ppef[y,4]*v4eg1[na]+beta*ppbg[e,2]*ppef[y,4]*v4ug1[na];
            bell1=bell1+beta*ppbg[e,1]*ppef[y,5]*v5eg1[na]+beta*ppbg[e,2]*ppef[y,5]*v5ug1[na];
            bell=bell+ppz[2,1]*bell1;

            bell2=bell2+beta*ppbb[e,1]*ppef[y,1]*v1eb1[na]+beta*ppbb[e,2]*ppef[y,1]*v1ub1[na];
            bell2=bell2+beta*ppbb[e,1]*ppef[y,2]*v2eb1[na]+beta*ppbb[e,2]*ppef[y,2]*v2ub1[na];
            bell2=bell2+beta*ppbb[e,1]*ppef[y,3]*v3eb1[na]+beta*ppbb[e,2]*ppef[y,3]*v3ub1[na];
            bell2=bell2+beta*ppbb[e,1]*ppef[y,4]*v4eb1[na]+beta*ppbb[e,2]*ppef[y,4]*v4ub1[na];
            bell2=bell2+beta*ppbb[e,1]*ppef[y,5]*v5eb1[na]+beta*ppbb[e,2]*ppef[y,5]*v5ub1[na];
	        bell=bell+ppz[2,2]*bell2;
        endif;
    
    elseif a1==a[1];
 
        if zz==zg;

            bell1=bell1+beta*ppgg[e,1]*ppef[y,1]*v1eg1[1]+beta*ppgg[e,2]*ppef[y,1]*v1ug1[1];
            bell1=bell1+beta*ppgg[e,1]*ppef[y,2]*v2eg1[1]+beta*ppgg[e,2]*ppef[y,2]*v2ug1[1];
            bell1=bell1+beta*ppgg[e,1]*ppef[y,3]*v3eg1[1]+beta*ppgg[e,2]*ppef[y,3]*v3ug1[1];
            bell1=bell1+beta*ppgg[e,1]*ppef[y,4]*v4eg1[1]+beta*ppgg[e,2]*ppef[y,4]*v4ug1[1];
            bell1=bell1+beta*ppgg[e,1]*ppef[y,5]*v5eg1[1]+beta*ppgg[e,2]*ppef[y,5]*v5ug1[1];
 	        bell=bell+ppz[1,1]*bell1;

            bell2=bell2+beta*ppgb[e,1]*ppef[y,1]*v1eb1[1]+beta*ppgb[e,2]*ppef[y,1]*v1ub1[1];
            bell2=bell2+beta*ppgb[e,1]*ppef[y,2]*v2eb1[1]+beta*ppgb[e,2]*ppef[y,2]*v2ub1[1];
            bell2=bell2+beta*ppgb[e,1]*ppef[y,3]*v3eb1[1]+beta*ppgb[e,2]*ppef[y,3]*v3ub1[1];
            bell2=bell2+beta*ppgb[e,1]*ppef[y,4]*v4eb1[1]+beta*ppgb[e,2]*ppef[y,4]*v4ub1[1];
            bell2=bell2+beta*ppgb[e,1]*ppef[y,5]*v5eb1[1]+beta*ppgb[e,2]*ppef[y,5]*v5ub1[1];
            bell=bell+ppz[1,2]*bell2;
 
        else;

            bell1=bell1+beta*ppbg[e,1]*ppef[y,1]*v1eg1[1]+beta*ppbg[e,2]*ppef[y,1]*v1ug1[1];
            bell1=bell1+beta*ppbg[e,1]*ppef[y,2]*v2eg1[1]+beta*ppbg[e,2]*ppef[y,2]*v2ug1[1];
            bell1=bell1+beta*ppbg[e,1]*ppef[y,3]*v3eg1[1]+beta*ppbg[e,2]*ppef[y,3]*v3ug1[1];
            bell1=bell1+beta*ppbg[e,1]*ppef[y,4]*v4eg1[1]+beta*ppbg[e,2]*ppef[y,4]*v4ug1[1];
            bell1=bell1+beta*ppbg[e,1]*ppef[y,5]*v5eg1[1]+beta*ppbg[e,2]*ppef[y,5]*v5ug1[1];
            bell=bell+ppz[2,1]*bell1;
        
            bell2=bell2+beta*ppbb[e,1]*ppef[y,1]*v1eb1[1]+beta*ppbb[e,2]*ppef[y,1]*v1ub1[1];
            bell2=bell2+beta*ppbb[e,1]*ppef[y,2]*v2eb1[1]+beta*ppbb[e,2]*ppef[y,2]*v2ub1[1];
            bell2=bell2+beta*ppbb[e,1]*ppef[y,3]*v3eb1[1]+beta*ppbb[e,2]*ppef[y,3]*v3ub1[1];
            bell2=bell2+beta*ppbb[e,1]*ppef[y,4]*v4eb1[1]+beta*ppbb[e,2]*ppef[y,4]*v4ub1[1];
            bell2=bell2+beta*ppbb[e,1]*ppef[y,5]*v5eb1[1]+beta*ppbb[e,2]*ppef[y,5]*v5ub1[1];
            bell=bell+ppz[2,2]*bell2;
        endif;

    else;

        if zz==zg;

            bell1=bell1+beta*ppgg[e,1]*ppef[y,1]*lininter(a,v1eg1,a1)+beta*ppgg[e,2]*ppef[e,1]*lininter(a,v1ug1,a1);
            bell1=bell1+beta*ppgg[e,1]*ppef[y,2]*lininter(a,v2eg1,a1)+beta*ppgg[e,2]*ppef[e,2]*lininter(a,v2ug1,a1);
            bell1=bell1+beta*ppgg[e,1]*ppef[y,3]*lininter(a,v3eg1,a1)+beta*ppgg[e,2]*ppef[e,3]*lininter(a,v3ug1,a1);
            bell1=bell1+beta*ppgg[e,1]*ppef[y,4]*lininter(a,v4eg1,a1)+beta*ppgg[e,2]*ppef[e,4]*lininter(a,v4ug1,a1);
            bell1=bell1+beta*ppgg[e,1]*ppef[y,5]*lininter(a,v5eg1,a1)+beta*ppgg[e,2]*ppef[e,5]*lininter(a,v5ug1,a1);
    	    bell=bell+ppz[1,1]*bell1;

            bell2=bell2+beta*ppgb[e,1]*ppef[y,1]*lininter(a,v1eb1,a1)+beta*ppgb[e,2]*ppef[y,1]*lininter(a,v1ub1,a1);
            bell2=bell2+beta*ppgb[e,1]*ppef[y,2]*lininter(a,v2eb1,a1)+beta*ppgb[e,2]*ppef[y,2]*lininter(a,v2ub1,a1);
            bell2=bell2+beta*ppgb[e,1]*ppef[y,3]*lininter(a,v3eb1,a1)+beta*ppgb[e,2]*ppef[y,3]*lininter(a,v3ub1,a1);
            bell2=bell2+beta*ppgb[e,1]*ppef[y,4]*lininter(a,v4eb1,a1)+beta*ppgb[e,2]*ppef[y,4]*lininter(a,v4ub1,a1);
            bell2=bell2+beta*ppgb[e,1]*ppef[y,5]*lininter(a,v5eb1,a1)+beta*ppgb[e,2]*ppef[y,5]*lininter(a,v5ub1,a1);
	        bell=bell+ppz[1,2]*bell2;

        else;
 
            bell1=bell1+beta*ppbg[e,1]*ppef[y,1]*lininter(a,v1eg1,a1)+beta*ppbg[e,2]*ppef[y,1]*lininter(a,v1ug1,a1);
            bell1=bell1+beta*ppbg[e,1]*ppef[y,2]*lininter(a,v2eg1,a1)+beta*ppbg[e,2]*ppef[y,2]*lininter(a,v2ug1,a1);
            bell1=bell1+beta*ppbg[e,1]*ppef[y,3]*lininter(a,v3eg1,a1)+beta*ppbg[e,2]*ppef[y,3]*lininter(a,v3ug1,a1);
            bell1=bell1+beta*ppbg[e,1]*ppef[y,4]*lininter(a,v4eg1,a1)+beta*ppbg[e,2]*ppef[y,4]*lininter(a,v4ug1,a1);
            bell1=bell1+beta*ppbg[e,1]*ppef[y,5]*lininter(a,v5eg1,a1)+beta*ppbg[e,2]*ppef[y,5]*lininter(a,v5ug1,a1);
	        bell=bell+ppz[2,1]*bell1;

            bell2=bell2+beta*ppbb[e,1]*ppef[y,1]*lininter(a,v1eb1,a1)+beta*ppbb[e,2]*ppef[y,1]*lininter(a,v1ub1,a1);
            bell2=bell2+beta*ppbb[e,1]*ppef[y,2]*lininter(a,v2eb1,a1)+beta*ppbb[e,2]*ppef[y,2]*lininter(a,v2ub1,a1);
            bell2=bell2+beta*ppbb[e,1]*ppef[y,3]*lininter(a,v3eb1,a1)+beta*ppbb[e,2]*ppef[y,3]*lininter(a,v3ub1,a1);
            bell2=bell2+beta*ppbb[e,1]*ppef[y,4]*lininter(a,v4eb1,a1)+beta*ppbb[e,2]*ppef[y,4]*lininter(a,v4ub1,a1);
            bell2=bell2+beta*ppbb[e,1]*ppef[y,5]*lininter(a,v5eb1,a1)+beta*ppbb[e,2]*ppef[y,5]*lininter(a,v5ub1,a1);
	        bell=bell+ppz[2,2]*bell2;

        endif;  
    endif;

    retp(bell);
endp;

proc lininter(xd,yd,x);
  local j;
  j=sumc(xd.<=x');
  retp(yd[j]+(yd[j+1]-yd[j]).*(x-xd[j])./(xd[j+1]-xd[j]));
endp;

proc rrate(x,y,z);
    retp(alpha*x*y^(alpha-1)*z^(1-alpha)-delta);
endp;

proc wage(x,y,z);
    retp((1-alpha)*x*y^(alpha)*(z)^(-alpha));
endp;


proc equivec1(p);
    local n,x;
    n=rows(p);
    p=diagrv(p,diag(p)-ones(n,1));
    p=p[.,1:n-1]~ones(n,1);
    x=zeros(n-1,1)|1;
    retp((x'*inv(p))');
endp;

proc bilina(a0,k0,i,z0); /* bilinear interpolation */
    local n0,n1,n2,m0,m1,m2,phi,f,d,z1,z2,y,pol;
       if z0==zg;
       if i<=typ1[1];
                pol=aopt1eg;
       elseif typ1[1]<i and i<=typ1[2];
            pol=aopt2eg;
       elseif typ1[2]<i and i<=typ1[3];
            pol=aopt3eg;
       elseif typ1[3]<i and i<=typ1[4];
            pol=aopt4eg;
       elseif typ1[4]<i and i<=typ1[5];
            pol=aopt5eg;
       elseif typ1[5]<i and i<=typ1[6];
            pol=aopt1ug;
       elseif typ1[6]<i and i<=typ1[7];
            pol=aopt2ug;
       elseif typ1[7]<i and i<=typ1[8];
            pol=aopt3ug;
       elseif typ1[8]<i and i<=typ1[9];
            pol=aopt4ug;
       elseif typ1[9]<i and i<=typ1[10];
            pol=aopt5ug;
       endif;
        else;
       if i<=typ1[1];
                pol=aopt1eb;
       elseif typ1[1]<i and i<=typ1[2];
            pol=aopt2eb;
       elseif typ1[2]<i and i<=typ1[3];
            pol=aopt3eb;
       elseif typ1[3]<i and i<=typ1[4];
            pol=aopt4eb;
       elseif typ1[4]<i and i<=typ1[5];
            pol=aopt5eb;
       elseif typ1[5]<i and i<=typ1[6];
            pol=aopt1ub;
       elseif typ1[6]<i and i<=typ1[7];
            pol=aopt2ub;
       elseif typ1[7]<i and i<=typ1[8];
            pol=aopt3ub;
       elseif typ1[8]<i and i<=typ1[9];
            pol=aopt4ub;
       elseif typ1[9]<i and i<=typ1[10];
            pol=aopt5ub;
       endif;
        endif;

    if a0<amin1; a0=amin1; endif;
    if a0>amax1; a0=amax1; endif;
    if k0<kmin; k0=kmin;  endif;
    if k0>kmax; k0=kmax; endif;

    m0=(a0-amin1)/(amax1-amin1)*(na-1)+1;
    m2=floor(m0);
    m1=m0-m2;
    n0=(k0-kmin)/(kmax-kmin)*(nk-1)+1;
    n2=floor(n0);
    n1=n0-n2;
    phi=zeros(4,1);
    f=zeros(4,1);
    d=zeros(4,1);
    
    /* bilinear interpolation */
        if m2<na and n2<nk;
            d=k[n2+1]|k[n2]|a[m2+1]|a[m2];
            z1=(k0-d[2])/(d[1]-d[2]);
            z2=(a0-d[4])/(d[3]-d[4]);
            phi[1]=(1-z1)*(1-z2);
            phi[2]=z1*(1-z2);
            phi[3]=z1*z2;
            phi[4]=(1-z1)*z2;
            f[1]=pol[m2,n2];
            f[2]=pol[m2,n2+1];
            f[3]=pol[m2+1,n2+1];
            f[4]=pol[m2+1,n2];
            retp(phi'*f);
        elseif n2==nk and m2==na;
            retp(pol[m2,n2]);
        elseif n2==nk and m2<na;  /*linear interpolation */
            retp((1-m1)*pol[m2,nk]+m1*pol[m2+1,nk]);
        elseif n2<nk and m2==na;  /*linear interpolation */
            retp((1-n1)*pol[na,n2]+n1*pol[na,n2+1]);
     endif;
endp;

proc matroot(x,n);
	local s,u,i,j,y;
	{s,u}=eigv(x);
	s=diagrv(eye(rows(x)),s^(1/n));
	y=u*s*inv(u);
	i=0;
	do until i==rows(x);
		i=i+1;
		j=0;
		do until j==rows(x);
			j=j+1;
			if y[i,j]<0; y[i,j]=0; endif;
		endo;
			y[i,.]=y[i,.]/sumc(y[i,.]');
	endo;
	retp(y);
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


proc rf2(x,y);
    local rf2,a0,k0,c0,cex0,r,zz,w,k1,a1,s1,cex,e1,h1;
    local c1,pp1,x01,zz1,a2,r1,w1;
    rf2=0;
    a0=x;
    k0=y;
    if s==1; zz=zg; else; zz=zb; endif;
    a1=bilina(a0,k0,(e-1)*5+h,zz);
    r=rrate(zz,k0,nn0[s]*hz[s]); 
    w=wage(zz,k0,nn0[s]*hz[s]); 
    k1=exp(gam0+gam1*ln(k0));
    if e==1;
        c0=(1+r)*a0+w*hz[s]*ef[h]-a1;
    else;
        c0=(1+r)*a0+w0-a1;          
    endif;
    cex0=uc(c0);
    /* computation of E_t+1 */
    s1=0;
    cex=0;
    do until s1==2;
        s1=s1+1;
        if s1==1;
           zz1=zg;
        else;
           zz1=zb;
        endif;
        e1=0;
        do until e1==2;
           e1=e1+1;
           h1=0;
           do until h1==5;
              h1=h1+1;
              r1=rrate(zz1,k1,nn0[s1]*hz[s1]);
              w1=wage(zz1,k1,nn0[s1]*hz[s1]); 
              a2=bilina(a1,k1,(e1-1)*5+h1,zz1);
                if e1==1;
                    c1=(1+r1)*a1+w1*hz[s1]*ef[h1]-a2;
                else;
                    c1=(1+r1)*a1+w0-a2;          
                endif;
              pp1=prob(s,e,h,s1,e1,h1);
              cex=cex+beta*pp1*uc(c1)*(1+r1);
  if ismiss(cex)>0; "ismiss>0"; wait; endif; 
           endo;
        endo;         
    endo;                                    
    rf2=1-cex/cex0;
    retp(abs(rf2));    
endp;


proc acold(s1,e1,h1);
    local acr;
    if s1==1;
        if e1==1;
            if h1==1;
              acr=aopt1eg;
            elseif h1==2;
              acr=aopt2eg;
            elseif h1==3;
              acr=aopt3eg;
            elseif h1==4;
              acr=aopt4eg;
            elseif h1==5;
              acr=aopt5eg;
            endif;
        else;
            if h1==1;
              acr=aopt1ug;
            elseif h1==2;
              acr=aopt2ug;
            elseif h1==3;
              acr=aopt3ug;
            elseif h1==4;
              acr=aopt4ug;
            elseif h1==5;
              acr=aopt5ug;
            endif;
        endif;
    else;
        if e1==1;
            if h1==1;
              acr=aopt1eb;
            elseif h1==2;
              acr=aopt2eb;
            elseif h1==3;
              acr=aopt3eb;
            elseif h1==4;
              acr=aopt4eb;
            elseif h1==5;
              acr=aopt5eb;
            endif;
        else;
            if h1==1;
              acr=aopt1ub;
            elseif h1==2;
              acr=aopt2ub;
            elseif h1==3;
              acr=aopt3ub;
            elseif h1==4;
              acr=aopt4ub;
            elseif h1==5;
              acr=aopt5ub;
            endif;
        endif;
    endif;
    retp(acr);
endp;

/* probability for z_t+1=z[s1], e_t+1=e1, type=h1 for given s,e,h */
proc prob(s,e,h,s1,e1,h1);
    local ppr;
    if s==1;
        if s1==1;
            if h==1;
                ppr=ppgg1[e,e1];
            elseif h==2;
                ppr=ppgg2[e,e1];
            elseif h==3;
                ppr=ppgg3[e,e1];
            elseif h==4;
                ppr=ppgg4[e,e1];
            elseif h==5;
                ppr=ppgg5[e,e1];
            endif;
        else;
            if h==1;
                ppr=ppgb1[e,e1];
            elseif h==2;
                ppr=ppgb2[e,e1];
            elseif h==3;
                ppr=ppgb3[e,e1];
            elseif h==4;
                ppr=ppgb4[e,e1];
            elseif h==5;
                ppr=ppgb5[e,e1];
            endif;
        endif;
    else;
        if s1==1;
            if h==1;
                ppr=ppbg1[e,e1];
            elseif h==2;
                ppr=ppbg2[e,e1];
            elseif h==3;
                ppr=ppbg3[e,e1];
            elseif h==4;
                ppr=ppbg4[e,e1];
            elseif h==5;
                ppr=ppbg5[e,e1];
            endif;
        else;
            if h==1;
                ppr=ppbb1[e,e1];
            elseif h==2;
                ppr=ppbb2[e,e1];
            elseif h==3;
                ppr=ppbb3[e,e1];
            elseif h==4;
                ppr=ppbb4[e,e1];
            elseif h==5;
                ppr=ppbb5[e,e1];
            endif;
        endif;
    endif;
    ppr=ppr*ppz[s,s1]*ppef[h,h1];
    retp(ppr);
endp;



/* marginal utility */
proc(1)=uc(x);

  local c;
  c=x;
  if c<0;
    c=miss(1,1);  @ set c to Gauss missing value code if c<0 @
  else;
    c=c^(-sigma);
  endif;
 retp(c);
 
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
   tol=1e-10;
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

