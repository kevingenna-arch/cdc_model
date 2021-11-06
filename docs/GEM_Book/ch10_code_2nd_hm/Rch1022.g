@ ----------------------------------------- Rch1022.g -----------------------------


    last change: 25.10.2008

    author: Burkhard Heer

    computes the OLG model
    with the help of Krussell/Smith algorithm

    the model: OLG model with 70 generations, all equal measure
                *  2 stochastic individual productivities z[1], z[2] with 
                    transition probability probz
                *  2 permanent efficiency types eps[1],eps[2]
                *  2 stochastic aggregate productivities s[1], s[2] with 
                    transition probability probs
                *  endogenous labor
                *  lump-sum pensions

                *  aggregate state variables: K,s
                   individual state variables: a,z
                
                * algorithm: Krussell/Smith for OLG models, 
                  literature: see also Carroll and Young, mimeo

                * projected dynamics for K', N':
                    ln(K')= d[1] + d[2] * (s==1) + d[3]*ln(K) + d[4]* (s==1)*ln(K)
                    ln(N')= f[1] + f[2] * (s==1) + f[3]*ln(K') + f[4]* (s==1)*ln(K')

                * R^2 of the regression: practically 1

                * convergence: after approximately 10-12 hours on my Pentium 4
                
                * optimal policy functions etc are saved in order to simulate the 
                   economy and compute the correlations of the Ginis with income

results:
betak~d: 
     0.060964217      0.060546182 
     0.016128068      0.016159328 
      0.90620528       0.90696958 
   -0.0061691310    -0.0065087839 
betan~f: 
      -1.2694608       -1.2698890 
     0.023033914      0.023368544 
     -0.16890999      -0.16822869 
    0.0075210028     0.0071940122 

------------------------------------------------------------------------------------ @

new;
clear all;
library pgraph,nlsys;
nlset;
nint=100;
neg=-1e30;
Macheps=1e-30;

#include olg1.src;
#include ch8_toolbox.src;
GraphSettings;
_plwidth=7;
_MNR_Print=0;
_MNR_Global=1;
_MNR_QR=1;


cls;

firstrun=1; /* 0 - first run, saves steady state values, 
               1 - else */


/* parameters */
load efage=ef1;

/* periods */
nage=70;
Rage=46;    /* first period of retirement */
nr=nage-Rage+1; /* number of retirement years */
nw=Rage-1;     /* number of working years */
mass=ones(nage,1)/nage;     /* measure of all generations normalized to one, 
                                all generations have equal size */


/* final goods production */
alpha=0.35;
delta=0.08;

/* preferences */
sigma=2.0;
beta=0.99;
gam=0.28;

/* households types */
nz=2;   /* idiosyncratic shock, very persistent */
nj=2;   /* productivity types */
eps=(0.57|1.43);
probz=(0.98~0.02|0.02~0.98);
sigmaz=0.08;
z=zeros(nj,1);

x0=0.9|1.1;
bounds=(0~10|0~100);

    {crit,x}=FixvMN(x0,bounds,&prod1);
    z[1]=x[1];
    z[2]=x[2];
"efficiency types: ";
eps;
"stochastic productivityies: ";
z;
wait;
/* calibration: Mean efficiency=1 in steady state */

ef1=eps[1]*z[1]*efage;
ef2=eps[2]*z[1]*efage;
ef3=eps[1]*z[2]*efage;
ef4=eps[2]*z[2]*efage;

_plegctl=1;
_plegstr="j=1,z=0.727\000j=2,z=0.727\000j=1,z=1.27\000j=2,z=1.27";
_plctrl={0,10,15,0};
_pltype={6,6,6,3};
_pstype={0,2,8,0};
_psymsiz={0,3,3,0};

title("efficiency-age profile");
tw=seqa(20,1,nw);
xlabel("age");
ylabel("e(s,j,z)");
xy(tw,ef1~ef2~ef3~ef4);
wait;

/* aggregate shock */
/* a period is calibrated as one year; the expected duration of 
a business cycle is equal to six years */
s=(0.98|1.02);
probs=(2/3~1/3|1/3~2/3);
ns=2;

/* government */
zeta=0.3;     /* pension replacement ratio */
lbar0=0.3;      /* average working time */

/* initialization */
taubbar=(nage-rage+1)/nw*zeta;
taubold=taubbar;
bigl=0.2;       /* aggregate effective labor supply */
lbar=0.3;       /* average hourly labor supply */
lbarold=lbar;
rbar=0.02;
kbar=(alpha/(delta+rbar))^(1/(1-alpha));
"kbar: " kbar;
wbar=(1-alpha)*kbar^alpha;
"wbar: " wbar;
ybar=kbar^(alpha);
pen=zeta*lbar0*wbar;
wait;

/* initial values for procedures rfold, rfyoung */
r0=rbar; r1=rbar;
w0=wbar; 
taub0=taubbar; taub1=taubbar;
a0=0; it=0; iz=0; ia=0; ij=0;

@ --- computational parameters --- @

/* individual asset grid */
amin1=0;
amax1=15;
na=50;     /* grid on assets for value function */
agrid=seqa(amin1,(amax1-amin1)/(na-1),na);   /* individual asset grid */




psic=0.00001;       /* small constant */
psil=0.001;          /* minimum labor supply */
lmax=0.9;
maxit=50;          /* maximum number of iterations over steady state k */
tolk=0.001;         /* tolerance for equilibrium capital stock */
psi=0.7;            /* updating parameter */
kritold=100;

/* policy functions in steady state */
ass=zeros(nj*na*nz,nage);
css=zeros(nj*na*nz,nage);
lss=zeros(nj*na*nz,nw);
vss=zeros(nj*na*nz,nage);


@ ---------------------------------------------------------------------------

Step 1: computation of the steady state 

-------------------------------------------------------------------------- @

bsec=hsec;
krit=tolk+1;
kold=kbar;


if firstrun==0;

iq=0;
do until iq==maxit or krit<tolk;
    iq=iq+1;
    {a,c,l}=getvaluess();

    /* computation of aggregate employment and savings */
    ga=zeros(nj*na*nz,nage);
    bigl=0; /* total employment */
    bigc=0; /* total consumption */
    bigk=0; /* total assets k+i */
    lbar=0; /* average labor supply */
    
    /* initialization: assets at age 1 are 0 */
    ia=1; 
    ij=0;
    do until ij==nj;
    ij=ij+1;
    iz=0;
    do until iz==nz;
        iz=iz+1;
        ga[(ij-1)*na*nz+(iz-1)*na+ia,1]=mass[1]/(nz*nj);    
    endo;
    endo;

    it=0;
    do until it==nage-1;
        it=it+1;
        ia=0;
        do until ia==na;
            ia=ia+1;
            ij=0;
            do until ij==nj;
            ij=ij+1;
            iz=0;
            do until iz==nz;
                iz=iz+1;
                mass0=ga[(ij-1)*na*nz+(iz-1)*na+ia,it];
                a1=a[(ij-1)*na*nz+(iz-1)*na+ia,it];
                if a1<=amin1;
                    ia1=1;
                    iz1=0;
                    do until iz1==nz;
                        iz1=iz1+1;
                        ga[(ij-1)*na*nz+(iz1-1)*na+ia1,it+1]=ga[(ij-1)*na*nz+(iz1-1)*na+ia1,it+1]+probz[iz,iz1]*mass0;    
                    endo;
                elseif a1>=amax1;
                    ia1=na;
                    iz1=0;
                    do until iz1==nz;
                        iz1=iz1+1;
                        ga[(ij-1)*na*nz+(iz1-1)*na+ia1,it+1]=ga[(ij-1)*na*nz+(iz1-1)*na+ia1,it+1]+probz[iz,iz1]*mass0;    
                    endo;
                else;
                    ia1=minc(sumc(agrid.<a1)+1|na);
                    is1=(a1-agrid[ia1-1])/(agrid[ia1]-agrid[ia1-1]);
                    iz1=0;
                    do until iz1==nz;
                        iz1=iz1+1;
                        ga[(ij-1)*na*nz+(iz1-1)*na+ia1,it+1]=ga[(ij-1)*na*nz+(iz1-1)*na+ia1,it+1]+is1*probz[iz,iz1]*mass0;    
                        ga[(ij-1)*na*nz+(iz1-1)*na+ia1-1,it+1]=ga[(ij-1)*na*nz+(iz1-1)*na+ia1-1,it+1]+(1-is1)*probz[iz,iz1]*mass0;    
                    endo;
                endif;
            endo;   /* iz */
            endo;   /* ij */
        endo;
    endo;

    grossysum=0;
    ysum=0;    
    gag=zeros(nage,1);
    glg=zeros(nage,1);
    gyg=zeros(nage,1);
    fa=zeros(na,1);
    totalpen=sumc(mass[Rage:nage])*pen;
    totalwage=0;

    it=0;
    do until it==nage;
        it=it+1;
        ia=0;
        do until ia==na;
            ia=ia+1;
            ij=0;
            do until ij==nj;
            ij=ij+1;
            iz=0;
            do until iz==nz;
                iz=iz+1;
                mass0=ga[(ij-1)*na*nz+(iz-1)*na+ia,it];
                c1=c[(ij-1)*na*nz+(iz-1)*na+ia,it];
                if it>nw;
                    w0=pen*eps[ij];
                    w0=pen;
                else;
                    w0=(1-taubbar)*wbar*z[iz]*eps[ij]*efage[it]*l[(ij-1)*na*nz+(iz-1)*na+ia,it];
                endif;
                grossysum=grossysum+rbar*agrid[ia]*mass0;
                if it<=nw;
                    grossysum=grossysum+wbar*z[iz]*eps[ij]*efage[it]*l[(ij-1)*na*nz+(iz-1)*na+ia,it]*mass0;
                endif;                   
                fa[ia]=fa[ia]+mass0;
                bigk=bigk+mass0*agrid[ia];

                if it<=nw;
                    lbar=lbar+mass0*l[(ij-1)*na*nz+(iz-1)*na+ia,it];
                    bigl=bigl+mass0*l[(ij-1)*na*nz+(iz-1)*na+ia,it]*z[iz]*eps[ij]*efage[it];
                    totalwage=totalwage+mass0*l[(ij-1)*na*nz+(iz-1)*na+ia,it]*z[iz]*eps[ij]*efage[it]*wbar;
                endif;
                gag[it]=gag[it]+mass0*agrid[ia];
                gyg[it]=gyg[it]+(w0+rbar*agrid[ia])*mass0;
                if it<=nw;
                    glg[it]=glg[it]+mass0*l[(ij-1)*na*nz+(iz-1)*na+ia,it];
                endif;
            endo;   /* iz */
            endo;   /* ij */
        endo;
    endo;

    /* aggregation */
    lbar=lbar/sumc(mass[1:nw]);
    "gross income: " grossysum;
    "average labor supply: " lbar;
    "bigk: " bigk;
    "L: " bigl;
    "taub: " taubbar;
    xtemp=ef1|ef2|ef3|ef4;
    gtemp=mass[1:nw]|mass[1:nw]|mass[1:nw]|mass[1:nw];
    gtemp=gtemp/sumc(gtemp);
    giniwage=ginid(xtemp,gtemp);
    "gini wage distribution: " giniwage;
    giniwealth=ginid(agrid,fa);
    "gini wealth: " giniwealth;

    kbar=bigk/bigl;
    ybar=kbar^(alpha);

    rbarnew=alpha*kbar^(alpha-1)-delta;
    "rbarnew: " rbarnew;
    krit=abs((rbarnew-rbar)/maxc(1|rbar));
    rbar=psi*rbar+(1-psi)*rbarnew;
    wbarnew=(1-alpha)*kbar^alpha;    
    krit=krit|abs((wbarnew-wbar)/maxc(1|wbar));
    "wbarnew: " wbarnew;
    wbar=psi*wbar+(1-psi)*wbarnew;
    pen=zeta*lbar0*wbar;
    
    
    /* update of taub from the social security budget */
    taubnew=totalpen/totalwage;   
    "taubnew: " taubnew;  
    krit=krit|abs((taubnew-taubbar)/maxc(1|taubbar));
    taubbar=psi*taubbar+(1-psi)*taubnew;
        
    krit=maxc(krit);
    if kritold<krit;
        wait;
    endif;
    kritold=krit;
endo;



    ?etstr(hsec-bsec);  wait;
    title("distribution of capital");
    xy(agrid,fa);
    wait;
    title("age-wealth profile");
    xy(seqa(1,1,nage),gag./mass);
    wait;
    title("age-labor profile");
    xy(seqa(1,1,nage),glg./mass);
    wait;

    save ass,css,vss,lss;
    save bigl,bigk;
    save ga;

else;
    load ass,css,vss,lss;
    load bigl,bigk;
    load ga;
endif;


load nw,ef1,ef2,ef3,ef4,nage,gyg,mass,glg,gag;
"Figure 10.5: " wait;
GraphSettings;
_plwidth=7;
_ptitlht=0.2;
begwind;
window(2,2,0);
title("Capital stock");
xlabel("Age");
ylabel("");
xy(seqa(20,1,nage),gag./mass);
nextwind;
title("Labor supply");
xy(seqa(20,1,nw),glg[1:nw]./mass[1:nw]);
nextwind;
title("Total income");
xy(seqa(20,1,nage),gyg./mass);
nextwind;
title("Efficiency-age profiles");
_plegctl={1 6 40 0.0};
_ptitlht=0.2;
_plegstr="e_1,z_1\000e_2,z_1\000e_1,z_2\000e_2,z_2";
_plegstr="i=1,z=0.727\000i=2,z=0.727\000i=1,z=1.273\000i=2,z=1.273";
_plctrl={0,10,15,0};
_pltype={6,3,6,3};
_pstype={0,2,8,0};
_psymsiz={0,3,3,0};
xy(seqa(20,1,nw),ef1~ef2~ef3~ef4);
endwind;
save nw,ef1,ef2,ef3,ef4,nage,gyg,mass,glg,gag;
wait;

@ -------------------------------------------------------------------------------

Step 2: Computation of the dynamics with Krussell/Smith

------------------------------------------------------------------------------ @

/* aggregate capital stock */
kmin=0.8*bigk;
kmax=1.2*bigk;
nk=7;     /* grid */
kgrid=seqa(kmin,(kmax-kmin)/(nk-1),nk);  

/* aggregate labor */
nn=2;
nmin=0.9*bigl;
nmax=1.1*bigl;
ngrid=seqa(nmin,(nmax-nmin)/(nn-1),nn);   

/* income distribution */
/* net income: net wage income plus interest income */
ynetmin=0;
ynetmax=2;
ny=500;
ynetgrid=seqa(ynetmin,(ynetmax-ynetmin)/(ny-1),ny); 

/* total income: net wage income plus interest income plus pensions */
ytotmin=0;
ytotmax=2;
ytotgrid=seqa(ytotmin,(ytotmax-ytotmin)/(ny-1),ny); 

/* simulation parameter */
ndiscard=100;   /* discard the first ndiscard periods */
nsim=1000+ndiscard;  /* number of periods */
psi1=0.9;   /* update of d,f */
nq=20;      /* number of iterations over dynamics K', N' */
tolbeta=0.001;
kritbeta=1+tolbeta;

/* computation of the time series */
@

load ktsim, ntsim, bigk, bigl, nsim;
load loptnew,aoptnew,valuenew,coptnew;
load d,f,betakt,betant;
load ga;


@


/* policy functions over business cycle */
aopt=zeros(nj*na*nz*nk*ns,nage);  /* optimal next period equity */
copt=zeros(nj*na*nz*nk*ns,nage); /* optimal consumption */
lopt=zeros(nj*na*nz*nk*ns,nw);   /* optimal labor supply */
value=zeros(nj*na*nz*nk*ns,nage);   /* value function */

aoptnew=aopt;
coptnew=copt;
loptnew=lopt;
valuenew=value;


@ --------------------------------------------------------------------------

Step 2.1.: initialization of the dynamics for K' and N'

log(K') = d' * ( 1 z==1 log(K) (z==1)*log(K) )'

log(N') = f' * (1 z==1 K' (z==1)*log(K') )'

see: Carroll, Porapokkarm and Young

-------------------------------------------------------------------------- @

d=zeros(4,1);
d[3]=0.9;
d[1]=(1.0-d[3])*ln(bigk);

f=zeros(4,1);
f[3]=0;
f[1]=ln(bigl);   


@ -------------------------------------------------------------------------

Step 2.3: initialization of the optimal policy functions

-------------------------------------------------------------------------- @

ik=0;
do until ik==nk;
    ik=ik+1;
        is=0;
        do until is==ns;
            is=is+1;
            aopt[(ik-1)*ns*na*nz*nj+(is-1)*na*nz*nj+1:(ik-1)*ns*na*nz*nj+is*na*nz*nj,.]=ass;
            lopt[(ik-1)*ns*na*nz*nj+(is-1)*na*nz*nj+1:(ik-1)*ns*na*nz*nj+is*na*nz*nj,.]=lss;
            copt[(ik-1)*ns*na*nz*nj+(is-1)*na*nz*nj+1:(ik-1)*ns*na*nz*nj+is*na*nz*nj,.]=css;
            value[(ik-1)*ns*na*nz*nj+(is-1)*na*nz*nj+1:(ik-1)*ns*na*nz*nj+is*na*nz*nj,.]=vss;
        endo;    
endo;

aoptnew=aopt;
coptnew=copt;
loptnew=lopt;
valuenew=value;


/* initialization for rfold, rfyoung, rfyoung1, getvalue */
k0=bigk; k1=k0; 
n0=bigl; n1=n0;
s0=1; s1=1;
pen0=pen; pen1=pen;
z0=1; eps0=1;
w0=wbar; w1=wbar;
r0=rbar; r1=rbar;

betant=zeros(nq,4);
betakt=zeros(nq,4);

@ -------------------------------------------------------------------------

Step 2.3: computation of the optimal policy functions given d,f

-------------------------------------------------------------------------- @

iq=0;
bsec=hsec;
do while iq<nq and kritbeta>tolbeta; 
    iq=iq+1;
    
    if iq==5;
        nsim=3000; psi1=0.8;
    endif;


    if round(iq/5)==iq/5;
        cls;
    endif;

if firstrun<=1;

    {a,c,l}=getvalue();

    save aoptnew,loptnew;
else;
    load aoptnew,loptnew;
endif;


?etstr(hsec-bsec);  
if iq<1;
wait; 
endif;
cls;

sshock=rndu(nsim,1);
/* time series of aggregate variables K, N, S */
kt=zeros(nsim,1);   /* predicted value */
nt=zeros(nsim,1);
st=zeros(nsim,1);
ktsim=zeros(nsim,1); /* simulated values */
ntsim=zeros(nsim,1);
ytsim=zeros(nsim,1);
giniynetsim=zeros(nsim,1);
giniytotsim=zeros(nsim,1);
sharey20sim=zeros(nsim,1);
sharey40sim=zeros(nsim,1);
sharey60sim=zeros(nsim,1);
sharey80sim=zeros(nsim,1);
sharey95sim=zeros(nsim,1);
sharey100sim=zeros(nsim,1);

kt[1]=bigk; k0=bigk;
ktsim[1]=bigk;
nt[1]=bigl;
ntsim[1]=bigl;
if sshock[1]<0.5;
    st[1]=1; is=1;
else;
    st[1]=2; is=2;
endif;

ga0=ga; /* distribution in the first period: steady state distribution */

isim=1;
do until isim==nsim or (k0<kmin) or (k0>kmax);
    isim=isim+1;
    gaynet=zeros(ny,1);
    gaytot=zeros(ny,1);
    k0=ktsim[isim-1];
    is=st[isim-1];
    k1=f'*(1|(is==2)|ln(k0)|(is==2)*ln(k0));
    k1=exp(k1);
    kt[isim]=k1;
    if sshock[isim]<probs[is,1];
        is1=1;
    else;
        is1=2;
    endif;
    st[isim]=is1;
    n1=f'*(1|(is1==2)|ln(k1)|(is1==2)*ln(k1));
    n1=exp(n1); 
    nt[isim]=n1;
    n0=nt[isim-1];

    r0=rate(k0,n0,s[is]);
    w0=wage(k1,n1,s[is]);
    pen0=zeta*lbar0*w0;
    /* contribution rate in period 0: balanced social security budget */
    taub0=sumc(mass[Rage:nage])*pen0/(w0*n0);   

    /* dynamics of the distribution */
    ga1=zeros(nj*na*nz,nage);  /* next-period distribution */
    nsum=0;         /* total labor supply */
    bigasum=0;         /* aggregate wealth */

    /* initialization: assets at age 1 are 0 */
    ia=1;
    ij=0;
    do until ij==nj;
    ij=ij+1;
    iz=0;
    do until iz==nz;
        iz=iz+1;
        ga1[(ij-1)*na*nz+(iz-1)*na+ia,1]=mass[1]/(nz*nj);    
    endo;
    endo;

    it=0;
    do until it==nage-1;
        it=it+1;
        ia=0;
        do until ia==na;
            ia=ia+1;
            ij=0;
            do until ij==nj;
            ij=ij+1;
            iz=0;
            do until iz==nz;
                iz=iz+1;
                mass0=ga0[(ij-1)*na*nz+(iz-1)*na+ia,it];
                a1=alin(is,k0,ia,iz,ij,it);
                ynet0=r0*agrid[ia];
                if it<=nw;   /* labor supply in period isim-1 */
                    labor=llin(is,k0,ia,iz,ij,it);
                    nsum=nsum+mass0*labor*efage[it]*eps[ij]*z[iz];
                    ynet0=ynet0+labor*efage[it]*eps[ij]*z[iz]*w0;
                    ytot0=ynet0;
                else;
                    ytot0=ynet0+pen0;
                endif;

                iy1=minc(sumc(ynetgrid.<ynet0)+1|ny);
                iy0=(ynet0-ynetgrid[iy1-1])/(ynetgrid[iy1]-ynetgrid[iy1-1]);
                if iy1==1;
                    gaynet[1]=gaynet[1]+mass0;
                else;
                    gaynet[iy1-1]=gaynet[iy1-1]+(1-iy0)*mass0;
                    gaynet[iy1]=gaynet[iy1]+iy0*mass0;
                endif;
                
                iy1=minc(sumc(ytotgrid.<ytot0)+1|ny);
                iy0=(ytot0-ytotgrid[iy1-1])/(ytotgrid[iy1]-ytotgrid[iy1-1]);
                gaytot[iy1-1]=gaytot[iy1-1]+(1-iy0)*mass0;
                gaytot[iy1]=gaytot[iy1]+iy0*mass0;
                


                bigasum=bigasum+a1*mass0; 
                if a1<=amin1;
                    ia1=1;
                    iz1=0;
                    do until iz1==nz;
                        iz1=iz1+1;
                        ga1[(ij-1)*na*nz+(iz1-1)*na+ia1,it+1]=ga1[(ij-1)*na*nz+(iz1-1)*na+ia1,it+1]+probz[iz,iz1]*mass0;    
                    endo;
                elseif a1>=amax1;
                    ia1=na;
                    iz1=0;
                    do until iz1==nz;
                        iz1=iz1+1;
                        ga1[(ij-1)*na*nz+(iz1-1)*na+ia1,it+1]=ga1[(ij-1)*na*nz+(iz1-1)*na+ia1,it+1]+probz[iz,iz1]*mass0;    
                    endo;
                else;
                    ia1=minc(sumc(agrid.<a1)+1|na);
                    is1=(a1-agrid[ia1-1])/(agrid[ia1]-agrid[ia1-1]);
                    iz1=0;
                    do until iz1==nz;
                        iz1=iz1+1;
                        ga1[(ij-1)*na*nz+(iz1-1)*na+ia1,it+1]=ga1[(ij-1)*na*nz+(iz1-1)*na+ia1,it+1]+is1*probz[iz,iz1]*mass0;    
                        ga1[(ij-1)*na*nz+(iz1-1)*na+ia1-1,it+1]=ga1[(ij-1)*na*nz+(iz1-1)*na+ia1-1,it+1]+(1-is1)*probz[iz,iz1]*mass0;    
                    endo;
                endif;
            endo;   /* iz */
            endo;   /* ij */
        endo;   /* ia */
    endo;   /* it */
    
    it=nage;
        ia=0;
        do until ia==na;
            ia=ia+1;
            ij=0;
            do until ij==nj;
            ij=ij+1;
            iz=0;
            do until iz==nz;
                iz=iz+1;
                mass0=ga0[(ij-1)*na*nz+(iz-1)*na+ia,it];
                a1=alin(is,k0,ia,iz,ij,it);
                ynet0=r0*agrid[ia];
                if it<=nw;   /* labor supply in period isim-1 */
                    labor=llin(is,k0,ia,iz,ij,it);
                    nsum=nsum+mass0*labor*efage[it]*eps[ij]*z[iz];
                    ynet0=ynet0+labor*efage[it]*eps[ij]*z[iz]*w0;
                    ytot0=ynet0;
                else;
                    ytot0=ynet0+pen0;
                endif;

                iy1=minc(sumc(ynetgrid.<ynet0)+1|ny);
                iy0=(ynet0-ynetgrid[iy1-1])/(ynetgrid[iy1]-ynetgrid[iy1-1]);
                if iy1==1;
                    gaynet[1]=gaynet[1]+mass0;
                else;
                    gaynet[iy1-1]=gaynet[iy1-1]+(1-iy0)*mass0;
                    gaynet[iy1]=gaynet[iy1]+iy0*mass0;
                endif;
                
                iy1=minc(sumc(ytotgrid.<ytot0)+1|ny);
                iy0=(ytot0-ytotgrid[iy1-1])/(ytotgrid[iy1]-ytotgrid[iy1-1]);
                gaytot[iy1-1]=gaytot[iy1-1]+(1-iy0)*mass0;
                gaytot[iy1]=gaytot[iy1]+iy0*mass0;
            endo;
        endo;
    endo;                

    ktsim[isim]=bigasum;
    ntsim[isim-1]=nsum;
    ytsim[isim-1]=s[is]*ktsim[isim-1]^alpha*ntsim[isim-1]^(1-alpha);
    "iq~isim: " iq~isim;
    "K~N~bigk~bigl: " bigasum~nsum~bigk~bigl;

    y11=ginid(ynetgrid,gaynet);
    giniynetsim[isim-1]=y11;
    y11=ginid(ytotgrid,gaytot);
    giniytotsim[isim-1]=y11;
    gtotsum=cumsumc(gaytot);
    gtotsum1=ytotgrid.*gaytot;
    gtotsum1=cumsumc(gtotsum1);
    nshare20=sumc(gtotsum.<=0.2);
    sharey20sim[isim-1]=gtotsum1[nshare20]/gtotsum1[ny];
    nshare40=sumc(gtotsum.<=0.4);
    sharey40sim[isim-1]=(gtotsum1[nshare40]-gtotsum1[nshare20])/gtotsum1[ny];
    nshare60=sumc(gtotsum.<=0.6);
    sharey60sim[isim-1]=(gtotsum1[nshare60]-gtotsum1[nshare40])/gtotsum1[ny];
    nshare80=sumc(gtotsum.<=0.8);
    sharey80sim[isim-1]=(gtotsum1[nshare80]-gtotsum1[nshare60])/gtotsum1[ny];
    nshare95=sumc(gtotsum.<=0.95);
    sharey95sim[isim-1]=(gtotsum1[nshare95]-gtotsum1[nshare80])/gtotsum1[ny];
    sharey100sim[isim-1]=1-gtotsum1[nshare95]/gtotsum1[ny];


    ga0=ga1;

endo;   /* isim */


t1=seqa(1,1,isim-2);
ttsim=minc(isim-2|200);
t2=seqa(1,1,ttsim);
yt1=ytsim[t1];
incomet1=giniytotsim[t1];
incomet2=sharey20sim[t1];
incomet3=sharey40sim[t1];
incomet4=sharey60sim[t1];
incomet5=sharey80sim[t1];
incomet6=sharey95sim[t1];
incomet7=sharey100sim[t1];
{yt1}=hpfilter(ln(yt1),100);
"standard deviation of y: " stdc(yt1);
corin=corrx(incomet1~incomet2~incomet3~incomet4~incomet5~incomet6~incomet7~yt1);
"correlation of total gini, 20%, 40%, 60%, 80%, 95%, top 5 percentile income share";
"with income";
corin=corin[8,1:7];
corin;

if iq<1 or iq==nq;
wait;
title("simulted and predicted time series");
xlabel("time");
ylabel("capital stock");
xy(t1,bigk*ones(isim-2,1)~ktsim[t1]);
wait;

ylabel("employment");
xy(t1,bigl*ones(isim-2,1)~ntsim[t1]);
wait;

ylabel("simulated factor income shares");
xy(t1,sharey20sim[t1]~sharey40sim[t1]~sharey60sim[t1]~sharey80sim[t1]~sharey95sim[t1]~sharey100sim[t1]);
wait;

ylabel("simulated gini income");
xy(t1,giniynetsim[t1]~giniytotsim[t1]);
wait;
save giniytotsim,t1;

ylabel("simulated production");
xy(t1,ytsim[t1]);
wait;
save ytsim;

endif;

ndiscard1=round(isim/10);
nsim1=isim-ndiscard1-1;
xk1=st[ndiscard1+1:isim-1]-1;
xk1a=st[ndiscard1+1:isim-1];
xk2=ktsim[ndiscard1+1:isim-1];
xk3=xk1.*ln(xk2);
xk=ones(nsim1,1)~xk1~ln(xk2)~xk3;
y1=selif(xk2,xk1a);
y1=meanc(ln(y1));
y2=selif(xk2,xk1);
y2=meanc(ln(y2));

yk=ktsim[ndiscard1+2:isim];
yk=ln(yk);
yn=ntsim[ndiscard1+1:isim-1];
yn=ln(yn);
betak=inv(xk'*xk)*xk'*yk;
"y1~y2~bigk: " exp(y1)~exp(y2)~bigk;
"betak~d: " betak~d;


if iq<1;
wait;

"regression: capital stock";

 { vnam,m,b,stb,vc,stderr,sigma1,cx,rsq,resid,dwstat } = OLS(0,yk,xk);

wait;
endif;

betan=inv(xk'*xk)*xk'*yn;
"betan~f: " betan~f;
betant[iq,.]=betan';
betakt[iq,.]=betak';

kritbeta=maxc(abs(d-betak)|abs(f-betan));


if iq<4;
    /* update d[2] and d[4] and choose d[1] and d[3]
     so that mean capital stocks in simulation are equal to mean capital stock of dynamics */
    d[3]=psi1*d[3]+(1-psi1)*betak[3];
    d[4]=psi1*d[4]+(1-psi1)*betak[4];
/*
    d[1]=y1*(1-d[3]);
    d[2]=y2*(1-d[3]-d[4])-d[1];
*/
    d[1]=ln(bigk)*(1-d[3]);
    d[2]=ln(bigk)*(1-d[3]-d[4])-d[1];

else;
    d=psi1*d+(1-psi1)*betak;
endif;  


f=psi1*f+(1-psi1)*betan;

@
    d=psi1*d+(1-psi1)*betak;
@
"isim: " isim;
"new d~f: ";
"kritbeta:" kritbeta;
d~f;



if iq<1 or iq==nq;
wait;
"regression: employment";

 { vnam,m,b,stb,vc,stderr,sigma1,cx,rsq,resid,dwstat } = OLS(0,yn,xk);

wait;

endif;

aopt=aoptnew;
copt=coptnew;
lopt=loptnew;
value=valuenew;


endo;


t1=seqa(1,1,isim-2);
title("simulted and predicted time series");
xlabel("time");
ylabel("capital stock");
xy(t1,bigk*ones(isim-2,1)~ktsim[t1]);
wait;

ylabel("employment");
xy(t1,bigl*ones(isim-2,1)~ntsim[t1]);
wait;

save ktsim, ntsim, bigk, bigl, nsim;
save loptnew,aoptnew,valuenew,coptnew;
save d,f,betakt,betant;
save ga;
save giniytotsim,t1,ytsim;




@ --------------------------------------------------------------------------------


procedures


---------------------------------------------------------------------------------- @

proc utility(x,y);
    retp( ( (x^(gam).*y^(1-gam))^(1-sigma) )/(1-sigma));
endp;

proc uc(x,y);
    retp( gam*x^(gam*(1-sigma)-1).*y^((1-gam)*(1-sigma))  );
endp;

proc ul(x,y);
    retp( (1-gam)*x^(gam*(1-sigma)).*y^((1-gam)*(1-sigma)-1) );
endp;

proc wage(x,y,z);
    retp((1-alpha)*z*x^(alpha)*y^(-alpha));
endp;

proc rate(x,y,z);
    retp(alpha*z*x^(alpha-1)*y^(1-alpha)-delta);
endp;

proc(3)=getvaluess();
    local c, x1, aopt0, bounds, x0, amax2, wtemp, l;
    local rf1,a1test, crit;

    /* values for procedures rfold, rfyoung */
    r0=rbar; r1=rbar;
    w0=wbar; 
    taub0=taubbar; taub1=taubbar;
    

    ia=0;
    do until ia==na;
	    ia=ia+1;
	    a0=agrid[ia];
        ij=0;
        do until ij==nj;
        ij=ij+1;
        iz=0;
        do until iz==nz;
            iz=iz+1;
            c=(1+rbar)*a0+pen*eps[ij]; 
            c=(1+rbar)*a0+pen; 
	        ass[(ij-1)*na*nz+(iz-1)*na+ia,nage]=0;
	        css[(ij-1)*na*nz+(iz-1)*na+ia,nage]=c;
	        vss[(ij-1)*na*nz+(iz-1)*na+ia,nage]=utility(c,1);
        endo;   /* z */
        endo;   /* eps */
    endo; 	/* ia */

    /* computation of the decision rules for the retired */
    it=nage;
    do until it==nw+1;
        "iq~it~krit: " iq~it~krit;
        it=it-1;
        ia=0;
        do until ia==na;
	        ia=ia+1;
	        a0=agrid[ia];
            ij=0;
            do until ij==nj;
            ij=ij+1;
            iz=0;
            do until iz==nz;
                iz=iz+1;
                x1=rfoldss(0);

                if x1[1]>0;     /* randlösung? */
                    aopt0=0;
                else;
                    amax2=pen*eps[ij]+(1+r0)*a0-psic;   /* maximum capital stock for c=0; */
                    amax2=pen+(1+r0)*a0-psic;
                    /* bounds=(0~amax2); */
                    x0=amax2/2;
                    {aopt0,crit}=FixVMN1(x0,&rfoldss);
                endif;
                ass[(ij-1)*na*nz+(iz-1)*na+ia,it]=aopt0;
                c=pen*eps[ij]+(1+r0)*a0-aopt0;
                c=pen+(1+r0)*a0-aopt0;
                css[(ij-1)*na*nz+(iz-1)*na+ia,it]=c;
                vss[(ij-1)*na*nz+(iz-1)*na+ia,it]=utility(c,1);
            endo;   /* iz */
            endo;   /* ij */
        endo;
    endo;


    it=nw+1;
    do until it==1;
        it=it-1;
        "iq~it~krit: " iq~it~krit;

@
if it<nw;
    aopt1=zeros(na,nz);
    copt1=aopt1;    
    lopt1=aopt1;
        ia=0;
        do until ia==na;
	        ia=ia+1;
	        a0=agrid[ia];
            iz=0;
            do until iz==nz;
                iz=iz+1;
                aopt1[ia,iz]=ass[(iz-1)*na+ia,it+1];
                copt1[ia,iz]=css[(iz-1)*na+ia,it+1];
                lopt1[ia,iz]=lss[(iz-1)*na+ia,it+1];
            endo;
        endo;
    
        title("consumption");
        xy(agrid,copt1);
        wait;
        title("labor");
        xy(agrid,lopt1);
        wait;
        title("a'");
        xy(agrid,aopt1);
        wait;
endif;
@

        ia=0;
        do until ia==na;
	        ia=ia+1;
	        a0=agrid[ia];
            ij=0;
            do until ij==nj;
            ij=ij+1;
            iz=0;
            do until iz==nz;
                iz=iz+1;
                    c=gam*((1-taubbar)*z[iz]*eps[ij]*efage[it]*w0 + (1+r0)*a0);
                    l=1-(1-gam)/gam*c/((1-taubbar)*z[iz]*eps[ij]*efage[it]*w0);
                    if l<=0;    /* corner solution for l?*/
                        x1=rfyoung1ss(0);
                        l=0;
                        c=(1+r0)*a0;
                    else;
                        x1=rfyoungss(0);
                    endif;

                    
    
                    if x1[1]>0;     /* corner solution for a? */
                        aopt0=0;
                    else;      
                        wtemp=(1-taub0)*z[iz]*eps[ij]*efage[it]*w0;                  
                        amax2=(wtemp*lmax+(1+r0)*a0-psic);   /* c=0; */
                        /* bounds=(0~amax2); */

                        /* looking for a good initial value for a' */
                        x0=ass[(ij-1)*na*nz+(iz-1)*na+ia,it+1];
                        if ia>1;
                            x1=ass[(ij-1)*na*nz+(iz-1)*na+ia-1,it];
                            if abs(rfyoungss(x1))<abs(rfyoungss(x0));
                                x0=x1;
                            endif;
                        endif;

                        x1=ass[(ij-1)*na*nz+(iz-1)*na+ia,it];
                        if abs(rfyoungss(x1))<abs(rfyoungss(x0));
                            x0=x1;
                        endif;

                        do while rfyoungss(x0)==miss(1,1);
                            x0=x0/2;
                        endo;       
                        {aopt0,crit}=FixVMN1(x0,&rfyoungss);

                        if crit[1]>0 and iq==maxit;
                            "crit: "; crit;
                            wait;
                        endif;

                        c=gam*((1-taub0)*z[iz]*eps[ij]*efage[it]*w0 + (1+r0)*a0-aopt0);
                        l=1-(1-gam)/gam*c/((1-taub0)*z[iz]*eps[ij]*efage[it]*w0 );
                    
                        if l<=0;
                            l=0;
                            {aopt0,crit}=FixVMN1(x0,&rfyoung1ss);
                            c=(1+r0)*a0-aopt0;
                        endif;


                        if crit[1]>0 and iq>3;
                            "crit: "; crit;
                            wait;
                        endif;

                    endif;  /* a'=0? */

                    ass[(ij-1)*na*nz+(iz-1)*na+ia,it]=aopt0;
                    css[(ij-1)*na*nz+(iz-1)*na+ia,it]=c;
                    lss[(ij-1)*na*nz+(iz-1)*na+ia,it]=l;
                    vss[(ij-1)*na*nz+(iz-1)*na+ia,it]=utility(c,1-l);

                endo;   /* iz */
                endo;   /* ij */
            endo;
        endo;

    retp(ass,css,lss);
endp;


/* first-order condition for retired agent: steady state */
proc rfoldss(x);
    local a1,c0,rf,iz1,cd,c1,factor;
    a1=x;   /* optimal next-period assets */
    c0=pen*eps[ij]+(1+r0)*a0-a1;
    c0=pen+(1+r0)*a0-a1;
    if c0<=0;
        retp(miss(1,1));
    endif;
    
    rf=uc(c0,1);
    factor=beta*(1+r1);
    iz1=0;
    do until iz1==nz;
        iz1=iz1+1;
/* linear interpolation between two grid points of a */
        cd=css[(ij-1)*na*nz+(iz1-1)*na+1:(ij-1)*na*nz+iz1*na,it+1];
        c1=lininter(agrid,cd,a1);
        rf=rf-factor*probz[iz,iz1]*uc(c1,1);
    endo;
    retp(rf);
endp;

/* first order condition for young agent with l>0 */
proc rfyoungss(x);
    local a1,c0,l0,rf,iz1,cd,c1,l1,factor;
    a1=x;   /* optimal next-period assets */
    c0=gam*((1-taub0)*z[iz]*eps[ij]*efage[it]*w0 + (1+r0)*a0-a1);
    if c0<=0;
        retp(miss(1,1));
    endif;
    l0=1-(1-gam)/gam*c0/((1-taub0)*z[iz]*eps[ij]*efage[it]*w0);
/*
    "rfyoung: a0~a1~l0~c0: " a0~a1~l0~c0;
*/

    rf=uc(c0,1-l0);
    factor=beta*(1+r1);
    iz1=0;
    do until iz1==nz;
        iz1=iz1+1;
        cd=css[(ij-1)*na*nz+(iz1-1)*na+1:(ij-1)*na*nz+iz1*na,it+1];
        c1=lininter(agrid,cd,a1);
        if c1<=0; retp(miss(1,1)); endif;
        if it<nw;
            l1=1-(1-gam)/gam*c1/((1-taub1)*z[iz1]*eps[ij]*efage[it+1]*w0);
            if l1<0;
                l1=0;
            endif;
        else;
            l1=0;
        endif;
        rf=rf-factor*probz[iz,iz1]*uc(c1,1-l1);
    endo;

    retp(rf);
endp;

/* first order conditon for young agent with l=0 */
proc rfyoung1ss(x);
    local a1,c0,l0,rf,iz1,cd,c1,l1,factor;
    a1=x;   /* optimal next-period assets */
    c0=(1+r0)*a0-a1 ;
    if c0<=0;
        retp(miss(1,1));
    endif;

    l0=0;


    rf=uc(c0,1-l0);
    factor=beta*(1+r1);
    iz1=0;
    do until iz1==nz;
        iz1=iz1+1;
        cd=css[(ij-1)*na*nz+(iz1-1)*na+1:(ij-1)*na*nz+iz1*na,it+1];
        c1=lininter(agrid,cd,a1);
        if c1<=0; retp(miss(1,1)); endif;
        if it<nw;
            l1=1-(1-gam)/gam*c1/((1-taub1)*z[iz1]*eps[ij]*efage[it+1]*w0);
            if l1<0;
                l1=0;
            endif;
        else;
            l1=0;
        endif;
        rf=rf-factor*probz[iz,iz1]*uc(c1,1-l1);
    endo;

    retp(rf);
endp;

proc lininter(xd,yd,x);
   local j,nx,y;
   nx=rows(xd);
   j=sumc(xd.<=x);
   if j>=nx;    /* extrapolation */
      y= yd[nx]+(yd[nx]-yd[nx-1]).*(x-xd[nx])./(xd[nx]-xd[nx-1]);
   elseif j<=1;
      y=yd[1]+(yd[2]-yd[1]).*(x-xd[1])./(xd[2]-xd[1]);
   else;
      y=yd[j]+(yd[j+1]-yd[j]).*(x-xd[j])./(xd[j+1]-xd[j]);
   endif; 
   retp(y);
endp;


/* computes gini for discrete distribution of n households */
proc ginidn(x);
    local gini,y,n,i;
    n=rows(x);
    x=maxc(x'|zeros(1,n));
    y=sortc(x,1);
    y=y/sumc(y);
    gini=0;
    i=0;
    do until i==n;
        i=i+1;
        gini=gini+(2*i-n-1)/n*y[i];
    endo;
    retp(gini);
endp;



/* computes the gini for distribution where x has measure g(x) */
proc ginid(x,g);
    local xmean,ng,f,gini,i,y;	
    ng=rows(x);
    x=maxc(x'|zeros(1,ng));
	xmean=x'*g;
    y=(x~g);
    y=sortc(y,1);
    x=y[.,1];
    g=y[.,2];
    f=zeros(ng,1);  /* accumulated frequency */
    f[1]=g[1]*x[1]/xmean;
	gini=1-f[1]*g[1];
    i=1;
    do until i==ng;
        i=i+1;
        f[i]=f[i-1]+g[i]*x[i]/xmean;
        gini=gini-(f[i]+f[i-1])*g[i];
	endo;
    retp(gini);
endp;



proc(3)=getvalue();
    local c, x1, aopt0, bounds, x0, amax2, wtemp, l;
    local crit, count;


    ik=0;
    do until ik==nk;
        ik=ik+1;
        k0=kgrid[ik];
            is=0;
            do until is==ns;
                is=is+1;
                s0=s[is];  
                n0=f'*(1|(is==2)|ln(k0)|(is==2)*ln(k0));
                n0=exp(n0);
                w0=wage(k0,n0,s0);
                r0=rate(k0,n0,s0);
    
                pen0=zeta*lbar0*w0;
                /* contribution rate in period 0: balanced social security budget */
                taub0=sumc(mass[Rage:nage])*pen0/(w0*n0);   
                k1=d'*(1|(is==2)|ln(k0)|(is==2)*ln(k0));
                k1=exp(k1);
                
                "n0~k0~k1: " n0~k0~k1; (1|(is==2)|ln(k0)|(is==2)*ln(k0)); 


    /* last period optimization problem */                
                ia=0;
                do until ia==na;
	                ia=ia+1;
	                a0=agrid[ia];
                    ij=0;
                    do until ij==nj;
                    ij=ij+1;
                    iz=0;
                    do until iz==nz;
                        iz=iz+1;
                        "iq~ik~is~iz~ia: " iq~ik~is~iz~ia;
                        count=(ik-1)*ns*na*nz*nj+(is-1)*na*nz*nj+(ij-1)*na*nz+(iz-1)*na+ia;
                        c=(1+r0)*a0+pen0; 
	                    aoptnew[count,nage]=0;
	                    coptnew[count,nage]=c;
	                    valuenew[count,nage]=utility(c,1);
                    endo;   /* z */
                    endo;   /* eps */
                endo; 	/* ia */

    /* computation of the decision rules for the retired */
                it=nage;
                do until it==nw+1;
                    it=it-1;
                    ia=0;
                    do until ia==na;
	                    ia=ia+1;
	                    a0=agrid[ia];
                        ij=0;
                        do until ij==nj;
                        ij=ij+1;
                        iz=0;
                        do until iz==nz;
                            iz=iz+1;
                            count=(ik-1)*ns*na*nz*nj+(is-1)*na*nz*nj+(ij-1)*na*nz+(iz-1)*na+ia;
                            x1=rfold(0);
                            if x1[1]>0;     /* randlösung? */
                                aopt0=0;
                            else;
                                amax2=pen+(1+r0)*a0-psic;   /* maximum capital stock for c=0; */
                                x0=amax2/2;
                                {aopt0,crit}=FixVMN1(x0,&rfold);
                            endif;
                            aoptnew[count,it]=aopt0;
                            c=pen+(1+r0)*a0-aopt0;
                            coptnew[count,it]=c;
                            valuenew[count,it]=utility(c,1);
                        endo;
                        endo;
                    endo;
                endo;


                it=nw+1;
                do until it==1;
                    it=it-1;
                    ia=0;
                    do until ia==na;
	                    ia=ia+1;
	                    a0=agrid[ia];
                        iz=0;
                        do until iz==nz;
                            iz=iz+1;
                            z0=z[iz];
                            count=(ik-1)*ns*na*nz*nj+(is-1)*na*nz*nj+(ij-1)*na*nz+(iz-1)*na+ia;
                            c=gam*((1-taub0)*z0*eps[ij]*efage[it]*w0 + (1+r0)*a0);
                            l=1-(1-gam)/gam*c/((1-taub0)*z0*eps[ij]*efage[it]*w0);
                            if l<=0;    /* corner solution for l?*/
                                x1=rfyoung1(0);
                                l=0;
                                c=(1+r0)*a0;
                            else;
                                x1=rfyoung(0);
                            endif;

                            if x1[1]>0;     /* corner solution for a? */
                                aopt0=0;
                            else;      
                                wtemp=(1-taub0)*z[iz]*eps[ij]*efage[it]*w0;                  
                                amax2=(wtemp*lmax+(1+r0)*a0-psic);   /* c=0; */

                                /* looking for a good initial value for a' */
                                x0=aopt[count,it+1];
                                if ia>1;
                                    x1=aopt[count-1,it];
                                    if abs(rfyoung(x1))<abs(rfyoung(x0));
                                        x0=x1;
                                    endif;
                                endif;

                                x1=aopt[count,it];
                                if abs(rfyoung(x1))<abs(rfyoung(x0));
                                    x0=x1;
                                endif;

                                do while rfyoung(x0)==miss(1,1);
                                    x0=x0/2;
                                endo;       

                                {aopt0,crit}=FixVMN1(x0,&rfyoung);

                                if crit[1]>0 and iq==maxit; 
                                    "crit: "; crit;
                                    wait;
                                endif;

                                c=gam*((1-taub0)*z[iz]*eps[ij]*efage[it]*w0 + (1+r0)*a0-aopt0);
                                l=1-(1-gam)/gam*c/((1-taub0)*z[iz]*eps[ij]*efage[it]*w0 );
                    

                                if l<=0;
                                    l=0;
                                    {aopt0,crit}=FixVMN1(x0,&rfyoung1);
                                    c=(1+r0)*a0-aopt0;
                                endif;


                                if crit[1]>0 and iq>3;
                                    "crit: "; crit;
                                    wait;
                                endif;

                            endif;  /* a'=0? */

                            aoptnew[count,it]=aopt0;
                            coptnew[count,it]=c;
                            loptnew[count,it]=l;
                            valuenew[count,it]=utility(c,1-l);
                        endo;   /* endo iz */
                    endo;   /* endo ia */
                endo;   /* it */

            endo;   /* is */
    endo;   /* ik */

    retp(aoptnew,coptnew,loptnew);
endp;


/* first-order condition for retired agent: dynamics */
proc rfold(x);
    local a1,c0,rf,iz1,cd,c1,factor;
    local in1,is1,n1,s1,w1,r1,pen1;


    a1=x;   /* optimal next-period assets */
    c0=pen+(1+r0)*a0-a1;
    if c0<=0;
        retp(miss(1,1));
    endif;
    
    rf=uc(c0,1);
    is1=0;
    do until is1==ns;
        is1=is1+1;
        n1=f'*(1|(is1==2)|ln(k1)|(is1==2)*ln(k1));
        n1=exp(n1);
        s1=s[is1];  
        w1=wage(k1,n1,s1);
        r1=rate(k1,n1,s1);    
        pen1=zeta*lbar0*w1;
        factor=beta*(1+r1);
        iz1=0;
        do until iz1==nz;
            iz1=iz1+1;
/* bi-linear interpolation between two grid points of a,K */
            c1=cbil(is1,k1,a1,iz1);
            rf=rf-factor*probz[iz,iz1]*probs[is,is1]*uc(c1,1);
        endo;
    endo;
    retp(rf);
endp;


/* first order conditon for young agent with l=0 */
proc rfyoung1(x);
    local a1,c0,l0,rf,iz1,cd,c1,l1,factor,in1;
    local is1,n1,s1,w1,r1,pen1,taub1;

    a1=x;   /* optimal next-period assets */
    c0=(1+r0)*a0-a1 ;
    if c0<=0;
        retp(miss(1,1));
    endif;

    l0=0;
    rf=uc(c0,1-l0);

    is1=0;
    do until is1==ns;
        is1=is1+1;
        n1=f'*(1|(is1==2)|ln(k1)|(is1==2)*ln(k1));
        n1=exp(n1);
        s1=s[is1];  
        w1=wage(k1,n1,s1);
        r1=rate(k1,n1,s1);    
        pen1=zeta*lbar0*w1;
        taub1=sumc(mass[Rage:nage])*pen1/(w1*n1);   
        factor=beta*(1+r1);
        iz1=0;
        do until iz1==nz;
            iz1=iz1+1;
/* bi-linear interpolation between two grid points of a,K */
            c1=cbil(is1,k1,a1,iz1);

            if c1<=0; retp(miss(1,1)); endif;
            if it<nw;
                l1=1-(1-gam)/gam*c1/((1-taub1)*z[iz1]*eps[ij]*efage[it+1]*w1);
                if l1<0;
                    l1=0;
                endif;
            else;
                l1=0;
            endif;
            rf=rf-factor*probz[iz,iz1]*probs[is,is1]*uc(c1,1-l1);
        endo;
    endo;
    
    retp(rf);
endp;

/* first order condition for young agent with l>0 */
proc rfyoung(x);
    local a1,c0,l0,rf,is1,n1,s1,w1,r1,pen1,taub1,iz1,cd,c1,l1,factor;

    a1=x;   /* optimal next-period assets */
    c0=gam*((1-taub0)*z[iz]*eps[ij]*efage[it]*w0 + (1+r0)*a0-a1);
    if c0<=0;
        retp(miss(1,1));
    endif;
    l0=1-(1-gam)/gam*c0/((1-taub0)*z[iz]*eps[ij]*efage[it]*w0);


    rf=uc(c0,1-l0);

    is1=0;
    do until is1==ns;
        is1=is1+1;
        n1=f'*(1|(is1==2)|ln(k1)|(is1==2)*ln(k1));
        n1=exp(n1);
        s1=s[is1];  
        w1=wage(k1,n1,s1);
        r1=rate(k1,n1,s1);    
        pen1=zeta*lbar0*w1;
        taub1=sumc(mass[Rage:nage])*pen1/(w1*n1);   
        factor=beta*(1+r1);
        iz1=0;
        do until iz1==nz;
            iz1=iz1+1;
/* bi-linear interpolation between two grid points of a,K */
            c1=cbil(is1,k1,a1,iz1);
            if c1<=0; retp(miss(1,1)); endif;
            if it<nw;
                l1=1-(1-gam)/gam*c1/((1-taub1)*z[iz1]*eps[ij]*efage[it+1]*w1);
                if l1<0;
                    l1=0;
                endif;
            else;
                l1=0;
            endif;
            rf=rf-factor*probz[iz,iz1]*probs[is,is1]*uc(c1,1-l1);
        endo;
    endo;

    retp(rf);
endp;


/* bi-linear interpolation: if x or y are outside the grid, extrapolation */
proc bilint(xd,yd,zd,x,y);
   local nx,ny,z,z1,z2,i1,i2,zd1;
   nx=rows(xd);
   ny=rows(yd);
   i1=sumc(xd.<=x);
   i2=sumc(yd.<=y);
   i1=maxc(1|i1);
   i2=maxc(1|i2);
   i1=minc(nx-1|i1);
   i2=minc(ny-1|i2);

   z1=(x-xd[i1])/(xd[i1+1]-xd[i1]);
   z2=(y-yd[i2])/(yd[i2+1]-yd[i2]);
   z=(1-z1)*(1-z2)*zd[i1,i2];
   z=z+z1*(1-z2)*zd[i1+1,i2];
   z=z+z1*z2*zd[i1+1,i2+1];
   z=z+(1-z1)*z2*zd[i1,i2+1];

   retp(z);
endp;

/* next-period consumption - bilinear interpolation at K', a' */
proc cbil(is11,k11,a11,iz11);
    local c,xd,yd,zd,ik11,count1;
    xd=kgrid;
    yd=agrid;
    zd=zeros(na,nk);
    ik11=0;
    do until ik11==nk;
        ik11=ik11+1;
        count1=(ik11-1)*ns*na*nz*nj+(is11-1)*na*nz*nj+(ij-1)*na*nz+(iz11-1)*na;
        zd[.,ik11]=coptnew[count1+1:count1+na,it+1];
    endo;
    zd=zd';
    c=bilint(xd,yd,zd,k11,a11);
    retp(c);
endp;    
    


proc alin(is11,k11,ia11,iz11,ij11,it1);
    local a1,xd,yd,count1,ik11;
    xd=kgrid;
    yd=zeros(nk,1);
    ik11=0;
    do until ik11==nk;
        ik11=ik11+1;
        count1=(ik11-1)*ns*na*nz*nj+(is11-1)*na*nz*nj+(ij11-1)*na*nj+(iz11-1)*na+ia11;
        yd[ik11]=aoptnew[count1,it1];
    endo;
    a1=lininter(xd,yd,k11);
    retp(a1);
endp;    
    
proc llin(is11,k11,ia11,iz11,ij11,it1);
    local a1,xd,yd,count1,ik11;
    xd=kgrid;
    yd=zeros(nk,1);
    ik11=0;
    do until ik11==nk;
        ik11=ik11+1;
        count1=(ik11-1)*ns*na*nz*nj+(is11-1)*na*nz*nj+(ij11-1)*na*nj+(iz11-1)*na+ia11;
        yd[ik11]=loptnew[count1,it1];
    endo;
    a1=lininter(xd,yd,k11);
    retp(a1);
endp;    


proc prod1(x);
    local z1,z2,y,varz;
    z1=x[1];
    z2=x[2];
    y=x;
    y[1]=x[1]+x[2]-2;   /* normalization to one */
    /* variance of log productivity */
    varz=( (ln(z1))^2+(ln(z2))^2 ) / 2;    
    y[2]=sigmaz-varz;
    retp(y);
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

