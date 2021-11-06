@ ------------------------------ Rch8_unc.g -----------------

date: October 25, 2008

author: Burkhard Heer

aggregate uncertainty

PRIOR TO run THE PROGAM YOU NEED TO ADJUST THE save PATH!!!!!!!

computation of the dynamics in the heterogenous-agent neoclassical
growth model with value function iteration


-   interpolation between kbar (aggregate capital stock):
    linear
    lininter.g

-   interpolation between k (individual capital stock):
   linear
   lininter.g

-   interpolation between a and k (optimal next period capital
stock)
    bilinear
    bilina.g

-   Maximization: golden section search method
    golden.g


input procedures: 

golden -- routine for golden section search

lininter -- linear interpolation 
bilina -- bilinear interpolation of policy function bellman 
bellman -- right-hand side of bellman equation

-------------------------------------------------------------------------@


new; 
clear all; 
cls; 
library pgraph,user; 
GraphSettings;
_plwidth=7;
h0=hsec;

/*
save path="c:\\Documents and Settings\\bheer\\My Documents\\buch\\gauss\\RCh7\\";
*/


@ ----------------------------------------------------------------

Step 1: Parameterization and Initialization

----------------------------------------------------------------- @

tol=0.00001;    /* stopping criterion for final solution */
tol1=1e-7;      /* stopping criterion for golden section search */
eps=0.05; 
alpha=0.36; 
beta=0.96; 
delta=0.1; 
sigma=1.5; 
r=0.04;
w=0.5; 
tau=0.1; 
rep=0.25;

nit=100;    /* number of iterations over value function */

nq=30; /* number of iterations over aggrgate capital stock dynamics */

/* calibration of z */

zg=1.03; 
zb=0.97;

ppz=(0.8~0.2|
    0.2~0.8);

ppgg=(0.9615~0.0385|
        0.9581~0.0419); 

ppgg1=equivec1(ppgg);

ug=ppgg1[2];    /* unemployed in good times */


ppbb=(0.9525~0.0475|
        0.3952~0.6048); 
ppbb1=equivec1(ppbb);

ub=ppbb1[2];    /* unemployed in bad times */

ppgb=((1-ub)/(1-ug)~1-(1-ub)/(1-ug)|0~1);

ppbg=(1~0|1-ug/ub~ug/ub);

amin1=0;                 /* asset grid individual agents */
amax1=12; 
na=101; 
a=seqa(amin1,(amax1-amin1)/(na-1),na);

kmin1=2.5; 
kmax1=5.5;       /* asset grid aggregate capital stock */
nk=10; 
k=seqa(kmin1,(kmax1-kmin1)/(nk-1),nk);

/* simulation parameters */

nt=3000; /* number of transition periods */

nt1=500; /* stationary region nt1<nt */

nh=5000; /* number of households for simulation */

kind00=(kmax1-kmin1)/2+kmin1;   /* initial capital stock of all households */

z0=zg; /* initial technology level */

/* initialization of the value function: */

veg=zeros(na,nk);    /* employed agents, good state */
veb=zeros(na,nk);    /* employed agents, bad state */
vug=zeros(na,nk);    /* unemployed agents, good state */
vub=zeros(na,nk);    /* unemployed agents, bad state */

i=0; 
do until i==nk;
    i=i+1;
    b=rep*(1-tau)*w;
    veg[.,i]=u((1-tau)*r*a+(1-tau)*w);
    veb[.,i]=u((1-tau)*r*a+(1-tau)*w);
    vug[.,i]=u((1-tau)*r*a+b);
    vub[.,i]=u((1-tau)*r*a+b);
endo;

veg=veg/(1-beta);                /* agents consume their income */
veb=veb/(1-beta); 
vug=vug/(1-beta); 
vub=vub/(1-beta);

copteg=zeros(na,nk);           /* optimal consumption, z=zg */
copteb=zeros(na,nk); 
coptug=zeros(na,nk); 
coptub=zeros(na,nk);

aopteg=zeros(na,nk);        /* optimal next-period assets, z=zg */
aopteb=zeros(na,nk);
aoptug=zeros(na,nk); /* optimal next-periodassets, z=zb */
aoptub=zeros(na,nk);


/* law of motion for capital: initial guess,
ln(kbar')=gam0+gam1*ln(kbar) */

gam0g=0.17; 
gam1g=0.89; 
gam0b=0.13; 
gam1b=0.90;

phi1=0.5;   /* updating of gam0,gam1 */

kritgam=1; 
tolgam=0.001;

gam0q=zeros(nq,2); 
gam1q=zeros(nq,2);


q=0;
do until q==nq; /* or kritgam<tolgam; */
    q=q+1;

    if q==5 or q==10 or q==15 or q==20 or q==25 or q==30 or q==35;
        cls;
    endif;

    @ ----------------------------------------------------------------

    Step 4: Iteration of the value function over K and a

    ----------------------------------------------------------------- @

    crit=tol+1;
    neg=-1e10;
    j=0;

    do until (crit<tol and j>50) or(j==nit);
        j=j+1;
    
        "iteration j~q: " j~q; 
        "time elapsed: " etstr(hsec-h0); 
        "error: " crit;

        voldeg=veg; 
        voldeb=veb; 
        voldug=vug; 
        voldub=vub;

        s=0;
        do until s==2;      /* iteration over z=z_g (s=1) and z=z_b (s=2) */
            s=s+1;
            if s==1;
                zz=zg; 
                gam0=gam0g; 
                gam1=gam1g;  
                nn0=1-ug;
            else;
                zz=zb; 
                gam0=gam0b; 
                gam1=gam1g; 
                nn0=1-ub;
            endif;

            m=0; 
            do until m==nk;     /* iteration over aggregate capital stock k */
                m=m+1;
                r=zz*rrate(k[m],nn0);
                w=zz*wage(k[m],nn0);
                tau=w*(1-nn0)*rep/(w*nn0+r*k[m]+w*(1-nn0)*rep);
                b=(1-tau)*rep*w;
                k1=exp(gam0+gam1*ln(k[m]));  /* next period aggregate capital stock */
                veg1=zeros(na,1);    /* value function of employed with k'=k1, z'=zg */
                veb1=zeros(na,1);
                vug1=zeros(na,1);    /* value function of unemployed with k'=k1, z'=zg */
                vub1=zeros(na,1);
                i=0;
                do until i==na;     /* iteration over next-period individual wealth a' */
                    i=i+1;
                    if k1<=kmin1;
                        veg1[i]=veg[i,1];
                        veb1[i]=veb[i,1];
                        vug1[i]=vug[i,1];
                        vub1[i]=vub[i,1];
                    elseif k1>=kmax1;
                        veg1[i]=veg[i,nk];
                        veb1[i]=veb[i,nk];
                        vug1[i]=vug[i,nk];
                        vub1[i]=vub[i,nk];
                    else;
                        veg1[i]=lininter(k,veg[i,.]',k1);
                        veb1[i]=lininter(k,veb[i,.]',k1);
                        vug1[i]=lininter(k,vug[i,.]',k1);
                        vub1[i]=lininter(k,vub[i,.]',k1);
                    endif;
                endo;

                e=0;    /* iteration over the employment status */
                do until e==2;  /* e=1 employed, e=2 unemployed */
                    e=e+1;
                    i=0;       
                    l0=0;
                    do until i==na;  /* iteration over asset grid a in period t */
                        i=i+1;
                        l=l0;
                        v0=neg;
                        ax=amin1; bx=amin1; cx=amax1;
                        do until l==na; /* iteration over a' in period t*1 */
                            l=l+1;
                            if e==1;
                                c=(1+(1-tau)*r)*a[i]+(1-tau)*w-a[l];
                            else;
                                c=(1+(1-tau)*r)*a[i]+b-a[l];
                            endif;
                            if c>0;
                                v1=bellman(a[i],a[l],e);
                                if v1>v0;
                                    if e==1;
                                        if zz==zg;
                                            veg[i,m]=v1;
                                        else;
                                            veb[i,m]=v1;
                                        endif;
                                    else;
                                        if zz==zg;
                                            vug[i,m]=v1;
                                        else;
                                            vub[i,m]=v1;
                                        endif;
                                endif;

                                if l==1;
                                   ax=a[1]; bx=a[1]; cx=a[2];
                                elseif l==na;
                                    ax=a[na-1]; bx=a[na]; cx=a[na];
                                else;
                                    ax=a[l-1]; bx=a[l]; cx=a[l+1];
                                endif;

                                v0=v1;
                                l0=maxc(0|l-2);
                            else;
                                l=na;   /* concavity of value function */
                            endif;
                        else;
                            l=na;       /* c<0 */
                        endif;
                    endo;   /* l=1,..,na */

                    if ax==bx;  /* boundary optimum, ax=bx=a[1]  */
                        bx=ax+eps*(a[2]-a[1]);
                        if bellman(a[i],bx,e)<bellman(a[i],ax,e);
                            if e==1;
                                if zz==zg;
                                    aopteg[i,m]=a[1];
                                else;
                                    aopteb[i,m]=a[1];
                                endif;
                            else;
                                if zz==zg;
                                    aoptug[i,m]=a[1];
                                else;
                                    aoptub[i,m]=a[1];
                                endif;
                            endif;
                        else;
                            if e==1;
                                if zz==zg;
                                    aopteg[i,m]=golden(&value1,ax,bx,cx,tol1);
                                else;
                                    aopteb[i,m]=golden(&value1,ax,bx,cx,tol1);
                                endif;
                            else;
                                if zz==zg;
                                    aoptug[i,m]=golden(&value1,ax,bx,cx,tol1);
                                else;
                                    aoptub[i,m]=golden(&value1,ax,bx,cx,tol1);
                                endif;
                            endif;
                        endif;
                    elseif bx==cx;  /* boundary optimum, bx=cx=a[n] */
                        bx=cx-eps*(a[na]-a[na-1]);
                        if bellman(a[i],bx,e)<bellman(a[i],cx,e);
                            if e==1;
                                if zz==zg;
                                    aopteg[i,m]=a[na];
                                else;
                                    aopteb[i,m]=a[na];
                                endif;
                            else;
                                if zz==zg;
                                    aoptug[i,m]=a[na];
                                else;
                                    aoptub[i,m]=a[na];
                                endif;
                            endif;
                        else;
                            if e==1;
                                if zz==zg;
                                    aopteg[i,m]=golden(&value1,ax,bx,cx,tol1);
                                else;
                                    aopteb[i,m]=golden(&value1,ax,bx,cx,tol1);
                                endif;
                            else;
                                if zz==zg;
                                    aoptug[i,m]=golden(&value1,ax,bx,cx,tol1);
                                else;
                                    aoptub[i,m]=golden(&value1,ax,bx,cx,tol1);
                                endif;
                            endif;
                        endif;
                    else;
                        if e==1;
                            if zz==zg;
                                aopteg[i,m]=golden(&value1,ax,bx,cx,tol1);
                            else;
                                aopteb[i,m]=golden(&value1,ax,bx,cx,tol1);
                            endif;
                        else;
                            if zz==zg;
                                aoptug[i,m]=golden(&value1,ax,bx,cx,tol1);
                        else;
                            aoptub[i,m]=golden(&value1,ax,bx,cx,tol1);
                        endif;
                    endif;
                endif;

               if e==1;
                    if zz==zg;
                        veg[i,m]=bellman(a[i],aopteg[i,m],e);
                    else;
                        veb[i,m]=bellman(a[i],aopteg[i,m],e);
                    endif;
                else;
                    if zz==zg;
                        vug[i,m]=bellman(a[i],aoptug[i,m],e);
                    else;
                        vub[i,m]=bellman(a[i],aoptug[i,m],e);
                    endif;
               endif;
            
            endo;   /* i=1,..na */
        endo;   /* e=1,2 */        
    endo;   /* m=1..nk */

    if zz==zg;
        copteg=(1+(1-tau)*r)*a+(1-tau)*w-aopteg;
        coptug=(1+(1-tau)*r)*a+b-aoptug;
    else;
        copteb=(1+(1-tau)*r)*a+(1-tau)*w-aopteb;
        coptub=(1+(1-tau)*r)*a+b-aoptub;
    endif;

    endo;   /* s=1,2 -- zz=zg,zb */

crit=meanc(abs(voldeg-veg)|abs(voldeb-veb)|abs(voldug-vug)|abs(voldub-vub));

"iteration q~j " q~j; 
"time elapsed: " etstr(hsec-h0); "error
value function: " crit;

endo;   /* j=1,..nit */

if q==nq;
    title("optimal next-period asset of the employed agent in good times");
    ylabel("a(e,a,zg,K)");
    xy(a,aopteg[.,1:nk]-a);
    wait;

    title("optimal next-period asset of the employed agent in bad times");
    ylabel("a(e,a,zb,K)");
    xy(a,aopteb[.,1:nk]-a);
    wait;
    
    save aunc=a,aoptebunc=aopteb,aoptubunc=aoptub,nkunc=nk;

    title("optimal next-period asset of the unemployed agent in good times");
    ylabel("a(u,a,zg,K)");
    xy(a,aoptug[.,1:nk]-a);
    wait;

    title("optimal next-period asset of the unemployed agent in bad times");
    ylabel("a(u,a,zb,K)");
    xy(a,aoptub[.,1:nk]-a);
    wait;


    title("optimal consumption of the employed agent in good times");
    ylabel("c(e,a,zg,K)");
    xy(a,copteg[.,1:nk]);
    wait;

    title("optimal consumption of the employed agent in bad times");
    ylabel("c(e,a,zb,K)");
    xy(a,copteb[.,1:nk]);
    wait;

    title("optimal consumption of the unemployed agent in good times");
    ylabel("c(u,a,zg,K)");
    xy(a,coptug[.,1:nk]);
    wait;

    title("optimal consumption of the unemployed agent in bad times");
    ylabel("c(u,a,zb,K)");
    xy(a,coptub[.,1:nk]);
    wait;
endif;

    save kcopteg=copteg,kcoptug=coptug,kaopteg=aopteg;

    save kaoptug=aoptug,kveg=veg,kvug=vug;

    save kcopteb=copteb,kcoptub=coptub,kaopteb=aopteb;

    save kaoptub=aoptub,kveb=veb,kvub=vub; 
    save kuc=k,auc=a;

    @ ----------------------------------------------------------------

    Step 3: Simulation

    ----------------------------------------------------------------- @


    kt=zeros(nt,2); /* aggregate capital stock K and technology z over time */
    ct=zeros(nt,1);
    yt=zeros(nt,1);


    kind0=zeros(nh,2); /* individual capital stock k and employment e in t-1 */

    kind1=zeros(nh,2); /* k,e in t */

    kpred=zeros(nt-1,3);  /* predicted versus simulated capital stock */

    if z0==zg;
        j0=round(ug*nh);    /* equilibrium number of unemployed */
    else;
        j0=round(ub*nh);
    endif;

    kind0[1:j0,.]=kind00*ones(j0,1)~2*ones(j0,1);
    kind0[j0+1:nh,.]=kind00*ones(nh-j0,1)~ones(nh-j0,1);
    kt1=meanc(kind0[.,1]); kt[1,.]=kt1~z0;
    yt[1]=kind00^alpha*(1-ug)^(1-alpha);
    ct[1]=yt[1]-delta*kt[1,1];

    t=1;

    do until t==nt; /* transition of distribution function */

        x0=rndu(1,1);
        t=t+1; 
        "q~t~k(t-1): " q~t~kt[t-1,.]; 
        "gam0~gam1: ";
        gam0g~gam0b~gam1g~gam1b;
        if kt[t-1,2]==zg; 
            j0=1; 
        else; 
            j0=2; 
        endif;
        if x0<ppz[j0,1];
            kt[t,2]=zg;
        else;
            kt[t,2]=zb;
        endif;
        i=0;
        xs=rndu(nh,1);
        do until i==nh;
            i=i+1;
            kind1[i,1]=bilina(kind0[i,1],kt[t-1,1],kind0[i,2],kt[t-1,2]);
            if kt[t-1,2]==zg and kt[t,2]==zg;
                if xs[i]<ppgg[kind0[i,2],1];
                    kind1[i,2]=1;
                else;
                    kind1[i,2]=2;
                endif;
            elseif kt[t-1,2]==zg and kt[t,2]==zb;
                if xs[i]<ppgb[kind0[i,2],1];
                    kind1[i,2]=1;
                else;
                    kind1[i,2]=2;
                endif;
            elseif kt[t-1,2]==zb and kt[t,2]==zg;
                if xs[i]<ppbg[kind0[i,2],1];
                    kind1[i,2]=1;
                else;
                    kind1[i,2]=2;
                endif;
            else;
                if xs[i]<ppbb[kind0[i,2],1];
                    kind1[i,2]=1;
                else;
                    kind1[i,2]=2;
                endif;
            endif;
        endo; /* i=1,..,nh */

         /* test of u=ug if z=zg or u=ub if z=zb */
         kind1=sortc(kind1,2);
         if kt[t,2]==zg; 
            nu=round(nh*ug); 
        else; 
            nu=round(nh*ub); 
        endif;
         x0=sumc(kind1[.,2].<=1);
         if x0<(nh-nu);
            kind1[1:nh-nu,2]=ones(nh-nu,1);
         endif;
         if x0>(nh-nu);
            kind1[nh-nu+1:nh,2]=2*ones(nu,1);
         endif;
         kt[t,1]=meanc(kind1[.,1]);
        if kt[t,2]==1;
            yt[t]=kt[t,1]^alpha*(1-ug)^(1-alpha);
        else;
            yt[t]=kt[t,1]^alpha*(1-ub)^(1-alpha);
        endif;
        ct[t]=yt[t]-delta*kt[t-1,1]-(kt[t,1]-kt[t-1,1]);
         kind0=kind1;

            if kt[t-1,2]==1;
                gam0=gam0g; 
                gam1=gam1g;  
            else;
                gam0=gam0b; 
                gam1=gam1g; 
            endif;
            kpred[t-1,2]=exp(gam0+gam1*ln(kt[t-1,1]));
            kpred[t-1,1]=kt[t,1]; kpred[t-1,3]=kt[t-1,1];

     endo;  /* t=1,..,nt */


    /* transition dynamics for z=zg and z=zb */
   ktx=kt[nt1:nt-1,.]~kt[nt1+1:nt,1];
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

    kind00=meanc(kt[.,1]);

    save kkt=kt;
    kritgam=abs(ln(gam0g)-ln(gam01g)|ln(gam1g)-ln(gam11g)|
    ln(gam0b)-ln(gam01b)|ln(gam1b)-ln(gam11b));
    "kritgam: " kritgam';
    gam0g=phi1*gam0g+(1-phi1)*gam01g;
    gam1g=phi1*gam1g+(1-phi1)*gam11g;
    gam0b=phi1*gam0b+(1-phi1)*gam01b;
    gam1b=phi1*gam1b+(1-phi1)*gam11b;
    save kgam0g=gam0g,kgam1g=gam1g; 
    save kgam0b=gam0b,kgam1b=gam1b;
    gam0q[q,.]=gam0g~gam0b;
    gam1q[q,.]=gam1g~gam1b;
    save gam0q,gam1q;
endo;   /* q=1,..,nq */

        "iteration complete: q " q;
        "gam0g~gam1g: " gam0g~gam1g;
        "gam0b~gam1b: " gam0b~gam1b;
        "meanc(k): " meanc(kt[.,1]);
        "mean prediction error, %" meanc(abs( (kpred[.,2]-kpred[.,1])./kpred[.,1] ) );
        wait;

        "Figure 8.9: " wait;
        ylabel("");
        title("");
        dista=sortc(kind1,1);
        xlabel("Individual wealth a");
        xy(dista[.,1],seqa(1/nh,1/nh,nh));
        wait;
        save dista,nh;
        xlabel("time");
        title("capital stock");
        xy(seqa(1,1,nt),kt[.,1]);
        wait;
        save nt,kt;
        xlabel("capital stock in t-1");
        title("prediction error");
        xy(kpred[.,3],kpred[.,1]-kpred[.,2]);
        wait;
        "standard deviation k~y~c: ";
        stdc(kt[nt1:nt,1]);
        stdc(yt[nt1:nt]);
        stdc(ct[nt1:nt]);
        wait;
        "hp-filtered standard deviation k~y~c";
{x0,ct1}=hptrend(ln(ct[200:1000]),100);
sc=stdc(ct1);
{x0,yt1}=hptrend(ln(yt[200:1000]),100);
sy=stdc(yt1);
{x0,kt1}=hptrend(ln(kt[200:1000,1]),100);
sk=stdc(kt1);
    sk~sy~sc;

@  ----------------------------  procedures -----------


u(x) -- utility function

value(a,e) -- returns the value of the value function for asset
                a and employment status e

value1(x) -- given a=a[i] and epsilon=e, returns the
             value of the bellman equation for a'=x

bellman -- value for the right-hand side of the Bellman equation

bilina -- bilinear interpolation of a'(a,K)

------------------------------------------------------- @

proc u(x);
   retp(x^(1-sigma)/(1-sigma));
endp;

proc rrate(x,y);
    retp(alpha*(y/x)^(1-alpha)-delta);
endp;

proc wage(x,y);
    retp((1-alpha)*(x/y)^(alpha));
endp;

/* value function in period t+1 */
/* x - a_t+1, y - epsilon_t+1, z - z_t+1: zg,zb */
proc value(x,y,z);
    if y==1;
        if z==zg;
            retp(lininter(a,veg1,x));
        else;
            retp(lininter(a,veb1,x));
        endif;
    else;
        if z==zg;
            retp(lininter(a,vug1,x));
        else;
            retp(lininter(a,vub1,x));
        endif;
    endif;
endp;

proc value1(x);
    retp(bellman(a[i],x,e));
endp;


proc bellman(a0,a1,y);
   local c,k1,zeta0,zeta1;
   if y==1;
      c=(1+(1-tau)*r)*a0+(1-tau)*w-a1;
   else;
      c=(1+(1-tau)*r)*a0+b-a1;
   endif;
   if c<0;
      retp(neg);
   endif;
   if a1>=a[na];
    if zz==zg;
        zeta0=ppgg[y,1]*veg1[na]+ppgg[y,2]*vug1[na];
        zeta1=ppgb[y,1]*veb1[na]+ppgb[y,2]*vub1[na];
        retp(u(c)+beta*(ppz[1,1]*zeta0+ppz[1,2]*zeta1));
    else;
        zeta0=ppbg[y,1]*veg1[na]+ppbg[y,2]*vug1[na];
        zeta1=ppbb[y,1]*veb1[na]+ppbb[y,2]*vub1[na];
        retp(u(c)+beta*(ppz[2,1]*zeta0+ppz[2,2]*zeta1));
    endif;
   endif;
   if a1==a[1];
       if zz==zg;
            zeta0=ppgg[y,1]*veg1[1]+ppgg[y,2]*vug1[1];
            zeta1=ppgb[y,1]*veb1[1]+ppgb[y,2]*vub1[1];
            retp(u(c)+beta*(ppz[1,1]*zeta0+ppz[1,2]*zeta1));
        else;
            zeta0=ppbg[y,1]*veg1[1]+ppbg[y,2]*vug1[1];
            zeta1=ppbb[y,1]*veb1[1]+ppbb[y,2]*vub1[1];
            retp(u(c)+beta*(ppz[2,1]*zeta0+ppz[2,2]*zeta1));
        endif;
   endif;
   if zz==zg;
        zeta0=ppgg[y,1]*value(a1,1,zg)+ppgg[y,2]*value(a1,2,zg);
        zeta1=ppgb[y,1]*value(a1,1,zb)+ppgb[y,2]*value(a1,2,zb);
        retp(u(c)+beta*(ppz[1,1]*zeta0+ppz[1,2]*zeta1));
    else;
        zeta0=ppbg[y,1]*value(a1,1,zg)+ppbg[y,2]*value(a1,2,zg);
        zeta1=ppbb[y,1]*value(a1,1,zb)+ppbb[y,2]*value(a1,2,zb);
        retp(u(c)+beta*(ppz[2,1]*zeta0+ppz[2,2]*zeta1));
    endif;
endp;


proc bilina(a0,k0,e0,z0); /* bilinear interpolation */
    local n0,n1,n2,m0,m1,m2,phi,f,d,z1,z2,y,pol;
    if e0==1;
        if z0==zg;
            pol=aopteg;
        else;
            pol=aopteb;
        endif;
    else;
        if z0==zg;
            pol=aoptug;
        else;
            pol=aoptub;
        endif;
    endif;
    if a0<amin1; a0=amin1; endif;
    if a0>amax1; a0=amax1; endif;
    if k0<kmin1; k0=kmin1;  endif;
    if k0>kmax1; k0=kmax1; endif;
    m0=(a0-amin1)/(amax1-amin1)*(na-1)+1;
    m2=floor(m0);
    m1=m0-m2;
    n0=(k0-kmin1)/(kmax1-kmin1)*(nk-1)+1;
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


proc lininter(xd,yd,x);
  local j;
  j=sumc(xd.<=x');
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

