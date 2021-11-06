@ ------------------------------ Rch92AK6.g ----------------------------

    Burkhard Heer, 20.5.2008

    algorithm 9.2, solves the Auerbach-Kotlikoff 6-period model, transition dynamics

    unexpected change in the replacement ratio of pensions from 0.3 to 0.2

    direct computation of the OLG model in section 9.1

    transition periods: nt > 6

    REVISED VERSION WITHOUT NLSYS

    THE USER CAN CHOOSE BETWEEN THE FOLLOWING UPDATING METHODS IN LINE 33:

    1 -- LINEAR UPDATE
    2 -- NEWTON'S METHOD
    3 -- BROYDEN, with initialization of Broyden matrix with Jacobi matrix
    4 -- BROYDEN, initialization of Broyden matrix with derivertive at final steady state

-------------------------------------------------------------------------@


new;

clear all;
cls;
library pgraph,user;
MachEps=1e-10;
GraphSettings;
_plwidth=7;
_pxsci=5;

_update=4;  /* 1 -- linear update, 2 -- Newton,  3 -- Broyden, with Jacobi matrix initialization
              4 --- Broyden, initialization with derivate at final steady state */
phi=0.8;    /* updating parameter for method 1 */

t=4;
tr=2;

bsec=hsec;

@ ---------------------------

 parameter

---------------------------- @

beta=0.9;         /* discount factor */
r=0.2;         /* initial value of the interest rate */
s=2;            /* coefficient of relative risk aversion */
alp=0.3;        /* production elasticity of capital */
rep0=0.3;        /* replacement ratio */
rep1=0.2;
del=0.4;          /* rate of depreciation */
tau=rep0/(2+rep0);    /* income tax rate */
gam=2;            /* disutility from working */

kmax=1;        /* upper limit of capital grid */
kinit=0;
na=101;          /* number of grid points on assets */
a=seqa(0,kmax/(na-1),na);   /* asset grid */
psi=0.001;      /* parameter of utility function */
tol=0.001;       /* percentage deviation of final solution */
tolk=0.0001;       /* percentage deviation of final solution for k_1 */
nq1=30;             /* number of iterations over k^1 */
nt=20;         /* number of transition periods */
nqt=20;         /* maximum number of iterations over transition of K_t,N_t */
tolt=0.0001;     /* tolerance with regard to K_t,N_t */

@ ---------------------------------

initialization

--------------------------------- @

nbar=0.2;
kbar=(alp/(r+del))^(1/(1-alp))*nbar;

kold=100;
nold=2;


/* agents' policy function in steady state */
aopt=zeros(6,1);  /* optimal asset */
copt=zeros(6,1);   /* optimal consumption */
nopt=0.3*ones(4,1);    /* optimal labor supply */
w=0;
r=0;
pen=0;
i=0;
wseq=zeros(6,1);
rseq=zeros(6,1);
penseq=zeros(6,1);
tauseq=zeros(6,1);



@ ---------------------------------------------

iteration of policy function, wealth distribution,..

in the old and new steady state

---------------------------------------------- @

its=0;
do until its==2;
    its=its+1;
    if its==1;
        rep=rep0;
    else;
        rep=rep1;
    endif;  

    tau=rep/(2+rep);    /* income tax rate */

    tt=1;
    x0=kbar|nbar;
    {fvp,xf}=FixVMN(x0,&dgetss);

    kbar=xf[1];
    nbar=xf[2];

    "kbar~nbar: " xf';
    if its==1;
        aoptold=aopt;
        noptold=nopt;
        kbarold=kbar;
        nbarold=nbar;
        tauold=tau;
        wold=w;
        rold=r;
        penold=pen;
        copt=zeros(6,1);
        copt[1]=(1-tau)*w*nopt[1]-aopt[1];
        copt[2]=(1-tau)*w*nopt[2]+(1+r)*aopt[1]-aopt[2];
        copt[3]=(1-tau)*w*nopt[3]+(1+r)*aopt[2]-aopt[3];
        copt[4]=(1-tau)*w*nopt[4]+(1+r)*aopt[3]-aopt[4];
        copt[5]=pen+(1+r)*aopt[4]-aopt[5];
        copt[6]=pen+(1+r)*aopt[5];
        coptold=copt;
    else;
        aoptnew=aopt;
        noptnew=nopt;
        kbarnew=kbar;
        nbarnew=nbar;
        taunew=tau;
        wnew=w;
        rnew=r;
        pennew=pen;
        copt[1]=(1-tau)*w*nopt[1]-aopt[1];
        copt[2]=(1-tau)*w*nopt[2]+(1+r)*aopt[1]-aopt[2];
        copt[3]=(1-tau)*w*nopt[3]+(1+r)*aopt[2]-aopt[3];
        copt[4]=(1-tau)*w*nopt[4]+(1+r)*aopt[3]-aopt[4];
        copt[5]=pen+(1+r)*aopt[4]-aopt[5];
        copt[6]=pen+(1+r)*aopt[5];
        coptnew=copt;
    endif;

    if its==2 and _update==4;   /* computation of the Initial Guess for Broyden matrix */
        broyinit=CDJac(&dgetss,kbar|nbar,2);
        broy=zeros(2*nt,2*nt);
        i=0;
        do until i==nt;
            i=i+1;
            broy[i,i]=broyinit[1,1];
            broy[i+nt,i+nt]=broyinit[2,2];
        endo;
        broyinv=inv(broy);
    endif;
    "its: " its;
    "rep: " rep;
    "kbar: " kbar;

endo;   /* its */

save aoptold,aoptnew,noptold,noptnew,coptold,coptnew;

"old ~ new steady state : ";
"capital stock: " kbarold~kbarnew;
"employment: " nbarold~nbarnew;
"tau: " tauold~taunew; 

bsec=hsec;

@ ------------------------------------------------------

computation of the transition

------------------------------------------------------ @

/* initial guesses for K_t,k^s_t, N_t,n^s_t */
kst=zeros(nt,6);
ktold=zeros(nt,1);
ktnew=zeros(nt,1);
nst=zeros(nt,6);
ntold=zeros(nt,1);
ntnew=zeros(nt,1);
ktold=seqa(kbarold,(kbarnew-kbarold)/(nt-1),nt);
ntold=seqa(nbarold,(nbarnew-nbarold)/(nt-1),nt);
xold=ktold|ntold;
rept=ones(nt,1)*rep1;
taut=ones(nt,1)*taunew;

ktold1=ktold;
q=0;
kritt=1+tolt;
do until q==nqt or (q>1 and kritt<tolk);
    q=q+1;
    "iteration over sequence capital stock/employment: " q;
    wt=(1-alp).*ktold^alp.*ntold^(-alp);
    rt=alp.*ktold^(alp-1).*ntold^(1-alp)-del;
    pent=rept.*(1-taut).*wt.*ntold*3/2;
    tt=nt;
    if _update==1;
        {ktnew,ntnew}=getkn();
        kritt=meanc(abs(ktnew-ktold));
        ktold=phi*ktold+(1-phi)*ktnew;
        ntold=phi*ntold+(1-phi)*ntnew;
    elseif _update==2;
        /* computation of the Jacobian */
        yold=f1(xold);
        dx=CDJac(&f1,xold,2*nt);
/*
    "computational time: ";
    ?etstr(hsec-bsec);  
    "jacobian: "; wait;
        dx; wait;
*/
        xnew=xold-inv(dx)*yold;
/*
        "xold~xnew: ";
        wait;
        xold~xnew; wait;
*/
        kritt=meanc(abs(xnew[1:nt]-xold[1:nt]));
        ktnew=xnew[1:nt]; ntnew=xnew[nt+1:2*nt];
        xold=xnew;

    elseif _update==3 or _update==4;
        if q==1 and _update==3;
                broy=CDJac(&f1,xold,2*nt);
                broyinv=inv(broy);
        endif;
        yold=f1(xold);
        dx=-broyinv*yold;
        xnew=xold+dx;
        kritt=meanc(abs(xnew[1:nt]-xold[1:nt]));
        ynew=f1(xnew);
        /* update of Broyden */
        s1=xnew-xold;
        s2=(ynew-yold)-broy*s1;
        broy=broy+s2*s1'/(s1'*s1);
        broyinv=inv(broy);
        xold=xnew;
        ktnew=xnew[1:nt]; ntnew=xnew[nt+1:2*nt];
    endif;        

    "iteration: " q;
    "crit K_t: " kritt;
endo;


"computational time: ";
    ?etstr(hsec-bsec);  
"convergence after iteration: " q;
wait;


GraphSettings;
@
"Figure 9.7: " wait;
kt1=(ones(4,1)*kbarold|ktnew);
title("");
ylabel("Aggregate capital stock");
xlabel("Transition period");
xy(seqa(-3,1,nt+4),kt1);
wait;


"Figure 9.3: " wait;
_plegctl=1;
_pltype=1|6;
_plegstr="Old steady state\000New steady state";
xtics(1,6,1,0);
xlabel("Generation");
title("");
ylabel("Individual capital stock");
xy(seqa(1,1,6),aoptold~aoptnew); 
wait;


"Figure 9.4: " wait;
xtics(1,4,1,0);
_plegctl={1 5 1.5 0.3};
ylabel("Labor Supply");
xy(seqa(1,1,4),noptold~noptnew);
wait;
@

load aopt2, nopt2;
"Figure 9.5: " wait;

begwind;
window(1,2,0);
xlabel("Age");
xtics(1,6,1,0);
_pltype=1|6;
_plegctl={1 7 3 0.01};
_plegstr="Old steady state\000Individual born in t=-2";
_ptitlht=0.3;
xlabel("Age");
ylabel("");
title("Individual capital stock");
xy(seqa(1,1,6),aoptold~aopt2); 
nextwind;
xtics(1,4,1,0);
_plegctl={1 7 1.5 0.29};
title("Individual labor supply");
xy(seqa(1,1,4),noptold~nopt2); 
endwind;

wait;
save nt,ktold,nbarold,ntold;

@ ------------------------------------------------------

Procedures

------------------------------------------------------  @

  proc dgetss(x);
    local y0;
    y0=getss(x);
    retp(y0-x);
  endp;

  proc getss(x);
    local kbar,nbar,cy,n0,k2,k3,k4,k5,k6,x0,fvp,xf,knew,nnew;
    kbar=x[1];
    nbar=x[2];
    w=(1-alp)*kbar^alp*nbar^(-alp);
    r=alp*kbar^(alp-1)*nbar^(1-alp)-del;
    pen=rep*(1-tau)*w*nbar*3/2;
    wseq=ones(6,1)*w;
    rseq=ones(6,1)*r;
    penseq=ones(6,1)*pen;
    tauseq=ones(6,1)*tau;
/* initial guess for (k^1,k^2,..,k^6,n^1,..,n^4: n==0.3, consumption=80% of income */
    cy=0.8;
    n0=0.3;
    k2=(1-cy)*(1-tau)*w*n0;
    k3=(1-cy)*((1-tau)*w*n0+r*k2)+k2;
    k4=(1-cy)*((1-tau)*w*n0+r*k3)+k3;
    k5=(1-cy)*((1-tau)*w*n0+r*k4)+k4;
    k6=k5/2;    
    x0=k2|k3|k4|k5|k6|ones(4,1)*n0;
    {fvp,xf}=FixVMN(x0,&rftr);
    aopt=zeros(6,1);
    aopt[2:6]=xf[1:5];
    nopt=xf[6:9];
    /* computation of the aggregate capital stock and employment nbar */
    knew=meanc(aopt);
    nnew=meanc(nopt)*2/3;
    retp(knew|nnew);
  endp;


/* first-order condtions of the household for wages wseq, interest rates rseq, pensions penseq */
proc rftr(x);
    local y,k2,k3,k4,k5,k6,n1,n2,n3,n4,c1,c2,c3,c4,c5,c6;
    y=x;
    k2=x[1];
    k3=x[2];
    k4=x[3];
    k5=x[4];
    k6=x[5];
    n1=x[6];
    n2=x[7];
    n3=x[8];
    n4=x[9];
    c1=(1-tauseq[1])*wseq[1]*n1-k2;
    c2=(1-tauseq[2])*wseq[2]*n2+(1+rseq[2])*k2-k3;
    c3=(1-tauseq[3])*wseq[3]*n3+(1+rseq[3])*k3-k4;
    c4=(1-tauseq[4])*wseq[4]*n4+(1+rseq[4])*k4-k5;
    c5=penseq[5]+(1+rseq[5])*k5-k6;
    c6=penseq[6]+(1+rseq[6])*k6;
    y[1]=gam*(c1+psi)-(1-tauseq[1])*wseq[1]*(1-n1);
    y[2]=gam*(c2+psi)-(1-tauseq[2])*wseq[2]*(1-n2);
    y[3]=gam*(c3+psi)-(1-tauseq[3])*wseq[3]*(1-n3);
    y[4]=gam*(c4+psi)-(1-tauseq[4])*wseq[4]*(1-n4);
/* in period 1: unexpected change, therefore, prior to period 1, all agents
    behave as in the old steady state */
    if tt<=0;
        c1=coptold[1];
        y[1]=n1-noptold[1];
        y[5]=k2-aoptold[2];
    else;
        y[5]=1/beta-uc(c2,1-n2)/uc(c1,1-n1)*(1+rseq[2]);
    endif;
    if tt<=-1;
        c2=coptold[2];
        y[2]=n2-noptold[2];
        y[6]=k3-aoptold[3];
    else;
        y[6]=1/beta-uc(c3,1-n3)/uc(c2,1-n2)*(1+rseq[3]);
    endif;
    if tt<=-2;
        c3=coptold[3];
        y[3]=n3-noptold[3];
        y[7]=k4-aoptold[4];
    else;
        y[7]=1/beta-uc(c4,1-n4)/uc(c3,1-n3)*(1+rseq[4]);
    endif;
    if tt<=-3;
        c4=coptold[4];
        y[4]=n4-noptold[4];
        y[8]=k5-aoptold[5];
    else;
        y[8]=1/beta-uc(c5,1)/uc(c4,1-n4)*(1+rseq[5]);
    endif;
    if tt<=-4;
        c5=coptold[5];
        y[9]=k6-aoptold[6];
    else;
        y[9]=1/beta-uc(c6,1)/uc(c5,1)*(1+rseq[6]);
    endif;
    retp(y);
endp;




proc f1(x);
    local y1,y2,y0;
    ktold=x[1:nt];
    ntold=x[nt+1:2*nt];
    wt=(1-alp).*ktold^alp.*ntold^(-alp);
    rt=alp.*ktold^(alp-1).*ntold^(1-alp)-del;
    pent=rept.*(1-taut).*wt.*ntold*3/2;
    tt=nt;
    {y1,y2}=getkn();
    y0=y1|y2;
    retp(y0-x);
endp;


proc(2)=getkn();
    local ktnew,ntnew,x0,fvp,xf,i,tp,aopt2,nopt2;
    tt=nt+1; 
    /* computation of the enty n^1_nt where k^1_nt=0 and k^2_nt=\bar k^2*/
    ktnew=zeros(nt,1);
    ntnew=zeros(nt,1);

    /* computation of the optimal allocation of an agent born in period tt */
    do until tt==-4;
        tt=tt-1;
        /* wage, interest rate and pensions over the lifetime of the individual */
        if tt>nt-5; /* agent is alive after period nt */
            wseq[1:nt-tt+1]=wt[tt:nt];
            wseq[nt-tt+2:6]=ones(6-nt+tt-1,1)*wnew;
            rseq[1:nt-tt+1]=rt[tt:nt];
            rseq[nt-tt+2:6]=ones(6-nt+tt-1,1)*rnew;
            penseq[1:nt-tt+1]=pent[tt:nt];
            penseq[nt-tt+2:6]=ones(6-nt+tt-1,1)*pennew;
            tauseq[1:nt-tt+1]=taut[tt:nt];
            tauseq[nt-tt+2:6]=ones(6-nt+tt-1,1)*taunew;
        elseif tt<1;
            wseq[1:1-tt]=ones(-tt+1,1)*wold;
            wseq[2-tt:6]=wt[1:4+tt+1];
            rseq[1:1-tt]=ones(-tt+1,1)*rold;
            rseq[2-tt:6]=rt[1:4+tt+1];
            penseq[1:1-tt]=ones(-tt+1,1)*penold;
            penseq[2-tt:6]=pent[1:4+tt+1];
            tauseq[1:1-tt]=ones(-tt+1,1)*tauold;
            tauseq[2-tt:6]=taut[1:4+tt+1];
        else;
            wseq[1:6]=wt[tt:tt+5];
            rseq[1:6]=rt[tt:tt+5];
            penseq[1:6]=pent[tt:tt+5];
            tauseq[1:6]=taut[tt:tt+5];
        endif;

        /* computation of allocation of the individual born in period tt */
        x0=aopt[2:6]|nopt[1:4];
        {fvp,xf}=FixVMN(x0,&rftr);
        aopt[1]=0;
        aopt[2:6]=xf[1:5];
        nopt[1:4]=xf[6:9];
        tp=tt-1; i=0;
        do until tp==nt or i==6;
            tp=tp+1; i=i+1;
            if tp>0;
                ktnew[tp]=ktnew[tp]+1/6*aopt[i];
                if i<=4;
                    ntnew[tp]=ntnew[tp]+1/6*nopt[i];
                endif;
            endif;
        endo;

        if tt==-2;
            save aopt2=aopt, nopt2=nopt;
        endif;
    endo;

    retp(ktnew,ntnew);
endp;



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

@ ----------------------------- Ch8_Toolbox.src -----------------------------

  This file contains Gauss procedures introduced in Chapter 8.

  ChebEval: Evaluate a Chebyshev polynomial

  ChebCoef: Obtain Chebyshev approxiamtion from a regression

  QuasiNewton: Function minimization using a quasi Newton method with BFGS update
               and line search

  QNStep     : Line search for QuasiNewton

  GradTest   : computes relative gradient (used to test whether the algorithm is near a minimizer)
  
  MinStep    : computes minimal step size

  ParTest    : computes relative change of parameter (used to test whether the algorithm converged)

  GSearch    : Function minimization using a genetic search algorithm  

  CDJac      : Central difference Jacobian
  
  FixvMN1    : Non-Linear Equations Solver
 
  LSolve     : Used by FixvMN1

  NRStep     : Line search for Newton-Raphson

----------------------------------------------------------------------------- @


@ ----------------------------------------- QuasiNewton ----------------------------------------------------

    Purpose: Find the minimizer of a user supplied function

    Usage: {x1,crit}=QuasiNewton(x0,&f)

    Input:  x0 : n times 1 vector, initial value
            &f : pointer to f, where f must return a scalar

    Output: x1 : n times 1 vector, the found solution
           crit: 5 times 1 vector, crit[1]:=return code: 0: normal termination
                                                         1: function evaluation at x not possible
                                                         2: step2<Mstep, and no further reduction in f possible,
                                   crit[2]:=the maximum absolute value of the scaled gradient at x1
                                   crit[3]:=the maximum relative change of x between the last two iterations
                                   crit[4]:=the value of f at x1
                                   crit[5]:=number of iterations

    Gobals:     _QN_Print = 1 print to screen, =0 do not print messages to screen, default
                _GradTol  = 1 use defaul value of MachEps^(1/3)
                _ParTol   = 1.e-7 (about 7 good digits in x), default

---------------------------------------------------------------------------------------------------------- @

declare matrix _QN_Print    != 0;      @ Do not print messages do screen          @          
declare matrix _QN_GradTol  != 1;      @ use default value for GradTol            @
declare matrix _QN_ParToL   != 1.e-7;  @ about 7 good digits in computation of x1 @

proc(2)=QuasiNewton(x0,&f);

    local df1,df2,x1,x2,f1,f2,crit, mStep, step1, step2,dx,dgrad,h,itn,typf,maxit,rc,ptol,gtol;

    local f:proc;
      

    /* initialize */
    crit=zeros(5,1);

    if _QN_GradTol; @ Gradient tolerance @
          gTol=MachEps^(1/3);
    else;
          gTol=_QN_GradTol;
    endif;

    ptol =MachEps^(2/3); @ Parameter tolerance used in computation of stepsize @
    h=eye(rows(x0));
    x1=x0;
    maxit=500;
    
    f1=f(x1);
    @df1=gradp(&f,x1);@
    df1=CDJac(&f,x1,1);
    typf=1;
    crit[2]=GradTest(df1,x1,f1,typf);

    /* Check whether initial value is already a good solution canditate */
    if crit[2]<1.e-3*gTol; 
             crit[1]=0; crit[5]=0; crit[3]=miss(1,1); crit[4]=f1;
             retp(x1,crit);
    endif;    

    /* iterate */
          itn=1;    
      crit[3]=1;

    do while itn<maxit;
       if _QN_Print; 
          locate 1,1;?"QN-Iteration #= " ftos(itn,"*.*lf",5,0);
          locate 2,1;?"gTol   := " ftos(crit[2],"*.*lf",15,8);
          locate 3,1;?"pTol   := " ftos(crit[3],"*.*lf",15,8);
          locate 4,1;?"f(x)   := " ftos(crit[4],"*.*lf",15,8);
       endif;

       dx=solpd(-df1',h);

       /* Reduce step size,if f cannot be computed at full Newton step */
    @   mStep=MinStep(x1,dx,ones(rows(x1),1),pTol);@
        mStep=pTol;
       step1=1;
@       do while (scalmiss(f(x1+step1*dx)) and (step1>mStep)); @
        do while (scalmiss(f(x1+step1*dx)) and (step1>pTol));
          step1=step1/2;
       endo;
       if (step1<mStep);
          crit[1]=1;
          retp(x1,crit);
       endif;
       dx=step1*dx;       
       {step2,rc}=QNStep(x1,dx,f1,df1,&f);
       if _QN_Print;
          locate 5,1; ?"MStep=   " ftos(MStep,"*.*lf",20,15);
          locate 6,1; ?"Step1=   " ftos(step1,"*.*lf",20,15);
          locate 7,1; ?"Step2=   " ftos(step2,"*.*lf",20,15);
       endif;
       dx=step2*dx;
       x2=x1+dx;
       f2=f(x2);
       crit[4]=f2;
       @df2=gradp(&f,x2);@
        df2=CDJac(&f,x2,1);
       typf=(typf+f2)/(itn+1);

       /* Check for convergence */       
       crit[2]=GradTest(df2,x2,f1,1.0);
       crit[3]=ParTest(x2,dx,ones(rows(x2),1));
       if ((crit[2] > gtol) and rc);
          crit[1]=2;
          retp(x2,rc);
       endif;
       if crit[2]<gtol;
          if rc; crit[1]=0; retp(x2,crit); endif;          
          if crit[3]<_QN_ParTol;
             crit[1]=0;
             retp(x2,crit);
          endif;
       endif;

       /* Update h */
      dgrad=df2-df1;
      h=h-(((h*dx)*(dx'h))/(dx'h*dx))+((dgrad'dgrad)/(dgrad*dx));
      df1=df2;
      x1=x2;
      f1=f2;
      itn=itn+1;
      crit[5]=itn;
 endo;
 retp(x1,2);

endp;

@ --------------------------------------- GradTest -------------------------------------------------------

    Purpose: Computes the relative gradient (the elasticity of f at x) and returns the
             maximum absolute value. This is one of the stopping criteria suggested
             by Dennis and Schnabel (1983).

    Usage: crit=GradTest(df,x,typf);

    Input: df := n times 1 vector, the gradient of f at x
            x := n times 1 vector, the point at which df is evaluated
           fx := scalar, the value of f at x
          typf:= scalar, the typical value of f at or near x

    Output: crit := scalar the relative gradient
----------------------------------------------------------------------------------------------------- @

proc(1)=GradTest(df,x,fx,typf);

    local crit, i, n;
    n=rows(x);
    crit=zeros(n,1);
    i=1;
    do until i>n;
       crit[i]=abs(df[i])*maxc(abs(x[i])|1.0);
       crit[i]=crit[i]/maxc(abs(fx)|abs(typf));
             i=i+1;
    endo;

  retp(maxc(crit));

endp;


@--------------------------------------- MinStep ----------------------------------------------------

  Purpose: Compute the minimal step size so that x1=x0+mstep*dx=0
           (in terms of the parameter tolerance criterium pTol)

  Usage: mstep=MinStep(x,dx,typx,pTol);

  Input:  x := n times 1 vector
         dx := n times 1 vector, change in x
        typx:= n times 1 vector, typical elements of x          
        pTol:= scalar


  Output: mstep:= scalar

---------------------------------------------------------------------------------------------------   @

proc(1)=MinStep(x,dx,typx,pTol);

   local temp, converge, i, n;

   n=rows(x);
   i=1;
   converge=0;

   do until i>n;

       temp=abs(dx[i])/maxc(abs(x[i])|abs(typx[i]));
       if (temp > converge);
            converge=temp;
       endif;
       i=i+1;

   endo;

   retp(PTol/converge);

endp;


@ ----------------------------------- ParTest -------------------------------------------------- 

   Purpose: Compute relative change in x

   Usage:   crit=ParTest(x,dx,typx);

   Input:   x := n times 1 vector
           dx := n times 1 vector, change in x
         typx := n times 1 vector, typical elements of x

   Output: crit := scalar

----------------------------------------------------------------------------------------------- @

 proc(1)=ParTest(x,dx,typx);

   local i, n, crit;

    n=rows(x);
    crit=zeros(n,1);
    i=1;
    do until i>n;
       crit[i]=abs(dx[i])/maxc(abs(x[i])|abs(typx[i]));
             i=i+1;
    endo;

   retp(maxc(crit));

endp;



@ --------------------------------------------- QNStep ----------------------------------------

 Purpose:  Find the step size s so that the Quasi-Newton algorithm
           always moves in the direction of a (local) minimum of
           f(x)


 Usage:    {s,rc}=GetStep4(x0,dx0,f0,df,&f);

 Input:    x0 := n times 1 vector, the initial point
           dx0:= n times 1 vector, the Newton direction
           f0 := f(x0)
           df := 1 times n vector, the gradient of f at x_0
           &f := pointer to the function f

 Output:   s  := admissible stepsize
          rc  := return code: rc=0=normal exit, rc=1=minimal stepsize reached
------------------------------------------------------------------------------------------------ @

      
proc(2)=QNStep(x0,dx0,f0,df,&f);

  local s, x1, s1, s2, smult, smin, smax, pTol, f1, f2, amat, bvec, ab, dfdx, nobs, rc, disc;
  local f:proc;


  /* Fixed parameters of the algorithm */
  smult=1.0e-4; smin=0.1; smax=0.5; pTol=MachEps^(2/3);
  
  /* Initialize */
  s=1.0;
  s1=1.0;
  s2=1.0;
  rc=0;
  amat=zeros(2,2); bvec=zeros(2,1);   
 dfdx=df*dx0;                        @ df(x0)*dx0                          @
 f1=f(x0+dx0);
 f2=f1;

 /* Try the full Newton step s=1 */
 if f1 <= f0 + smult*dfdx;
    retp(s1,rc);
 else;
    s=-dfdx/(2*(f1-f0-dfdx));
    if s<smin; s=smin; endif;
    if s>smax; s=smax; endif;
    x1=x0+s*dx0;
    f2=f(x1);
 endif;
 s2=s;

 /* Reduce s2 further unless f2 < f0 + s2*smult*dfdx */
 do while (f2 > (f0 + smult*s2*dfdx) );
    amat[1,1]=1/(s2^2);  amat[1,2]=-1/(s1^2); amat[2,1]=-s1/(s2^2); amat[2,2]=s2/(s1^2);
    bvec[1]=f2-s2*dfdx-f0; bvec[2]=f1-s1*dfdx-f0;
    ab=(amat*bvec)/(s2-s1);
    if (ab[1,1] == 0.0);
            s=-dfdx/(2*ab[2,1]);
      else;
          disc=(ab[2,1]^2)-3*ab[1,1]*dfdx;  
          if (disc < 0.0);
             s=s2*smax;
           elseif (ab[2,1] <= 0.0);
             s=(-ab[2,1]+sqrt(disc))/(3*ab[1,1]);
           else;
             s=-dfdx/(ab[2,1]+sqrt(disc));
           endif;
      endif;           
     if s < s2*smin; s=s2*smin; endif;
     if s > s2*smax; s=s2*smax; endif;
    @ if s<MinStep(x0,s*dx0,ones(rows(x0),1),pTol); retp(s,1); endif; @
      if s<pTol; retp(s,1); endif;
     s1=s2;
     s2=s;
     f1=f2;
     x1=x0+s2*dx0;
     f2=f(x1);

 endo;

 retp(s2,rc);

endp;
      
@ ----------------------------------------------- GSearch1 ----------------------------------------------------

   Purpose: function minimization using a genetic search algorithm

   Usage:   {x1,f1}=GSearch(NPar,Npop,Ngen,&f)

   Input:   Npar :=  integer, the number of elements in x
            Npop :=  integer, the size of the population
            Ngen :=  integer, the number of iterations (number of generations)
               &f:=  pointer to the function that evalautes the fitness of a canditate solution
                     (the function whose minimim is to be found)

   Output:    x1 := Npar times 1 vector, the best canditate solution
              f1 := scalar, f(x1)

------------------------------------------------------------------------------------------------------------ @
 

proc(2)=GSearch1(NPar,Npop,Ngen,&F);

  local Genes0, Genes1, Par, fit, fit1, fit2, i, i1, j1, j2, j3, j4, r1,
        p, p1, p2, c1, c2, lambda, m_ind, dpar, smin, mu1, mu2,
        fmin, fmax, fstdv, gnum, get_i, get_p, get_n, b, probco, probmut;

  local F:proc;

  /* Parameters that determine the behavior of GSearch */
  probco=0.95;     @ probability of crossover @

  mu1=0.15;        @ parameters of probability function for mutations to apply @
  mu2=0.33;
  b=2;             @ used in mutation operations @



  /* Initialize population */
  Genes0=zeros(npop,Npar+2);
  Genes1=Genes0;

  locate 2,2;
  ?"Initialize Population";
  i=1;
  do until i>npop;

     locate 2,25;
     ?ftos(i,"*.*lf",4,0);  
     Par=rndn(Npar,1);  
     fit=F(Par);     
     do while scalmiss(fit);     
       Par=rndn(Npar,1);
       Fit=F(Par);
     endo;
    Genes0[i,1:Npar]=Par';
    Genes0[i,NPar+1]=fit;
      i=i+1;
  endo;

  /* Print summary statistics */
   fmin=minc(Genes0[.,NPar+1]);
   fmax=maxc(Genes0[.,NPar+1]);
  fstdv=stdc(Genes0[.,NPar+1]);
  locate 2,2;  
  ?"Best Cromosome:     " ftos(fmin,"*.*lG",10,6);
  locate 3,2;
  ?"Worst Cromosome:    " ftos(fmax,"*.*lG",10,6);
  locate 4,2;
  ?"Standard Deviation: " ftos(fstdv,"*.*lG",10,6);
   
  /* Start Selection Process */

  gnum=1;

  do until gnum>Ngen;  @ Begin Loop over successive generations @
      Genes1=zeros(npop,NPar+2);
      locate 6,2;
      ?"Generation No.: " ftos(gnum,"*.*lf",3,0);
      i=1;
      i1=1;
      do until i>(npop/2); @ Begin Loop over selection from one generation @
      locate 6,25;?"Current population size=" ftos(i1,"*.*lf",3,0);
 
     /* Draw two pairs from Genes0 */   
      get_i=rndu(4,1);
      P=zeros(4,1);  @ holds indices of parents @
      P=1+round((npop-1)*get_i);
         
      /* Find the two fittest parents */
      Get_p=p~genes0[p,NPar+1];
      Get_p=sortc(Get_p,2);
      p1=get_p[1,1];
      p2=get_p[2,1];

      /* Decide whether crossover is to be performed */
      if rndu(1,1)<=probco;
         /* Find out which method is to be employed */
         m_ind=rndu(1,1);
         if (m_ind<1/3);  @ arithmetic crossover @
            lambda=rndu(1,1);
            C1=lambda*Genes0[p1,1:Npar]+(1-lambda)*Genes0[p2,1:Npar];
            C2=lambda*Genes0[p2,1:Npar]+(1-lambda)*Genes0[p1,1:Npar];
         endif;
         if ((m_ind>=1/3) and (m_ind<2/3)); @ Single-point crossover @
            lambda=1+round((Npar-2)*rndu(1,1));
            C1=Genes0[p1,1:lambda]~Genes0[p2,1+lambda:Npar];
            C2=Genes0[p2,1:lambda]~Genes0[p1,1+lambda:Npar];
         endif;
         if (m_ind>=2/3); @ Shuffle crossover @
            j3=1;
            C1=zeros(1,Npar);
            C2=C1;
            do until j3>Npar;
               lambda=rndn(1,1);
               if lambda >= 0.5; 
                  C1[j3]=Genes0[p2,j3];
                  C2[j3]=Genes0[p1,j3];
               else;
                  C1[j3]=Genes0[p1,j3];
                  C2[j3]=Genes0[p2,j3];
               endif;
               j3=j3+1;
             endo;
          endif;
      else; @ otherwise Children and Parents are identically @
         C1=Genes0[p1,1:Npar];
         C2=Genes0[p2,1:Npar];
      endif;

      /*  Mutation operations: on C1 */
      j4=1;
      do until j4>3;
         gosub compute_probmut;
         if rndu(1,1)<probmut;
           dpar=rndn(1,1);
           r1=rndu(1,1);
           if rndu(1,1)>0.5;
              C1[j4]=C1[j4]+dpar*(1-(r1^((1-gnum/Ngen)^b)));
           else;
              C1[j4]=C1[j4]-dpar*(1-(r1^((1-gnum/Ngen)^b)));
           endif;
         endif;
         j4=j4+1;
      endo;
      /*  Mutation operations: on C2 */
      j4=1;
      do until j4>3;
         gosub Compute_probmut;
         if rndu(1,1)<probmut;
           dpar=rndn(1,1);
           r1=rndu(1,1);
           if rndu(1,1)>0.5;
              C2[j4]=C1[j4]+dpar*(1-(r1^((1-gnum/Ngen)^b)));
           else;
              C2[j4]=C1[j4]-dpar*(1-(r1^((1-gnum/Ngen)^b)));
           endif;
         endif;
         j4=j4+1;
      endo;
      /* Evalutate Fittness of Children and select fittest pair
      ** for next generation */
       Fit1=F(C1'); 
       Fit2=F(C2');
       get_n=C1~Fit1|C2~Fit2|Genes0[p1,1:Npar+1]|Genes0[p2,1:Npar+1];
       get_n=packr(sortc(get_n,Npar+1));
       Genes1[i1:i1+1,1:Npar+1]=get_n[1:2,.];
        i=i+1;
        i1=i1+2;
      
      endo; @ end loop over selection from one generation @

     /* Apply Elitism */
    Genes1=sortc(Genes1,Npar+1);
    Genes0=sortc(Genes0,Npar+1);
        locate csrlin+1,2;
        ?"Best old chromosome: " Genes0[1,Npar+1];        
        locate csrlin,2;
        ?"Best new chromosome: " Genes1[1,Npar+1];
        
    if Genes1[1,Npar+1]>Genes0[1,Npar+1];
        Genes1[npop,1:Npar+1]=Genes0[1,1:Npar+1];        
    endif;

   /* Interchange Populations  */
   Genes0=Genes1;
   
   gnum=gnum+1;

endo;

smin=minindc(genes0[.,Npar+1]);

retp(genes0[smin,1:Npar]',genes0[smin,Npar+1]);

  compute_probmut:
   probmut=mu1 + mu2/gnum;
  return;

endp;



/* CDJac
**
** usage: Jac=CDJac(&f,x0,n)
**
**
** purpose: computes a central difference approximation of the Jacobian
**          matrix of a system of n non-linear functions y_i=f^i(x), where
**          x is a column vector of dimension m.
**
** input:  &f: pointer to the routine that returns the m-vector f(x)
**         x0: m vector x0, the point at which the derivatives are to be evaluated
**         n : the size of the vector f(x)
**
** output: Jac: n by m matrix of parital derivatives
**
** algorithm: based on (A.2.8) in Heer and Maussner, see also Dennis and Schnabel (1983), 
**            Algorithm A5.6.4.
*/

proc(1)=CDJac(&f,x0,n);

  local eps, i, j, m, h, df, f1, f2, x1, x2, temp;
  local f:proc;
  m=rows(x0);
  df=zeros(n,m);
  eps=MachEps^(1/3);
  x1=x0;
  x2=x0;
  i=1;

  do until i>m;
     if x0[i]<0;
        h=-eps*maxc(abs(x0[i])|1.0);
     else;
        h=eps*maxc(abs(x0[i])|1.0);
    endif;
    temp=x0[i];
    x1[i]=temp+h;
    x2[i]=temp-h;
        h=x1[i]-temp; @ Trick to increase precision slightly, see Dennis and Schnabel (1983), p. 99 @
    f1=f(x1);
    f2=f(x2);
    j=1;
    do until j>n;
       df[j,i]=(f1[j]-f2[j])/(2*h);
       j=j+1;
    endo;
    x1[i]=x0[i];
    x2[i]=x0[i];
    i=i+1;
 endo;
 retp(df);

endp;

/* FixvMN1: Solves a system of non-linear equations using a modified Newton Method
**
** Usage:   {fx,x}=FixVMN(x0,&F)
**
**  Input:  &F  := Pointer the vector valued function F(x), whose
**                 zero, F(x1)=0, is to be computed. 
**
**          x0  := k times 1 vector of starting values
** 
**
** 
**  Output: fx  := k times 1 vector, the value of f at the solution x1
**
**          x1  := k times 1 vector, the approximate solution to F(x1)=0
**
**
** Gobals: _MNR_Print=1 (0) do (not) print messages to the screen
**
**         _MNR_Gobal=1 (0) do (not) use line search
**
**         _MNR_QR=1 (0) do (not) use QR factorization to solve for the Newton step
*/

declare matrix _MNR_Global != 0;  @ Do not use GetStep in FixvMN             @
declare matrix _MNR_Print  != 0;  @ FixvMN print no messages to the screen   @
declare matrix _MNR_QR     != 0;  @ do not use QR factorization for solution @


proc(2)=FixVMN(x0,&f);

  local x1, x2, dx, df, dg, crit, crit1, crit2, stopc, itn, step1, step2, maxit;
  local f:proc;

  /* Initialize */
    maxit=5000;     @ stop after 5000 Iterations @
    stopc=1e-8;     @ stopping criterium @
       x1=x0;
      itn=1;     

  /* Start Iterations */
  itn=1;
  crit=1;crit1=1;crit2=1;
  if _MNR_Print; cls; endif; @ clear screen if output is printed to the screen @

  do until ((crit<stopc) or (itn > maxit));      @ start iterations @
      if _MNR_Print;
         locate 1,2;
         ?"Step No: " ftos(itn,"*.*lf",5,0) "Convergence criterion:" crit;
         if _MNR_Global;
           locate 2,2;
           ?"Minimization criterion: "  crit1 " crit1<crit2: " (crit1<crit2);
         endif;
      endif;
        df=gradp(&F,x1);
      if _MNR_Global;  dg=F(x1)'df; endif;

     if _MNR_QR;
        dx=LSolve(df,-f(x1)); @ use the QR-factorization @
     else;
        dx=-inv(df)*F(x1);  
     endif;    
         x2=x1+dx;
      step1=1;
      do while scalmiss(f(x1+step1*dx));
         step1=0.75*step1;
      endo;
      if _MNR_Print;
       locate 3,2;
       ?"Step1= " step1;
      endif;
      if step1<1.0e-16;
           ?"Not able to find x within admissible bounds";
           ?"Results may be missleading"; 
           ?"Press any key to continue.";
           wait;
           retp(F(x1),x1);
      endif;
      dx=step1*dx;
      if _MNR_Global;  step2=NRStep(x1,dx,dg,&f); else; step2=1; endif;

      if _MNR_Print;
         locate 4,2;
         ?"Step2= " step2;
      endif;

      if step2<0.0;
            ?"Stepsize cannot be further reduced."; 
            ?"Algorythm may stuck at a local minimum";
            ?"Results may be misleading";
            ?"Press any key to continue";
            wait;
            retp(f(x1),x1);
       endif;               
       x2=x1+step2*dx;
       crit=maxc(abs(F(x2)));
       if _MNR_Global;     crit2=crit1;   crit1=(F(x2)'F(x2))/2;    endif;
       x1=x2;
      itn=itn+1;
   endo;
   if itn >= maxit;
        ?"Maximum number of iterations exceeded! Returned solution may be missleading";
        ?"Press any key to continue";wait;
   endif;
retp(f(x1), x1);

endp;



/* LSolve:
**
**  Purpose: Solve a system of linear equations Ax=b using the QR-decomposition
**
**  Usage:  x=LSolve(A,b)
**
**  Input:  A: square matrix of dimension n
**          b: vector of dimension n
**
**  Output: x: vector of dimension n
*/   

proc(1)=LSolve(a,b);

 local r, d, q;
 {q,r}=qqr(a); 
retp(qrsol(q'b,(q'q)*r));

endp;

/* NRStep:
**
** Purpose:  Find the step size s so that the Newton-Raphson algorithm
**           always moves in the direction of a (local) minimum of (1/2)(f(x)'f(x))
**
**
** Usage:    s=NRStep(x0,dx0,dg,&f);
**
** Input:    x0 := n times 1 vector, the initial point
**           dx0:= n times 1 vector, the Newton direction
**           dg := 1 times n vector, the gradient of (1/2)(f'f) at x_0
**           &f := pointer to the function whose zero x solves the system of equations
**
** Output:   s  := admissible stepsize
*/

      
proc(1)=NRStep(x0,dx0,dg,&f);

  local s, x1, s1, s2, smult, smin, smax, stol, tol, g0, g1, g2, amat, bvec, ab, dgdx, disc;
  local f:proc;


  /* Fixed parameters of the algorithm */
  smult=1.0e-4; smin=0.1; smax=0.5; stol=1.e-11;
  
  /* Initialize */
  s1=1.0;
  amat=zeros(2,2); bvec=zeros(2,1);
   g0=(1/2)*(f(x0)'f(x0));
 dgdx=dg*dx0;                        @ dg(x0)*dx0                          @
 g1=(1/2)*(f(x0+dx0)'f(x0+dx0));


 /* Try the full Newton step s=1 */
 if g1 <= g0 + smult*dgdx;
    retp(s1);
 else;
    s=-dgdx/(2*(g1-g0-dgdx));
    if s<smin; s=smin; endif;
    if s>smax; s=smax; endif;
    x1=x0+s*dx0;
    g2=(1/2)*(f(x1)'f(x1));
 endif;
 s2=s;

 /* Reduce s2 further unless g2 < g0 + s2*smult*dgdx */
 do while (g2 > (g0 + smult*s2*dgdx) );

    amat[1,1]=1/(s2^2);  amat[1,2]=-1/(s1^2); amat[2,1]=-s1/(s2^2); amat[2,2]=s2/(s1^2);
    bvec[1]=g2-s2*dgdx-g0; bvec[2]=g1-s1*dgdx-g0;
    ab=(amat*bvec)/(s2-s1);

    if (ab[1,1] == 0.0);
            s=-dgdx/(2*ab[2,1]);
      else;
          disc=(ab[2,1]^2)-3*ab[1,1]*dgdx;  
          if (disc < 0.0);
             s=s2*smax;
           elseif (ab[2,1] <= 0.0);
             s=(-ab[2,1]+sqrt(disc))/(3*ab[1,1]);
           else;
             s=-dgdx/(ab[2,1]+sqrt(disc));
           endif;
      endif;               

    if s < s2*smin; s=s2*smin; endif;
    if s > s2*smax; s=s2*smax; endif;

    tol=sqrt((s*dx0)'(s*dx0))/(1+sqrt(x0'x0));
    if tol < stol; retp(-1.0); endif;
    s1=s2;
    s2=s;
    g1=g2;
    x1=x0+s2*dx0;
    g2=(1/2)*(f(x1)'f(x1));

 endo;

retp(s2);

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

