@ ------------------------------- Benchmark_LL.g ------------------------------------

    31 December 2007
    11 October 2011, last revision
    Alfred Maußner

    Compute impulse responses and second moments
    from the Benchmark real business cycle model of Heer and Maußner, 2nd. ed.,
    using log-linear approximations of the model's policy functions.

    External procedures:
    
    SolveLA from SolveLA.src
    Graphsettings from Tools.src
    HPFilter from Filter.src
    
    
    Note: If you are using the program without having included the required
          procedures into a Gauss library (e.g., the user.lcg) the program will
          terminate with error messages. To cirumvent this problem you can
          uncomment the lines 27-29.
          
------------------------------------------------------------------------------- @

new; @ clear memory @

library pgraph;

@ #include SolveLa.src; @
@ #include Tools.src;   @
@ #include Filter.src;  @

/* The following code definces a Gauss structure 
** to provide variable names and other information required .
*/

struct Variable 
{string varname;
 string vartype;    @ u=control variable, x=state variable, l=costate variable @
 scalar varpos;     @ the index of the variable in ut, xt, or lamt             @
 scalar crosscorr;  @ =1 the program compute cross correlations with this variable @
 scalar relsx;      @ =1 the program computes standard deviations relative to the standard deviation of this variable @
 scalar varprint;   @ =1 print second moments for this variable @
 scalar varplot;    @ =1 plott impulse responses of this variable @
 scalar plotno;     @ the number of the graphics panel in which the variable is plotted, up to four different panels are supported @
};

/* The next part must be modified by the user:
**
** _nu=number of controls
** _nx=number of states
** _nl=number of costates
** _ny=number of other variables of interest to the researcher
** _nz=number of shocks
**
** nobs=length of simulated time series
** nobs1=length of impulse response function
** nofs=number of simulations from whicht to compute averages of second moments
**
** HPL=0 if simulated data are to be passed through the Hodrick-Prescott filter
**       you must set HPL equal to the value of the smoothing parameter lambda
**
** _IR=1 if impulse response functions are to be computed
** _Moments=1 if time series moments are to be computed 
**
** _Scale=1 if all panels in graphic output share the same scale (=0 each panel will be scaled individually)
**
** outfile=string, the name of the output file
*/ 
_nu=6;
_nx=1;
_nl=1;
_nz=1;
_ny=0;

nobs=80;
nobs1=30;
nofs=300;
HPL=1600;

_Moments=0;
_IR=1;
_scale=1;
_nocolor=1;

outfile="Benchmark_LL.txt";

/********************************************************************************/

_nvar=_nu+_nx+_nl+_ny+_nz;

struct Variable Var;
Var=reshape(Var, _nvar, 1);

/* The next part must be modified by the user:
**
** You must put in the names of all variables in your model,
** and indicate whether the respective variable is a control,
** a state or a costate variable. Thus, there must be nvar statements of Var[i].
**
** Furthermore you must indicate for which variables you want to
** plott impulse responses, and for which variables you want
** to compute and output second moments.
**
** The order in which the variables appear below is also the
** order in which results will be ouptput.
**
** for example:
*/

Var[1].varname="Output";
Var[1].vartype="u";
Var[1].varpos=3;
Var[1].varprint=1;
Var[1].crosscorr=1;
Var[1].varplot=1;
Var[1].plotno=3;

Var[2].varname="Consumption";
Var[2].vartype="u";
Var[2].varpos=1;
Var[2].varprint=1;
Var[2].crosscorr=0;
Var[2].varplot=1;
Var[2].plotno=3;

Var[3].varname="Investment";
Var[3].vartype="u";
Var[3].varpos=4;
Var[3].varprint=1;
Var[3].crosscorr=0;
Var[3].varplot=1;
Var[3].plotno=3;

Var[4].varname="Hours";
Var[4].vartype="u";
Var[4].varpos=2;
Var[4].varprint=1;
Var[4].crosscorr=0;
Var[4].varplot=1;
Var[4].plotno=2;

Var[5].varname="Real Wage";
Var[5].vartype="u";
Var[5].varpos=5;
Var[5].varprint=1;
Var[5].crosscorr=0;
Var[5].varplot=1;
Var[5].plotno=2;

Var[6].varname="Real Rate of Interest";
Var[6].vartype="u";
Var[6].varpos=6;
Var[6].varprint=1;
Var[6].crosscorr=0;
Var[6].varplot=0;
Var[6].plotno=1;

Var[7].varname="Capital Stock";
Var[7].vartype="x";
Var[7].varpos=1;
Var[7].varprint=0;
Var[7].crosscorr=0;
Var[7].varplot=1;
Var[7].plotno=4;


Var[8].varname="Marginal Utility of Consumption";
Var[8].vartype="l";
Var[8].varpos=1;
Var[8].varprint=0;
Var[8].crosscorr=0;
Var[8].varplot=0;
Var[8].plotno=1;

Var[9].varname="Productivity";
Var[9].vartype="z";
Var[9].varpos=1;
Var[9].varplot=1;
Var[9].plotno=1;

/***************************************************************************************/

Cu=zeros(_nu,_nu);
Cxl=zeros(_nu,_nx+_nl);
Cz=zeros(_nu,_nz);
Dxl=zeros(_nx+_nl,_nx+_nl);
Fxl=Dxl;
Du=zeros(_nx+_nl,_nu);
Fu=Du;
Dz=zeros(_nx+_nl,_nz);
Fz=Dz;
_Rho=zeros(_nz,_nz);
_Sigma=zeros(_nz,_nz);

/* The next part must be modified by the user:
**
** First, you must insert the parameters of your model.
**
** Second, it is advisible to compute derived parameters
** and elasticities that appear in the matrices.
**
** Third, you must set up the matrices.
**
**
** for example:
*/

@ First Part @
a=1.005;
alpha=0.27;
beta=0.994;
eta=2.0;
nstar=0.13;
delta=0.011;
rhoZ=0.90;
sigmaz=0.0072;

@ Second Part @

yk=(a^eta-beta*(1-delta))/(alpha*beta);
ck=yk+(1-a-delta);
iy=(a+delta-1)/yk;
cy=ck/yk;
xsi=1.0-(beta*(a^(-eta)))*(1.0-delta);
theta = (1.0-alpha)*(yk/ck)*(1.0-nstar)*(1.0/nstar);

@ Third Part @

Cu[1,1]=-eta;       Cu[1,2]=-(theta*(1.0-eta)*nstar)/(1.0-nstar);
Cu[2,1]=1.0-eta;    Cu[2,2]=-(theta*(1.0-eta)-1.0)*(nstar/(1.0-nstar));                        Cu[2,5]=-1.0;
                    Cu[3,2]=alpha;                                                             Cu[3,5]=1.0;                                         Cu[3,5]=1.0;
                    Cu[4,2]=1.0-alpha;                                                                       Cu[4,6]=-1.0;                                        Cu[4,6]=-1.0;
                    Cu[5,2]=alpha-1.0;                                   Cu[5,3]=1.0;
Cu[6,1]=cy;                                                              Cu[6,3]=-1.0;  Cu[6,4]=iy; 

Cxl[1,2]=1.0;
Cxl[2,2]=1.0;
Cxl[3,1]=alpha;
Cxl[4,1]=1.0-alpha;
Cxl[5,1]=alpha;

Cz[3,1]=1.0;
Cz[4,1]=-1.0;
Cz[5,1]=1.0;

Dxl[1,1]=a;
Dxl[2,2]=1.0;

Fxl[1,1]=delta-1.0;
Fxl[2,2]=-1.0;

Du[2,6]=-xsi;

Fu[1,4]=a+delta-1.0;

_Rho[1,1]=RhoZ;
_Sigma[1,1]=SigmaZ;


/************************************************************************************/

/* Solve for the policy functions  */
output file=^outfile reset;

{Lxx,Lxz,Llx,Llz,Lux,Luz}=SolveLA(Cu,Cxl,Cz,Dxl,Fxl,Du,Fu,Dz,Fz,_rho,_nx);

output off;

/* The next part must be modified by the user:
**
** Set ut the matrices Lyx and Lyz that define
** the policy function for those variables of interest
** that are neither in ut, xt, or lamt.
**
** for example if there are no additional variables:
*/

Lyx=0;
Lyz=0;

/****************************************************************************************************************/
/****************************************************************************************************************/

/* Compute Impulse Responses */
if _IR;
struct PString {string name;};
struct PString pstr;

for i (1,_nz,1);
    {ut,xt,lamt,yt,zt}=Impulse(nobs1,i);
    fig="Figure" $+ ftos(i,"*.*lf",1,0) $+ ".tkf";    
    Plott(fig);
endfor;
endif;

/* Simulate the model and compute and output second moments */
if not _Moments; goto ende; endif;

{sx,rx}=RBCRun(nobs,nofs,HPL);

/* Output Results to the file outfile */

sml=GetMaxStringLength;
ccr=GetCrossCorrelations;
csx=GetRelSx;
if not ismiss(ccr); nccr=rows(ccr); else; nccr=0; endif;
if not ismiss(csx); ncsx=rows(csx); else; ncsx=0; endif;

cls;
outwidth 250;
output on;   
?"";
?"Second moments from simulated data:";
?"";
for i (1,_nvar,1);
    if Var[i].varprint;
        vi=GetVarIndex(i);
        Ausstring=strtruncpad(Var[i].varname,sml) $+ ftos(sx[vi],"*.*lf",8,2);
         for j (1,ncsx,1);
                Ausstring=Ausstring $+ ftos(sx[vi]/sx[csx[j,2]],"*.*lf",8,2);
         endfor;
         for j (1,nccr,1);
                Ausstring=Ausstring $+ ftos(rx[vi,ccr[j,2]],"*.*lf",8,2);
         endfor;        
         Ausstring=Ausstring $+ ftos(rx[vi,_nvar+vi],"*.*lf",8,2);
        ?Ausstring;
    endif;
endfor;
?"";
?"Column 1: Variable name";
?"Column 2: Standard deviation";
for i (1,ncsx,1);
        Ausstring="Column " $+ ftos(2+i,"*.*lf",1,0) $+ ": Standard deviation relative to standard deviation of variable " $+ Var[csx[i,1]].varname;
        ?Ausstring;
endfor;
for i (1,nccr,1);
        Ausstring="Column " $+ ftos(2+ncsx+i,"*.*lf",1,0) $+ ": Cross correlation with " $+ Var[ccr[i,1]].varname;
        ?Ausstring;
endfor;
?"Column " ftos(3+ncsx+nccr,"*.*lf",1,0) ": First order autocorrelation";

output off;
ende:
end;

/* This procedure uses the policy matrices from the main part
** simulates the model and returns averages of second moments
*/

proc(2)=RBCRun(nobs,nofs,hpl);

   local xt, zt, ut, lamt, st, itn, epst, i, t,
         sx, rx, mom1,mom2, yt;
   
   sx=zeros(_nvar,1);   
   rx=zeros(2*_nvar,2*_nvar);
   
  /* Simulation step */
  DosWin;
  itn=1;cls; locate 5,5;?"Simulation Nr. ";
  do until itn>nofs;
  locate 6,5;ftos(itn,"*.*lf",4,0);

  /* Initialize vectors */

  ut=zeros(_nu,nobs);
if _nx>0;  xt=zeros(_nx,nobs+1); endif;
lamt=zeros(_nl,nobs);
  if _ny>0; yt=zeros(_ny,nobs); endif;
  zt=zeros(_nz,nobs+1);
epst=rndn(_nz,nobs+1);

  /* Compute time path */
  zt[.,1]=_Sigma*epst[.,1];
  t=1;
  do until t>nobs;
    zt[.,t+1]=_Rho*zt[.,t]+_Sigma*epst[.,t+1]; 
    if _nx>0; 
         xt[.,t+1]=Lxx*xt[.,t]+Lxz*zt[.,t];
         ut[.,t]  =Lux*xt[.,t]+Luz*zt[.,t];
       lamt[.,t]  =Llx*xt[.,t]+Llz*zt[.,t];
        if _ny>0;
            yt[.,t]  =Lyx*xt[.,t]+Lyz*zt[.,t];
        endif;
    else;
        ut[.,t]=Luz*zt[.,t];
      lamt[.,t]=Llz*zt[.,t];
        if _ny>0;
            yt[.,t]=Lyz*zt[.,t];
        endif;
    endif;
      t=t+1;
  endo;

  /* Filtering */
  st=zeros(nobs,_nvar);
  for i (1,_nu,1);
          if HPL/=0; st[.,i]=HPFilter(ut[i,.]',HPL); else; st[.,i]=ut[i,.]';endif;
  endfor;
  for i (1,_nx,1);
          if HPL/=0; st[.,_nu+i]=HPFilter(xt[i,1:nobs]',HPL); else; st[.,_nu+i]=xt[i,1:nobs]';endif;
  endfor;
  for i (1,_nl,1);
          if HPL/=0; st[.,_nu+_nx+i]=HPFilter(lamt[i,.]',HPL); else; st[.,_nu+_nx+i]=lamt[i,.]';endif;
  endfor;
  for i (1,_ny,1);
         if HPL/=0; st[.,_nu+_nx+_nl+i]=HPFilter(yt[i,.]',HPL); else; st[.,_nu+_nx+_nl+i]=lamt[i,.]';endif; 
  endfor;
  for i (1,_nz,1);
        if HPL/=0;  st[.,_nu+_nx+_nl+_ny+i]=HPFilter(zt[i,1:nobs]',HPL); else; st[.,_nu+_nx+_nl+_ny+i]=zt[i,1:nobs]';endif;
  endfor;

  /* Computation of second moments */
  mom1=vcx(st);
  mom2=corrx(trimr(st,1,0)~trimr(st,0,1));  
  sx=sx+sqrt(diag(mom1));
  rx=rx+mom2;
  
  itn=itn+1;

endo;
DosWinCloseAll;

retp((sx/nofs)*100,rx/nofs);

endp;


/* This procedure computes impulse reponse functions. It uses
** the policy functions in Lxx, Lxz, ... .
** Impulses are to a one standard deviation shock in shock # i.
**
** The shock hit in period 2.
*/
proc(5)=Impulse(nobs1,i);

    local ut,xt,lamt,zt, yt, t, evec;

    if _nx>0; xt=zeros(_nx,nobs1+1); else; xt=0; endif;
    ut=zeros(_nu,nobs1);
  lamt=zeros(_nl,nobs1);
    zt=zeros(_nz,nobs1+1);
    if _ny>0; yt=zeros(_ny,nobs1); else; yt=0; endif;
    evec=zeros(_nz,1);
    evec[i]=1;

    zt[.,2]=_Sigma*evec;

    for t (2,nobs1,1);
        zt[.,t+1]=_Rho*zt[.,t];
        if _nx>0;
            xt[.,t+1]=Lxx*xt[.,t]+Lxz*zt[.,t];
            ut[.,t]  =Lux*xt[.,t]+Luz*zt[.,t];
        lamt[.,t]  =Llx*xt[.,t]+Llz*zt[.,t];
        if _ny>0;
            yt[.,t]  =Lyx*xt[.,t]+Lyz*zt[.,t];
        endif;
        else;
            ut[.,t]=Luz*zt[.,t];
          lamt[.,t]=Llz*zt[.,t];
            if _ny>0;
                yt[.,t]=Lyz*zt[.,t];
            endif;
        endif;
    endfor;
    
    retp(real(ut),real(xt),real(lamt),real(yt),real(zt));
endp;

/* Plott: This procedure takes ut, xt, and lamt from the
** main program and plots the impulse respones.
*/
proc(0)=Plott(fig);

    local pmax, vmax, pmat, ord, ymin, ymax, ymaxall, yminall, i, j, k, t, pout;        

    {pmax,vmax}=GetMaxPlottNumber;
    ord=(pmax|rows(ut')|maxc(vmax));    
    ymax=zeros(pmax,1);
    ymin=zeros(pmax,1);
    k=zeros(pmax,1);    
    pstr=reshape(pstr, pmax, 1);

    pmat=arrayinit(ord,0);
    for i (1,_nvar,1);
        if Var[i].varplot;
            j=var[i].plotno;
            k[j]=k[j]+1;
            if Var[i].vartype=="u";
                pmat[j,.,k[j]]=ut[Var[i].varpos,.]';
                pstr[j].name=pstr[j].name $+ Var[i].varname $+ "\000";                                        
            endif;
            if Var[i].vartype=="x";
                pmat[j,.,k[j]]=xt[Var[i].varpos,1:cols(xt)-1]';
                pstr[j].name=pstr[j].name $+ Var[i].varname $+ "\000";                        
            endif;
            if Var[i].vartype=="l";
                pmat[j,.,k[j]]=lamt[Var[i].varpos,.]';
                pstr[j].name=pstr[j].name $+ Var[i].varname $+ "\000";                        
            endif;
            if Var[i].vartype=="z";
                pmat[j,.,k[j]]=zt[Var[i].varpos,1:cols(zt)-1]';
                pstr[j].name=pstr[j].name $+ Var[i].varname $+ "\000";                        
            endif;
            if Var[i].vartype=="y";
                pmat[j,.,k[j]]=yt[Var[i].varpos,.]';
                pstr[j].name=pstr[j].name $+ Var[i].varname $+ "\000";
            endif;
        endif;
    endfor;
    pmat=pmat*100;
    for i (1,pmax,1);
        ymax[i]=maxc(maxc(arraytomat(pmat[i,.,.])));
        ymin[i]=minc(minc(arraytomat(pmat[i,.,.])));
    endfor;
    ymaxall=maxc(maxc(ymax));
    yminall=minc(minc(ymin));
    GraphSettings;
    _paxht=0.20;
    _pnumht=0.18;
    _ptek=fig;
    _plwidth=7;    
    if _nocolor; _pcolor={0,0,0,0,0,0,0}; _pltype={6,3,5,4,2,1}; endif;
    j=rows(ut');
    if pmax==1;
        _plegstr=pstr[1].name;    
        t=seqa(1,1,j);
        _pline=(1~3~0~0~maxc(t)~0~1~0~7);
        SetLegend(j,yminall,ymaxall);
        xlabel("Period");
        ylabel("Percentage Deviation from Trend");        
        pout=arraytomat(pmat[1,.,.]);
        xy(t,pout[.,1:vmax[1]]);
    endif;
    if pmax==2;
        begwind;
        window(2,1,0);
        t=seqa(1,1,j);
        _pline=(1~3~0~0~maxc(t)~0~1~0~7);
        if _scale; SetLegend(j,yminall,ymaxall); else; SetLegend(j,ymin[1],ymax[1]); endif;
        _plegstr=pstr[1].name;
        xlabel("Period");
        ylabel("Percentage Deviation from Trend");
        pout=arraytomat(pmat[1,.,.]);        
        if _scale;  scale((0|t),(yminall|ymaxall)); endif;
        xy(t,pout[.,1:vmax[1]]);
        nextwind; 
        if not _scale; SetLegend(j,ymin[2],ymax[2]); endif;
        _plegstr=pstr[2].name;
        pout=arraytomat(pmat[2,.,.]);
        xy(t,pout[.,1:vmax[2]]);
        endwind;
    endif;
    if pmax==3;
        begwind;
        window(2,2,0);
        t=seqa(1,1,j);
        _pline=(1~3~0~0~maxc(t)~0~1~0~7);
        if _scale; SetLegend(j,yminall,ymaxall); else; SetLegend(j,ymin[1],ymax[1]); endif;
        _plegstr=pstr[1].name;
        xlabel("Period");
        ylabel("Percentage Deviation from Trend");
        pout=arraytomat(pmat[1,.,.]);
        if _scale; scale((0|t),(ymin|ymax)); endif;
        xy(t,pout[.,1:vmax[1]]);
        nextwind;        
        if _scale; SetLegend(j,yminall,ymaxall); else; SetLegend(j,ymin[2],ymax[2]); endif;
        _plegstr=pstr[2].name;
        pout=arraytomat(pmat[2,.,.]);
        xy(t,pout[.,1:vmax[2]]);
        nextwind;
        _plegstr=pstr[3].name;
        if _scale; SetLegend(j,yminall,ymaxall); else; SetLegend(j,ymin[3],ymax[3]); endif;
        xlabel("Period");
        ylabel("Percentage Deviation from Trend");
        pout=arraytomat(pmat[3,.,.]);
        xy(t,pout[.,1:vmax[3]]);
        endwind;
    endif;    

   if pmax==4;
        begwind;
        window(2,2,0);
        t=seqa(1,1,j);
        _pline=(1~6~0~0~maxc(t)~0~1~0~2);
        if _scale; SetLegend(j,yminall,ymaxall); else; SetLegend(j,ymin[1],ymax[1]); endif;
        _plegstr=pstr[1].name;
        xlabel("Period");
        ylabel("Percentage Deviation from Trend");
        pout=arraytomat(pmat[1,.,.]);
        if _scale;  scale((0|t),(ymin|ymax)); endif;
        xy(t,pout[.,1:vmax[1]]);
        nextwind;
        if not _scale; SetLegend(j,ymin[2],ymax[2]); endif;
        _plegstr=pstr[2].name;
        pout=arraytomat(pmat[2,.,.]);
        xy(t,pout[.,1:vmax[2]]);
        nextwind;
        if not _scale; SetLegend(j,ymin[3],ymax[3]); endif;
        _plegstr=pstr[3].name;
        xlabel("Period");
        ylabel("Percentage Deviation from Trend");
        pout=arraytomat(pmat[3,.,.]);
        xy(t,pout[.,1:vmax[3]]);
        nextwind;
        if not _scale; SetLegend(j,ymin[4],ymax[4]); endif;
        _plegstr=pstr[4].name;
        xlabel("Period");
        ylabel("Percentage Deviation from Trend");
        pout=arraytomat(pmat[4,.,.]);
        xy(t,pout[.,1:vmax[4]]);
        endwind;
    endif;    


    
endp;



/* The remaining procedures are used to output data. They are
** of no special interest and therefore I do not comment them.
*/

proc(2)=GetMaxPlottNumber;
    local n0, n1, n2;
    n0=zeros(4,1);
    n1=0;
    for i (1,_nvar,1);
        if Var[i].varplot;
            n0[var[i].plotno]=n0[var[i].plotno]+1;
            n2=var[i].plotno;
            if n2>n1;
                n1=n2;
            endif;
        endif;
    endfor;
    n2=maxc(n0);    
    retp(n1,n0);
endp;


proc(1)=GetMaxStringLength;

    local s1, s2, i;

    s1=strlen(Var[1].varname);
    for i (2,_nvar,1);
        s2=strlen(Var[i].varname);
            if s2>s1; s1=s2; endif;
    endfor;

    retp(s1);

endp;

proc(1)=GetCrossCorrelations;

    local i, m;
    m=(0~0);
    for i (1,_nvar,1);
        if Var[i].crosscorr; 
            if Var[i].vartype=="u"; m=m|(i~Var[i].varpos); endif;
            if Var[i].vartype=="x"; m=m|(i~(_nu+Var[i].varpos)); endif;
            if Var[i].vartype=="l"; m=m|(i~(_nu+_nx+Var[i].varpos)); endif;
            if Var[i].vartype=="y"; m=m|(i~(_nu+_nx+_nl+Var[i].varpos)); endif;
            if Var[i].vartype=="z"; m=m|(i~(_nu+_nx+_nl+_ny+Var[i].varpos)); endif;
        endif;
    endfor;
    if rows(m)>1;
        retp(m[2:rows(m),.]);
    else;   
        retp(miss(1,1)~miss(1,1));
    endif;

endp;

proc(1)=GetRelSx;

    local i, m;
    m=(0~0);
    for i (1,_nvar,1);
        if Var[i].relsx; 
            if Var[i].vartype=="u"; m=m|(i~Var[i].varpos); endif;
            if Var[i].vartype=="x"; m=m|(i~(_nu+Var[i].varpos)); endif;
            if Var[i].vartype=="l"; m=m|(i~(_nu+_nx+Var[i].varpos)); endif;
            if Var[i].vartype=="y"; m=m|(i~(_nu+_nx+_nl+Var[i].varpos)); endif;
            if Var[i].vartype=="z"; m=m|(i~(_nu+_nx+_nl+_ny+Var[i].varpos)); endif;
        endif;
    endfor;
    if rows(m)>1;
        retp(m[2:rows(m),.]);
    else;
        retp(miss(1,1)~miss(1,1));
    endif;

endp;

            
proc(1)=GetVarIndex(i);

    local vi;

    if Var[i].vartype=="u"; vi=Var[i].varpos;             endif;
    if Var[i].vartype=="x"; vi=_nu+Var[i].varpos;          endif;
    if Var[i].vartype=="l"; vi=_nu+_nx+Var[i].varpos;       endif;
    if Var[i].vartype=="y"; vi=_nu+_nx+_nl+Var[i].varpos;    endif;
    if Var[i].vartype=="z"; vi=_nu+_nx+_nl+_ny+Var[i].varpos; endif;

    retp(vi);

endp;

    
proc(0)=SetLegend(j,ymin,ymax);

    external matrix _plegctl;
 
    if abs(ymin)>ymax;
        _plegctl=(1~7~floor(j/3)~ymin);
    else;
       _plegctl=(1~7~floor(j/3)~0.5*ymax);
    endif;       
    
retp;

endp;
