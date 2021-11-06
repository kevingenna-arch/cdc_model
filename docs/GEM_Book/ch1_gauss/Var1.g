@ -------------------------- Var1.g ---------------------------------

  Author:    Alfred Maussner
  Date:      31 December 2007
  Purpose:   Estimate Var in x=[y,c,i,h] using 
             HP-filtered German quartely data on output (y),
             private consumption (c), investment expenditure (i),
             and working hours (h)

--------------------------------------------------------------------- @

new;

/* Read Data from ascii file */
load data[120,37]=data.txt;


/* Restrict data to the period 75.i bis 89.iv  */
@data=data[61:120,.];@

/* Convert nominal to real data */
pbip=data[.,15]/100;
pinvest=(data[.,7]./data[.,4]);

 bt=data[.,22];                        @ Wohnbevoelkerung in 1000 Personen    @
  yt=data[.,2]-(data[.,23]./pbip);     @ Reales BIP zu Faktorkosten in Mrd DM @
 ytk=(yt./bt)*1e6;                     @ Reales BIP zu Faktorkosten in DM pro Kopf und Quartal @
  it=data[.,4];                        @ Reale Anlageinvestitionen in Mrd DM  @
 itk=(data[.,4]./bt)*1e6;              @ Reale Anlageinvestitionen in DM pro Kopf und Quartal @
  ct=data[.,3];                        @ Realer privater Verbrauch in Mrd DM  @
 ctk=(data[.,3]./bt)*1e6;              @ Realer privater Verbrauch in DM pro Kopf und Quartal @
 ytn=(data[.,2].*pbip)-data[.,23];     @ Nominales BIP zu Faktorkosten in Mrd DM @
  ht=data[.,17];                       @ Arbeitsstunden in Mill. Stunden @
 htk=(data[.,17]./bt)*1000;            @ Arbeitsstunden in Stunden pro Kopf @
 wt =(data[.,16]./data[.,36])*1000;    @ DM je Arbeitsstunde der abhaengig Beschaeftigten               @
 wt = (100*wt)./data[.,37];            @ Kaufkraft DM je Arbeitsstund der abhaenig Beschaeftigten       @


/* Filtering */
data1=ytk~ctk~itk~htk;
@data1=yt~ct~it~ht;@
data2=zeros(rows(data1),cols(data1));

for i (1,4,1);
    data2[.,i]=HPFilter(ln(data1[.,i]),1600);
endfor;


/* Specification of the Var */
cls;
_varTrend =0;        @ Do not include a linear trend in an estimation of a VAR @
_Eviews   =1;        @ compute AIC and SIC using the formulas given in the Eviews 4.0 help file @
_varEx    =0;        @ set to nobs x l matrix, if dummies are to be included in VAR procedures  @

crit=VarOrder(data2,8,1);

?"Input Order of Var:";
p=con(1,1);
p=2;        @ number of lags @
nvar=4;     @ number of variables in the VAR @
var="Output"|"Consump"|"Investm"|"Hours";  @ Names of variables @
nir=10;     @ number of periods for which impulse responses are to be computed @
nofs=10000;  @ number of simulations used to construct confidence bounds @

/* Estimation and identification of impulse respones via Cholesky factroization */
_Var_Quick=0;   @ Display estimation results @

output file=Var1.txt reset;

{b0,e0,ir0}=SVar1(data2,p,nir,var);

output off;

/* Compute confidence bounds for the estimated impulse responses */
_Var_Quick=1; @ Do not display and compute additional stuff @

{irl,iru}=GetBounds_SVar(data2,b0,e0,p,nir,nofs); 
   
GraphSettings;
_paxht=0.20;
_pnumht=0.18;
 scale((0|nir),(-0.002|0.02));
_pline=(1~6~0~0~30~0~1~0~2);
_pcolor={0,0,0};
_pltype={6,3,3};
_plwidth={7,4,4};
_ptitlht=0.20;
_ptek="Fig1.tkf";
t=seqa(1,1,nir);
begwind;
window(2,2,0);
title("Output");
xlabel("Period");
    ir=arraytomat(ir0[1,.,.]);
    irl=arraytomat(irl[1,.,.]);
    iru=arraytomat(iru[1,.,.]);
    xy(t,ir[1:nir,1]~irl[1:nir,1]~iru[1:nir,1]);
nextwind;
title("Consumption");    
    xy(t,ir[1:nir,2]~irl[1:nir,2]~iru[1:nir,2]);
nextwind;
title("Investment");
    xy(t,ir[1:nir,3]~irl[1:nir,3]~iru[1:nir,3]);
nextwind;
title("Hours");
   xy(t,ir[1:nir,4]~irl[1:nir,4]~iru[1:nir,4]);
endwind;

goto ende;

GraphSettings;
_ptek="Fig2.tkf";
t=seqa(1,1,39);
begwind;
window(2,2,1);
_plegstr="Output\000Consumption\000Investment\000Hours";
_plegctl=(1~5~7~0.01);
for k (1,cols(data2),1);
    text="Response to shock in variable no. " $+ ftos(k,"*.*lf",1,0);
    title(text);
    xlabel("Period");
    ir1=arraytomat(ir[k,.,.]);    
    xy(t,ir1[2:40,.]);
    nextwind;
    _plegctl=0;
endfor;
endwind;

ende:

end;


@ -------------------------- Procedures ------------------------------ @

@ ---------------------------------- SVar1 --------------------------------

   Author: Alfred Maussner
   
   Purpose: Estimation of impulse response functions from
            a p-order stationary vector autoregressive model via
            the Cholesky decomposition.
            
            The model is:       
         
            y[t] = c +phi1*y[t-1] + ... + phip*y[t-p]+eps[t]
         
            y[t] is a nvar x 1-vector with observations on nvar I(0)-variables
            at time t, p is the maximal lag considered.
         
         
   Remark:  The code draws on the Gauss-Code from var.g and m2var.prg,
            repectively, from Cochrane (1998).
             
             
   Usage:    {b,e,ir,irl,iru}=MyVarIC(y,p,nir,var)
  
  Input:
  
         y   := nobs x nvar-data matrix, whose nvar columns hold the nobs
                observations on the nvar variables.
                       
                    
        p    := scalar, the number of lagged differences on the right hand
                side of the VAR
                
                
       nir   := scalar, number of periods for which impulse responses 
                are computed
                
       var   := nvar x 1 vector of variable names (8 characters!)
                               
  Output:
  
         ir := nvar by nir by nvar matrix. The ith page of
               this three-dimensional array stores in its jth
               column the impulse response of the jth variable
               due to an orthogonal shock of one standard deviation in
               variable i.

        irl := nbar by nir by nvar matrix. The same as above,
               but the interpretation is now: the lower 2.5 percent
               percentil of the impulse response from 10,000 simulations.

        iru := nbar by nir by nvar matrix, the upper 97.5 percent percentil
                                        
------------------------------------------------------------------------------ @

proc(3)=SVar1(data,p,nir,var);

  local nov, nobs, y,  x, i, j, k, b, e, sigma, stderr, rsq, r, ysim,
        vshock, out, fmt, mask, left, head, fmat, lambda;

  external matrix _Var_Quick;
  
  /* Get important parameters */
  nov=cols(data);    @ number of variables       @
  nobs=rows(data);   @ number of observations    @


  /* right hand side of the VAR */
  y=data[1+p:nobs,.];
  
  x=ones(nobs-p,1);     @ constant @  
  for j (1,p,1);
     for i (1,nov,1);     
        x=x~data[1+p-j:nobs-j,i];
     endfor;   
  endfor;   
       
 /* Estimation: for each equation is
 ** x'y[.,j]/(x'x)=b[j]=y/x the OLS estimator. if y is a matrix,
 ** GAUSS returns the matrix b, whose jth column holds the
 ** b[j].
 */

     b = y/x;
 
     e = y-x*b;
 sigma = (e'*e)/(nobs-p);
stderr = (diag(invpd(x'*x))*diag(sigma)')^(1/2);
   rsq = 1 - diag(sigma)./(stdc(y)^2);

/* Transformation of the VAR(p) in a VAR(1) model to check stability */
if not _Var_Quick;
    fmat=Var_Amat(b,p);
    lambda=rev(sortc(eig(fmat),1));

    gosub prnt; @ print resulsts @

endif;

 /* orthogonal residuals */
 
 r = chol(sigma)';

 /* Iterate the model to get impulse responses,
 ** where e=r*v and variable k is displaced by one standard deviation */
     ysim = arrayinit(nov|nir|nov,0);

for k (1,nov,1);  @ loop over shocks @
   
         x   = zeros(p*nov,1);   @ right hand side, without constant @ 
    vshock   = zeros(nov,1);
 vshock[k]   = 1;                       
 ysim[k,2,.] =  (r*vshock)';     @ one standard deviation shock   @
 
   for i (3,nir,1);              @ iterate over the estimated equation @
     x[1+nov:p*nov]=x[1:(p-1)*nov];
     for j (1,nov,1);            @ neue verzoegerte Veraenderungen einfuegen@     
         x[j]=ysim[k,i-1,j];
     endfor;
     ysim[k,i,.] = (b[2:rows(b),.]'x)';     @ neue Veraenderungen aus dem VAR berechnen @
   endfor;   
   
endfor;

retp(b,e,ysim);

prnt:

?"Results of OLS estimation";
left="Const.";
head="      ";
j=1;
for j (1,nov,1);
   head=head $+ "           "$+var[j];
endfor;

for j (1,p,1);
    for i (1,nov,1);    
           left=left|"Phi"$+ftos(i,"*.*lf",1,0)$+ftos(j,"*.*lf",1,0);
    endfor;
endfor;
?"";
?"Equation: ";
$head;?"";
out=left~b;
 mask=0~ones(1,cols(b));
 fmt={"*.*s" 8 8};
 for i (1,nov,1); fmt=fmt|("*.*lf"~14~3);endfor;

 i=printfm(out,mask,fmt);
 ?"";
 ?"eigen values of the VAR(1)-representation:";?"";
 for i (1,p*nov,1);
   ?ftos(lambda[i],"*.*lf",10,4)$+"  Modulus: "$+ftos(abs(lambda[i]),"*.*lf",10,4);
 endfor;

return;

endp;    

@ ------------------------------- Var_Amat -------------------------------------------------

   Purpose: Builds the matrix of the AR(1) representation of a VAR(p) model:

            y[t] = c  +  Phi[1]*y[t-1] + ... + Phi[p]*y[t-p1].

            The AR(1) representation of this model is build from

            y[t]-mu = A*(y[t-1]-mu)

            A= Phi[1] Phi[2] .... Phi[p-1] Phi[p]
               In     0      .... 0        0
               .      .      .... .        .
               0      0      .... In       0

            The eigenvalues of A must be all inside the unit ball. 
             
   Usage:   A=Svar_Amat(b,p);

   Input:   b :=  1+p*K times K matrix, each column holds the
                  estimated coefficients of the k-th equation.

           p  := the number of lags (not lagged differences!) in the VEC

   Output: A  := p*K times p*K matrix (defined above)

--------------------------------------------------------------------------------------- @

proc(1)=Var_Amat(b,p);

  local k,i,a;

  k=cols(b);

  a=zeros(p*k,p*k);

  for i (1,k,1);  
     a[i,.]=b[2:1+k*p,i]';
  endfor;
  for i (1,p-1,1);  
     a[i*k+1:(i+1)*k,1+(i-1)*k:i*k]=eye(k);
  endfor;

retp(a);

endp;  


@ --------------------------- Resample ----------------------------------------

   Author: Alfred Maussner
   Date:   12. December 2007
   
   Purpose: Get a subsample from a given data matrix

   Usage:  e1=Resample(e);

   Input:  e:  nobs x nvar data matrix

   Output: e1: nobs x nvar data matrix whose jth column
               is a random draw with replacement
               form the jth column of e0.
   Remark: In a former version, this was not done row-wise.
           (See the article on monetary policy shocks in
           the handbook of macroeconomics, vol. IA, chapter 2, footnote 23
           on this kind of resampling)

--------------------------------------------------------------------------- @

proc(1)=Resample(e);

   local nobs, extract;

   nobs=rows(e);
   extract=trunc(rndu(nobs,1)*nobs)+ones(nobs,1);

 retp(e[extract,.]);

endp;


@ ----------------------------------- MakeData ------------------------------ 

    Author: Alfred Maussner
    Date:   12 December 2007

    Purpose: create artifical time series data from an estimated
             Var(p) model with shocks from resampled estimated
             residuals.

    Usage:   Data=MakeData(data,b,e,p)

    Input:  data : nobs by nvar matrix of observations
            b    : nvar*p+1 by nvar matrix of estimated coefficients
            e    : nobs by nvar matrix of new shocks
            p    : number of lags

------------------------------------------------------------------------------- @


proc(1)=MakeData(data,b,e,p);

  local ysim, nvar, nobs,xvec, j, t;
  
  nvar=cols(data);
  nobs=rows(data);
  
  ysim =  zeros(nobs,nvar);     @ simulated data               @
 
  ysim[1:p,.]=data[1:p,.];      @ the first p values are treated as presample values @

  xvec=zeros(1,1+nvar*p);
  xvec[1]=1;
  for j (1,p,1);
    xvec[2+(j-1)*nvar:1+nvar*j]=data[p-(j-1),.];
  endfor;
  for t (1,nobs-p,1);  
        ysim[p+t,.]=xvec*b+e[t,.];
        xvec[2+nvar:cols(xvec)]=xvec[2:cols(xvec)-nvar];
        xvec[2:1+nvar]=ysim[p+t,.];
  endfor;
   
 retp(ysim);
  
endp;

@ -------------------------- GetBounds_SVar -------------------------------

    Author: Alfred Maussner
    Date:   12 December 2007
   
    Purpose: Compute confidence bounds for the estimate impulse responses
             from SVar1.

    Usage: {irl,iru}=GetBounds_SVar(data,b,e,p,nir,nofs);

    Input:   data  : nobs by nvar matrix with nobs obseravations on nvar
                     variables.
             b     : 1+nvar*p by nvar matrix, estimated coefficients from
                     SVar1
             e     : nobs by nvar, residuals from SVar1
             p     : number of lags
           nir     : number of periods for which impulse responese are to
                     be computed
           nofs    : number of simulations 

    Output: irl    : nvar by nir by nvar array. page k stores the
                     lower 2.5 bound of the impulse response
                     of variables 1 through nvar due to a shock
                     in the equation of variable k

            iru    : same as irl, yet with upper 97.5 bound.

----------------------------------------------------------------------------- @

proc(2)=GetBounds_SVar(data,b,e,p,nir,nofs);

    local s, n_var, bounds, data1, b1, e1, ir1, il, iu, irl, iru, k, i, j, temp;
    external matrix _Var_Quick;

    n_var=cols(data);

    bounds=arrayinit(nofs|n_var|nir|n_var,0);

    for s (1,nofs,1);
                     e1=Resample(e);
                  data1=MakeData(data,b,e1,p);
            {b1,e1,ir1}=SVar1(data1,p,nir,0);  
        bounds[s,.,.,.]=ir1;
    endfor;

    il=trunc(0.025*nofs);
    iu=trunc(0.975*nofs);
    irl=arrayinit(n_var|nir|n_var,0);
    iru=arrayinit(n_var|nir|n_var,0);

    for k (1,nvar,1);                 @ loop over shocks    @
        for j (1,nvar,1);             @ loop over variables @
            for i (2,nir,1);          @ loop over periods   @
                temp=zeros(nofs,1);
            for s (1,nofs,1);         
                temp[s]=bounds[s,k,i,j];  @ store period i response of  variable j to shock k in simulation s for sorting @
            endfor;
            temp=sortc(temp,1);
            irl[k,i,j]=temp[il];
            iru[k,i,j]=temp[iu];
        endfor;
    endfor;
endfor;
    
retp(irl,iru);

endp;

@ ---------------------------- VarOrder -----------------------------

  Author: Alfred Maussner
  
  Purpose: Compute VAR-order selection criteria AIC, HQ, and SIC
  
  call: c=VarOrder(data,m,p);
  
  Input: data := nobs x nvar matrix, whose columns hold the nobs
                 observations for nvar variables.
              
            m := scalar, maximum var-length considered
    
            p := logical, p==1: print results, p==0 don't
         
  Output: c := m x 3 vector, whose first, second, and third
                     column return for lag length
                     i=1, 2, ..., m the AIC, HQ,and SIC criteria

  Globals: _varEx := nobs times k matrix, exogenous regressors, as, e.g., 
                                          dummy variables
          _EViews := 0,1, see below

          _varTrend:= 1 include linear time trend, 0 no linear time trend

                       
  Remarks:- The formulas are taken from  Lütkepohl (1993),
            p. 129 (equation 4.3.2 for AIC), p. 132 (equation
            4.3.8 for HQ, 4.3.9 for SIC).

          - These criteria can be employed to determine the
            order of a stationary or an unstationary VAR
            (See Luetkepohl (1993), p. 382f).

          - The AIC is not a consistent criterium, wheras the
            HQ and SIC are. Nevertheless, and in particular in
            small samples, the AIC may not be inferior to SIC and HQ

          - for 16 and more observations the Var order chosen
            by the three criteria obeys
           
            p(SC)<=p(HQ)<=p(AIC)

            (See Luetkepol, 1993, p.133 and p. 383)
           

          - Different form there, however, the number of
            estimated parameters K includes the constant and
            the trend term (when present).
            
          - The program was tested using the data in Table E1, p.498f in
            Luetkepohl (1993) and the results in Table 4.5, p. 130

            Using the same data, Eviews 4.0 output on AIC, HQ, and SIC
            is very different from the results in Luetkepohl (1993),
            Table 4.5, though the program manuel often refers to that book.
            Setting _Eviews=1, the program uses the formulas found
            in the Eviews 4.0 help file to compute AIC and SIC. This
            provides the output given by Eviews 4.0.

            The difference is twofold:

                  . Eviews includes the term nvar*(1+ln(2*pi))
                  . The number of estimated parameters is computed
                    including the constant and - if applicable - the
                    slope of a linear time trend parameter.

            See the program Test18a.g and the workfile Test18.wf1 to check these
            propositions.
            
          - The program uses the same number of observations
            for all p calculations of AIC, HQ, and SIC.
             
           
------------------------------------------------------------------ @

proc(1)=VarOrder(data,m,p);

 local i, j, k, k1, aic, sic, hq, nvar, nobs, x, y, b, e, sigma, logl;
 external matrix _varTrend, _Eviews, _varEx;

 aic=zeros(m+1,1);
 sic=aic;
 hq=aic;
 k1=0;

 nvar=cols(data);
 nobs=rows(data);
   
 @ set up the left hand side of the var @
 
 y=data[1+m:nobs,.];
    
 @ set up right hand side of the var @ 
    
 i=0;
 
 do until i>m;
 x=ones(rows(y),1); 
if _varTrend; x=x~seqa(0,1,rows(x)); endif;
if rows(_varEx)>1;
   x=x~_varEx[1+m:nobs,.];
  k1=cols(_varEx);
endif;

 k=1;
 do until k>nvar;
      j=1;
      do until j>i;
         x=x~data[1+m-j:nobs-j,k];
         j=j+1;
      endo;   
      k=k+1;
 endo;    

 /* Estimation: The OLS estimator of equation k is
 ** given by
 ** x'y[.,k]/(x'x)=b[k]=y[.,k]/x.
 ** if y is a matrix with nvar columns,
 ** the GAUSS command line
 ** b=y/x returns a cols(x) x cols(y) matrix b, whose
 ** k-th column holds the OLS coefficients of the
 ** k-th equation.
 */

 b = y/x;
 e = y-x*b;
 sigma = (e'e)/(rows(y));
 if _Eviews;
    aic[i+1]=ln(det(sigma))+nvar*(1+ln(2*pi));
    if _varTrend;
      k=nvar*(2+k1+nvar*i);
    else;
      k=nvar*(1+k1+nvar*i);
    endif;
 else;
    aic[i+1]=ln(det(sigma));
    if _varTrend;
       k=nvar*(1+k1+nvar*i);
    else;
       k=nvar*(k1+nvar*i);
    endif;
 endif;
  hq[i+1]=aic[i+1];
 sic[i+1]=aic[i+1];
 aic[i+1]=aic[i+1]+2*(k/rows(y));
  hq[i+1]=hq[i+1]+((2*k*ln(ln(rows(y))))/rows(y));
 sic[i+1]=sic[i+1] + (k*ln(rows(y)))/(rows(y));

 i=i+1;
 
 endo;
 if p;  
    i=minindc(aic)~minindc(hq)~minindc(sic);
    ?"Var Order Selection criteria:";
    ?"    AIC            HQ          SIC";
    ?"=====================================";
    ?"";
    j=1;
    do until j>m+1;
       if i[1]==j; k=ftos(aic[j],"*.*lf*",12,5);     else; k=ftos(aic[j],"*.*lf",12,5); endif;
       if i[2]==j; k=k $+ ftos(hq[j],"*.*lf*",12,5); else; k=k $+ ftos(hq[j],"*.*lf",12,5); endif;
       if i[3]==j; k=k $+ ftos(sic[j],"*.*lf*",12,5);else; k=k $+ ftos(sic[j],"*.*lf",12,5);endif;
       ?$k;
       j=j+1;
    endo;
  endif;
    
 retp(aic~hq~sic);
 
endp;

