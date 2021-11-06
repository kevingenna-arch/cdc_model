@ --------------------------------------- Ramsey2e.g ---------------------------------------------------------

   Alfred Maussner
   16. November 2003 
   Last revision: 05 October 2011

   Purpose: Solve the deterministic growth model (Ramsey2) from Chapter 6 of Heer/Maussner, 2nd edition,
            using projection methods. You may choose between least squares (in this case
            we minimize the least squares function integral instead of solving the respective system of
            first order conditions), Galerkin and Chebyshev collocation.
            
   External procedures:
            ChebCoef      from FUNCTION.src
            GSearch1      from SEARCH.src
            FixvMN1       from NLEQ.src                    
            ChebEval1     from FUNCTION.src 
            GraphSettings from TOOLS.src           

   Notes: For the program to run you must either include the names of the procedures
          in the files above in your Gauss user library (user.lcg) or in any other
          Gauss library file or use the statement #include statements to include the code in these files
          at the top of the current file (given that the files are in your current working directory)
          
          The *.src files are in Ch11_Gauss.zip.


------------------------------------------------------------------------------------------------------------- @

new;

/* Parameters */
alpha=0.27;
beta=0.994;
delta=0.011;
@delta=1;@
eta=2;
@eta=1;  @

/* Stationary Solution */
kstar=((alpha*beta)/(1-beta*(1-delta)))^(1/(1-alpha));
cstar=(kstar^alpha)-delta*kstar;

/* Parameters for the algorithm */
n=5;                @ degree of polynomial @
nint=100;           @ number of nodes in Gauss-Chebysehv quadratur @
ku=1.5*kstar;       @ upper bound of capital stock                 @
kl=0.5*kstar;       @ lower bound of capital stock                 @
d=(kl|ku);
kl_eu=0.8*kstar;
ku_eu=1.2*kstar;
n_eu=100;
method=1;        @ =1: Least Squares, 2: Galerkin, 3: Collocation @

/* Name of the output file */
outfile="d:\\tbk\\tb5_2\\g\\ch6\\Ramsey2e_1.txt";

/* Open Dos Window */
DosWinOpen("Ramsey2: Projection Methdos",0|0|15|1|1);
output file=^outfile on;
?"";
MyDate;
output off;

if (eta==1) and (delta==1);
   kvec=seqa(kl,(ku-kl)/99,100);
   cvec=(1-alpha*beta)*kvec^alpha;
   avec1=ChebCoef(&c0,n,100,d);  @ Initial Coefficients @
   asec=hsec;
   if method==1; goto LS; endif;
   if method==2; goto GA; endif;
   if method==3; goto CO; endif;

endif;
next:
/* Parameters for the search routine */
Npar=n;
Npop=50;
Ngen=100;

cls;
asec=hsec;
FError=0;
{avec1,fit}=GSearch1(Npar,Npop,Ngen,&f); @ you may want to use this procedure to find acceptable starting values @

LS:
cls;
if method==1;
  /* Refine Solution using our own quasi Newton routine */
  cls;
  _QN_Print=1;
  _QN_GradTol=1.e-10;
  {avec2,crit}=QuasiNewton(avec1,&f);
  output on;
  ?"Exit from Quasi Newton:";
  crit;
  ?"Parameters from the Least Squares Solution:";
  i=1;
  do until i>rows(avec2);
     ?ftos(avec2[i],"*.*lf",10,6);
     i=i+1;
  endo;
  output off; 
  
endif;
GA:
cls;
if method==2;
  /* Refine Solution using Galerkin */
  _MNR_Print=1;
  _MNR_Global=1;
  {avec2,crit}=FixvMN1(avec1,&Sys2);
  output on;
  ?"Exit from FixvMN1:";  
  crit;
  ?"Parameters from the Galerkin Solution:";
  i=1;
  do until i>rows(avec2);
     ?ftos(avec2[i],"*.*lf",10,6);
     i=i+1;
  endo;
  output off;
endif;
CO:
cls;
if method==3;
  _MNR_Print=1;
  _MNR_Global=1;
  {avec2,crit}=FixvMN1(avec1,&Sys1);
  output on;
  ?"Exit from FixvMN1:";
  crit;
  ?"Parameters from the Collocation Solution:";
  i=1;
  do until i>rows(avec2);
     ?ftos(avec2[i],"*.*lf",10,6);
     i=i+1;
  endo;
  output off;
endif;

output on;
?"";
?"Run Time: " etstr(hsec-asec);
?"";
output off;
plot:
if method==1;text="\000Least Squares";endif;
if method==2;text="\000Galerkin";     endif;
if method==3;text="\000Collocation";  endif;

/* Compare to true solution */
if (eta==1) and (delta==1);
   kvec =seqa(kl,(ku-kl)/99,100);
   cvec1=(1-alpha*beta)*(kvec^alpha);
   cvec2=ChebEval1(avec2,kvec,d);
   output on;
   ?"Fit:";
   if method==1;
     ?"Maximum of distance of least squares solution to analytic solution:";
   endif;
   if method==2;
     ?"Maximum of distance of Galerkin solution to analytic solution:";
   endif;
   if method==3;
     ?"Maximum of distance of collocation solution to analytic solution:";
   endif;
   load test;
   
   ftos(maxc(abs(cvec1-cvec2)),"*.*lf",10,7);
   ?"Maximum of distance to 2nd order solution";
   ftos(maxc(abs(test[.,2]-cvec2)),"*.*lf",10,7);
   output off;

   GraphSettings;
   xlabel("K");
   ylabel("C");
   _plwidth=7;
   _pcolor={0, 0};
   _pltype={6, 3};
   _pstype=8;
   _plctrl={0, 3};
   _psymsiz=2;   
   _plegstr="True Solution"$+text;
   _plegctl=(1~5~minc(kvec)~0.95*maxc(cvec2));
   xy(kvec,cvec1~cvec2);
 endif;

/* Compute Euler equation residuals */    
eu_kvec=seqa(kl_eu,(ku_eu-kl_eu)/(n_eu-1),n_eu);
eer=Euler(eu_kvec);
output on;
?"Maximum absolute value of Euler equation residual: ";
?ftos(maxc(abs(eer)),"*.*lf",15,10);
output off;

ende:

end;


/* --------------------------------------------------- Procedures ------------------------------------------------ */

/* c0(k): returns rows(k) times the stationary solution */
proc(1)=c0(k);
 local nobs, i, c;
 nobs=rows(k);
 c=zeros(nobs,1);
 i=1;
 do until i>nobs;
    c[i]=cvec[sumc(kvec.<=k[i])];
       i=i+1;
 endo;
 retp (c);
endp;

/* f(avec) computes the objective function, whose minimizer is the desired solution
**
** from the main program it uses nint, ku, kl, and d  */

proc(1)=f(avec);

   local i, sum, r, kt, k0, k1, c0, c1;

   i=1; sum=0;
   do until i>nint;
      kt = cos(((2*i - 1)/(2*nint))*pi);   @ i-th zero of T_m            @
      k0 = (1/2)*(kt+1)*(ku-kl) + kl;      @ transformed to k in [kl,ku] @
      c0 = ChebEval1(avec,k0,d);
      if c0<0; retp(miss(1,1)); endif;
      k1 = (k0^alpha) + (1-delta)*k0 - c0;
      if ((k1<kl) or (k1>ku)); retp(miss(1,1)); endif;      
      c1 = ChebEval1(avec,k1,d);
      if (c1<0); retp(miss(1,1));endif;
      r  = beta*((c1/c0)^(-eta))*(1-delta+alpha*(k1^(alpha-1)));
      r  = r-1;
      sum=sum+(r*r)*sqrt(1-(kt^2));

        i=i+1;
   endo;
   sum=sum*(pi*(ku-kl))/(2*nint);

  retp(sum);

endp;

/* Sys1(avec) returns the rhs of the system of equations, whose zero defines
** the collocation solution.
**
** From the main program it uses nint, ku, kl, and d
*/



proc(1)=Sys1(avec);

  local  p, r, i, kt, k0, k1, c0, c1;

    p=rows(avec);
    r=zeros(p,1);

   i=1; 
   do until i>p;
      kt   = cos(((2*i - 1)/(2*p))*pi);      @ i-th zero of T_p            @
      k0   = (1/2)*(kt+1)*(ku-kl) + kl;      @ transformed to k in [kl,ku] @
      c0   = ChebEval1(avec,k0,d);
      if c0<0; retp(miss(1,1)); endif;
      k1   = (k0^alpha) + (1-delta)*k0 - c0;
      if ((k1<kl) or (k1>ku)); retp(miss(1,1)); endif;      
      c1   = ChebEval1(avec,k1,d);
      if (c1<0); retp(miss(1,1)); endif;
      r[i] = beta*((c1/c0)^(-eta))*(1-delta+alpha*(k1^(alpha-1)));
      r[i] = r[i]-1;
         i = i+1;
   endo;

 retp(r);

endp;

/* Sys2(avec) returns the rhs of the system of equations, whose zero defines
** the Galerkin solution.
**
** From the main program it uses nint, ku, kl, and d 
*/

proc(1)=Sys2(avec);

   local i, j, p, sum, r, kt, k0, k1, c0, c1, t0, t1, t2;
   p=rows(avec);
   sum=zeros(p,1);

   j=1;
   t0=ones(nint,1);   
   kt=cos(((2*seqa(1,1,nint)-1)/(2*nint))*pi); @ the nint zeros of T_nint @   
   t1=kt;
   
   do until j>p;          
     i=1;

   do until i>nint;      
      k0    = (1/2)*(kt[i]+1)*(ku-kl) + kl;      @ transformed to k in [kl,ku] @
      c0    = ChebEval1(avec,k0,d);
      if c0<0; retp(miss(1,1)); endif;
      k1    = (k0^alpha) + (1-delta)*k0 - c0;
      if ((k1<kl) or (k1>ku)); retp(miss(1,1)); endif;      
      c1    = ChebEval1(avec,k1,d);
      if (c1<0); retp(miss(1,1));endif;
      r     = beta*((c1/c0)^(-eta))*(1-delta+alpha*(k1^(alpha-1)));
      r     = r-1;
      r     = r*t0[i]*sqrt(1-(kt[i]^2));
      sum[j]= sum[j]+r;
           i=i+1;
   endo;
   sum[j]=sum[j]*(pi*(ku-kl))/(2*nint);
   t2=2*(kt.*t1)-t0;
   t0=t1;
   t1=t2;
   j =j+1;
 endo;
   

 retp(sum);

endp;

/* eer=Euler(kvec): Computation of the Euler equation  residuals over the elements in kvec 
**
**                  The vector of coefficients avec2 and of bounds d must be given in the main program
*/

proc(1)=Euler(kvec);


    local n, eer, i, c0, c1, k1, rhs;

    n=rows(kvec);
    eer=zeros(n,1);
    
    for i (1,n,1);
        
        if (kvec[i]<d[1]) or (kvec[i]>d[2]); ?"error"; wait; endif;        
        c0=ChebEval1(avec2,kvec[i],d);        
        k1=(kvec[i]^alpha)+(1-delta)*kvec[i] - c0;        
        if (k1<d[1]) or (k1>d[2]); ?"error"; wait; endif;
        c1=ChebEval1(avec2,k1,d);
        rhs=beta*(c1^(-eta))*(1-delta+alpha*(k1^(alpha-1)));
        c1=rhs^(-1/eta);
        eer[i]=(c1/c0)-1;   
        
    endfor;
    
retp(eer);

endp;



