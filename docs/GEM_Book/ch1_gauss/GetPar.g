@ ------------------------------------ GetPar.g ----------------------------

   Alfred Maussner
   22 May 2000

   Last revision: 20 November 2008

   Purpose: Compute parameter values to calibrate the benchmark model
            of Heer and Maussner, 2nd edition.

 ------------------------------------------------------------------------  @

/* clear memory */
new;

/* Read data from file */
load data[120,37]=data.txt;

/* Restrict data to the period 75.i to 89.iv */
_short=1;  @ if zero, the full data set will be used @
if _short;
    data=data[61:120,.];
endif;

/* price index of GDP */
pbip=data[.,15]/100;
/* price index of investment expenditures */
pinvest=(data[.,7]./data[.,4]);

 bt=data[.,22];                        @ population in 1000     @
  yt=data[.,2]-(data[.,23]./pbip);     @ real GDP at factor costs in billions of mark @
 ytk=(yt./bt)*1e6;                     @ ... per capita @
  it=data[.,4];                        @ real investment expenditures in billions of mark  @
 itk=(data[.,4]./bt)*1e6;              @ ... per capita  @
  ct=data[.,3];                        @ real private consumption expenditiures in billions of mark @
 ctk=(data[.,3]./bt)*1e6;              @ ... per capita @
 ytn=(data[.,2].*pbip)-data[.,23];     @ nominal GDP at factor costs in billions of mark @
  ht=data[.,17];                       @ working hours in millions  @
 htk=(data[.,17]./bt)*1000;            @ ... per capita @
  kt=GetKapital;                       @ Compute a measure of quartely capital input  @
 ktk=(kt./bt)*1e6;                     @ ... per capita @
 ntk=(htk/(16*90));                    @ working hours as a fraction of  1440 hours per quarter @
 wt =(data[.,16]./data[.,36])*1000;    @ hourly compensation of workers in marks   @
 wt = (100*wt)./data[.,37];            @ real hourly compensation of workers @

@ output file=capital.txt reset;  kt; output off; @


/* nstar */
 n=meanc(ntk);
 
 yk=meanc(yt./kt);
 
/* Compute trend growth rate of output */
{trend1,resid1,b1}=LTREND(ln(ytk));  @ Reales Bip zu Faktorkosten pro Kopf und Quartal @
a=exp(b1[1,2]);

/* compute delta */
dreal=data[.,6]./pinvest;
delta=meanc(dreal./kt);

/* compute alpha */
alpha=GetAlpha;

/* different measursures of the real rate of return:
** King, Plosser and Rebelo: return on equity
*/
Dax=data[.,34]./pbip; @ Dax index deflated by price index of GDP @
rstar1=(Dax[rows(dax)]/Dax[1])^(1/(rows(dax)-1)); @ average quarterly growth rate of DAX @
rstar1=100*(rstar1-1);
Faz=data[.,35]./pbip; @ Faz index deflated by price index of GDP @
rstar2=(Faz[rows(faz)]/Faz[1])^(1/(rows(faz)-1));
rstar2=100*(rstar2-1);

/* ex post real rate of interest */
dp=ln(data[2:rows(data),15])-ln(data[1:rows(data)-1,15]);     @ quarterly rate of inflation @
z1=(1+(data[.,29]/100))^(1/4) - 1;                            @ nominal rate of interest on money market loans (three quarters)  @
rz1=z1[1:rows(z1)-1]-dp;
rstar3=100*meanc(rz1);

/* Solow-Residuum:
**
**  z=y/((a*h)^alpha)*(k^(1-alpha)).
**
**  where a is computed from the average quarterly growth rate of output per capita
**
*/

at=zeros(rows(yt),1);at[1]=1;
i=1;
do until i>rows(at)-1;
  at[i+1]=at[i]*a;
  i=i+1;
endo;

zt=yt./((((at.*ht))^alpha).*(kt^(1-alpha)));

zhut=(zt-meanc(zt))./meanc(zt);

{rho,stderr,twert,pwert,resid,rbar}=MyOls(trimr(zhut,1,0),trimr(zhut,0,1),"");

sigma=(resid'*resid)/(rows(resid)-1);
sigma=sqrt(sigma);

/* Momente berechnen */
yt=HPFilter(ln(ytk),1600);
it=HPFilter(ln(itk),1600);
ct=HPFilter(ln(ctk),1600);
ht=HPFilter(ln(htk),1600);
wt=HPFilter(ln(wt),1600);
x=(yt~it~ct~ht~wt);
vc=vcx(x);
cr=corrx(trimr(x,1,0)~trimr(x,0,1));

/* Daten fuer die Ausgabe zusammenstellen */
out="Output  "~"" |
    "Investme"~"nt" |
    "Consumpt"~"ion"|
    "Hours   "~""   |
    "Real Wag"~"e";

out=out~sqrt(diag(vc))*100;  @ Standard devistions@
out=out~cr[1:5,1];           @ Crosscorrelation with Ouput @
out=out~diag(cr[1:5,6:10]);  @ First order autocorrelation @

mask=(0~0~1~1~1);
fmt={"-*.*s" 8 8, "-*.*s" 8 8, "*.*lf" 8 2, "*.*lf" 8 2,"*.*lf" 8 2};

output file=Moments.txt reset;
i=printfm(out,mask,fmt);
?"Column 1: Variable name";
?"Column 2: Standard deviation";
?"Column 3: Cross correlation with output";
?"Column 4: Firt order autocorrelation";
output off;

output file=Parameters.txt reset;

?"Results:";
?"a     = " ftos(a,"*.*lf",6,4);
?"Alpha = " ftos(alpha,"*.*lf",6,2);
?"Delta = " ftos(delta,"*.*lf",6,3);
?"n     = " ftos(n,"*.*lf",6,2);
?"Rstar1= " ftos(rstar1,"*.*lf",8,4);
?"Rstar2= " ftos(rstar2,"*.*lf",8,4);
?"Rstar3= " ftos(rstar3,"*.*lf",8,4);
?"yk    = " ftos(yk,"*.*lf",8,4);
?"Rho   = " ftos(rho,"*.*lf",6,3);
?"Sigma = " ftos(Sigma,"*.*lf",8,4);

output off;
ende:
end;

/* -------------------------------- Procedures ------------------------- */

/* GetAlpha: Wage Share
**
** Lohnquote: Das reale Bruttolohneinkommen wird ergaenzt um ein hypothetisches
**            Lohneinkommen fuer die Selbstaendigen in Hoehe des durchschnittlichen
**            Lohneinkommens der Unselbstaendigen. Dieses Lohneinkommen wird auf
**            das Bruttoinlandsprodukt zu Faktorkosten bezogen.
*/

proc(1)=GetAlpha;

   local l, alpha;
       
       l=100*(data[.,16]./data[.,15]);
       l=l+(l./data[.,18]).*data[.,19]; 

  alpha=meanc(l./yt);

 retp(alpha);
endp;

/* GetKapital: Measure of the quarterly capital stock */
proc(1)=GetKapital;

    local data2, preisindex, dreal, irealnetto, kj, kq, kt, i, t, a;

    load data2[31,7]  =NettoAV.txt;      @ annual capital stock data          @
    if _short;
       data2=data2[16:31,.];             @ shorten to 1975-1989        @
    endif;

    preisindex=data[.,7]./data[.,4];     @ price index investment expenditures @
    dreal=data[.,6]./preisindex;         @ reale depreciation          @
    irealnetto=data[.,4]-dreal;          

    kj=data2[.,2]+data2[.,3]+data2[.,4]+data2[.,5]+data2[.,7];

   /* Berechnung des Kapitalstocks:
   **
   ** Da die Summe der Nettoinvestitionen eines Jahres in der Regel nicht
   ** der Differenz zweier benachbarter Jahreswerte des Kapitals entspricht,
   ** werden die Investitionen gleichmaessig nach oben oder unter skaliert:
   **
   ** a=(kj[t+1]-kj[t])/sumc(I[t,1:4])
   **
   ** Die ausgewiesenen Kapitalbestaende sind Jahresanfangswerte.
   */

   kq=zeros(rows(irealnetto)+1,1);kt=ones(rows(kq),1).*miss(1,1);
   i=1;t=1;
   do until i>rows(kj)-1;
 
            a=(kj[i+1]-kj[i])/sumc(irealnetto[4*(i-1)+1:4*(i-1)+4]);
       kq[t]=kj[i];
     kq[t+1]=kq[t]+a*irealnetto[4*(i-1)+1];
     kq[t+2]=kq[t+1]+a*irealnetto[4*(i-1)+2];kt[t]=kj[i];
     kq[t+3]=kq[t+2]+a*irealnetto[4*(i-1)+3];
     kq[t+4]=kq[t+3]+a*irealnetto[4*(i-1)+4];

           i=i+1;
           t=t+4;  
   endo;
   kq=kq[1:rows(kq)-1];
   kt=kt[1:rows(kt)-1];
   GraphSettings;
   title("Quarterly and annual (squares) stock of capital");
   _plctrl=(0|-1|0);
   xy(seqa(76,0.25,rows(kq)),kq~kt);

retp(kq);
endp;

@ ---------------------------------- LTREND -------------------------------- @
@                                                                            @
@ Erstellt am 20.09.1994 von Alfred Maussner                                 @
@                                                                            @
@ Zweck:  Schaetzt einen linearen Trend                                      @
@                                                                            @
@ Aufruf: {trend,resid,b}=LTREND(y);                                         @
@                                                                            @
@ y      n mal 1 Vektor mit der Zeitreihe, fuer die der Trend geschaetzt     @
@        werden soll.                                                        @
@ trend  n mal 1 Vektor mit dem geschaetzten linearen Trend                  @
@                                                                            @
@ resid  n mal 1 Vektor mit der Trendabweichung, resid=y-trend               @
@                                                                            @
@ b      2 mal 2 Matrix, deren erste Zeile die geschaetzten Koeffizienten    @
@        der Trendgeraden enthaelt und deren zweite Zeile die zugehoerigen   @
@        t-Werte                                                             @
@                                                                            @

proc (3)=ltrend(y);
 
 local x,zeit,b,stderr,twert,pwert,resid,rbar,trend;


 /* Aufstellen der Regressorenmatrix */
 zeit=seqa(0,1,rows(y));
 x=ones(rows(y),1)~zeit;

 /* OLS-Schaetzung der Gleichung x=a1 + a2*t */
 output file=TREND.OUT reset;screen off;
 {b,stderr,twert,pwert,resid,rbar}=MYOLS(y,x,"KONST"|"STEIGUNG");
 output off; screen on;
 trend=b[1]+b[2]*zeit;

 retp(trend,resid,(b'|twert'));

endp;

@ ----------------------------- MYOLS -------------------------------------- @
@                                                                            @
@ Erstellt am 26.05.1994 von Alfred Maussner                                 @
@ Zuletzt geandert am 13.04.2007                                               @
@                                                                            @
@ Zweck: Einfache OLS-Schaetzung (es wird unterstellt, dass es keine         @
@        fehlenden Daten gibt)                                               @
@                                                                            @
@ Aufruf: {b,stderr,twert,pwert,resid,rbar}=MYOLS(y,x,var)                   @
@                                                                            @
@ b      := k mal 1 Vektor der geschaetzten Koeffizienten                    @
@ stderr := k mal 1 Vektor der geschaetzten Standardfehlter                  @
@ twert  := k mal 1 Vektor der t-Werte                                       @
@ pwert  := k mal 1 Vektor der p-Werte                                       @
@ resid  := n mal 1 Vektor der Residuen                                      @
@ y      := n mal 1 Vektor mit den n Beobachtungen fuer die unabhaengige     @
@                   Variable                                                 @
@ x      := n mal k Matrix mit den n Beobachtungen fuer die k unabhaengigen  @
@                   Variablen (inklusive einer Konstante, falls gewuenscht)  @
@ var    := k mal 1 Vektor mit den Namen der unabhaengigen Variablen.        @
@                   Falls Skalar 0 uebergeben wird, wird die Ausgabe der Er- @
@                   gebnisse unterdrueckt                                    @
@                                                                            @

proc (6)=MYOLS(y,x,var);

  local m,b,stderr,twert,rbar,dwstat,resid,sig,n,k,out,fmt,mask,i,fwert,probf,
        rsq,pwert,vc;


  /* Ueberpruefung auf Missing Values */

  if ismiss(y~x);
   ?"Die Daten enthalten Missing Values! Weiter mit beliebiger Taste!";
   wait;
   retp(-1,-1,-1,-1,-1,-1);
  endif;

  @ Anzahl der Regressoren und Anzahl der Beobachtungen @

  k=cols(x);
  n=rows(x);

  @ Eingabe pruefen @

  if rows(y) /= n;
   ?"y und x haben nicht dieselbe Zahl von Beobachtungen!";
   ?"Weiter mit beliebiger Taste!";
   wait;
   retp(-1,-1,-1,-1,-1,-1,-1,-1);
  endif;

  if n<k;
   ?"Es gibt zuwenige Beobachtungen um die Gleichung zu schaetzen!";
   ?"Weiter mit beliebiger Taste!";
   wait;
   retp(-1,-1,-1,-1,-1,-1,-1,-1);
  endif;

  /* Schaetzungen */

  m=invpd(moment(x,0));                        @ OLS-Schaetzer fuer beta @
  b=m*x'*y;
  resid=(y-x*b);                               @ Residuenvektor          @
  
  sig=sumc(resid^2)/(n-k);
  vc=sig*m;
  stderr=sqrt(diag(vc));                       @ Standardfehler          @

  /* R^2, Rbar^2 und der F-Test werden nur bei einer Regression mit einer
  ** Konstanten berechnet */

  if x[.,1] == ones(rows(x),1);
    rsq=1-(n-k)*sig/sumc( (y-meanc(y))^2);     @ R^2                     @rsq=real(rsq);
    rbar=1-(n-1)*(1-rsq)/(n-k);                @ Rbar^2                  @rbar=real(rbar);
    fwert=(n-k)*rsq/( (k-1)*(1-rsq) );         @ F-Test                  @fwert=real(fwert);
    probf=cdffc(fwert,k-1,n-k);
  else;
    rbar=miss(1,1);
  endif;
  twert=b./stderr;                             @ t- und p-Werte          @
  pwert=2*cdftc(abs(twert),n-k);


  /* Ausgabe der Ergebnisse */

  if strlen(var[1]) > 0;
     cls;
     ?"Ergebnisse der OLS-Schaetzung";
     ?"";
     ?"Variable         b            SE-b         t-Wert          p-Wert";
     ?"=======================================================================";
     ?"";
     fmt={"-*.*s" 8 8,"%*.*lf" 14 4,"%*.*lf" 14 4, "%*.*lf" 14 4,"%*.*lf" 14 4};
     out=var~b~stderr~twert~pwert;
     out=var~real(b)~real(stderr)~real(twert)~real(pwert);
     mask={0 1 1 1 1};
     i=printfm(out,mask,fmt);
     ?"";
     if not ismiss(rbar);
       ?"R^2 = " ftos(rsq,"%*.*lf",6,4) ", Rbar^2 = " ftos(rbar,"%*.*lf",6,4);
       ?"F(" ftos(k-1,"%*.*lf",4,0) ", " ftos(n-k,"%*.*lf",4,0) ") = "
         ftos(fwert,"%*.*lf",6,4) ", Prob = " ftos(probf,"%*.*lf",6,4);
     endif;
     ?"";?"Covariance matrix is stored in vc.fmt in the current directory";
     save vc;
  endif;

  /*  Uebergabe der Ergebnisse */
  retp(b,stderr,twert,pwert,resid,rbar);

endp;

/* HPFilter: returns the cyclical component of a time series x
**           as defined by the Hodrick-Prescott filter.
**
** usage: c=HPFilter(x,mu);
**
** Input:   x : nobs x 1, vector with observations
**
**         mu : scalar, the filter weight (i.e. mu=1600 for quarterly data)
**
** Output   c : nobs x 1, vector, the cyclical component of x
*/

proc(1)=HPFilter(x,mu);

local m, i, nobs;

  nobs=rows(x);
  
  m=zeros(nobs,3);

  /* Belegen der Matrix */
  i=1;
  do until i>nobs;
   m[i,1]=mu;
   m[i,2]=-4*mu;
   m[i,3]=1+6*mu;
        i=i+1;
  endo;
  m[1,1]=0;         m[2,1]=0;
  m[1,2]=0;         m[2,2]=-2*mu;     m[nobs,2]=m[2,2];
  m[1,3]=1+mu;      m[2,3]=1+5*mu;
  m[nobs,3]=m[1,3]; m[nobs-1,3]=m[2,3];
 
 retp(x-bandsolpd(x,m));  @ The command bandlospd provides the solution of equation (A.4.2) @

endp;   


