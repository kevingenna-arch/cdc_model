%   Ch7_US_debt_transition.m
%   computes the steady state and the transition in Ch. 7.5 of Heer (2018)
%
%   first run Ch7_US_debt.m
%
%   author: Burkhard Heer
%   date: July 20, 2018
%
% first run case_beta_calib=1, case_growth=1, year=2015, case_pen=1 for the benchmark case
% case_productivity=1
% so that the parameter beta gets calibrated and saved


clear all     % clear variables and functions from memory
close all     % closes all figure windows
clc           % clears the command window and homes the cursor.


disp('This program computes the solution for the TRANSITION DYNAMICS');
disp(' ' );
disp('in the Chapter 7.5 of Heer, Public Economics-The Macroeconomic Perspective.');
disp(' ' );
disp('First you need to run Ch7_US_debt.m in the');
disp('same directory to save starting values!!');
disp(' ');
disp('Next you need to run this program with policy=0');
disp(' ');
disp('run time with MATLAB R2013 on my pc: 33 minutes');
disp(' ');
disp('Next you may run the other policies by setting');
disp('policy=1 or policy=2 or policy=3 in line 98');
disp('run time for policy=2 and policy=3: 4 hours');
disp(' ');
disp('Hit any key when ready.....');
pause;
tic;        % computational time

% define global variables
def_global_USdebt
def_global_USdebt_transition


fhandle_function1 = @rftr;

run = 1;        % 1 - during first run, 0 - afterwards

% loading of the demographic data from survival_probs_US.xls";
if run==1
    survivalprobs=xlsread('survival_probs_US.xls',2,'F22:AI37');
    survivalprobs=flipud(survivalprobs);
    popgrowth=xlsread('survival_probs_US.xls',2,'F41:AI43');
    timespan=linspace(1950,2095,30);
    nrate=[timespan' popgrowth'];
    save('Ch7_US_debt','survivalprobs','popgrowth','nrate','-append');
else
    load('Ch7_US_debt','survivalprobs','popgrowth','nrate');
end

% ----------------------------------- Step 1 -----------------------------------  
%
% set model type
% =============
%
year1=2010;				% 2010 or 2050: population parameters
%year1=2050;
yearinitial=year1;

movavperiods=4;			% number of moving average periods for computation of population growth and suvival probs

case_tau=0;             % 1 -- labor income tax, 2 -- capital income tax, 3 -- consumption tax adjusts to balance budget
                        % 0 -- extra taxes are transfered lump-sum
case_tauc=1;            % 0 -- consumption tax=10%, 1 -- consumption tax=5%
case_repl=1;            % 1 -- replacement ratio = benchmark, 0 -- reduction by 10 percentage points
case_retirement=1;      % 1 -- retirement at age 65, 0 -- retirement at age 70

case_productivity=0;	% 1 -- all agents have productivity equal to one, 0 -- hump-shaped age-productivity profile
case_growth=1;			% 1 -- growth, 0 -- no growth
periodlength=5;			% 5 -- 5 years, 1 -- 1 year
%periodlength=1;
case_pen=1;				% 1 -- calibration such that replacement ratio is equal to empirical one			
						% 0 -- pension contribution rate equal to empirical one
case_level_pen=0;		% 1 -- level of pension remains at 2010 level, 0 -- repl ratio remains constant
case_level_pen1=0;
                        
case_UNscen=1;          % 1 -- medium variant, 2 -- low variant, 3 -- high variant for population projection UN (2015)
case_beta_calib=1;		% 1 -- benchmark case, the calibration of beta

benchmark=1;			% 1 -- benchmark case: saves the amount of transfers 

save_results=1;         % 1 -- save results
maxit=50;               % maximum number of iterations over aggregate capital stock
phi=0.9;                % updating parameter
tol = 1e-10;            % tolerance

% computational parameters transitionrfinalyear=2095;	
finalyear=2095;	% must be in {2050,2055,..,2095} 
				% 0 -- benchmark, transfers are constant and labor income tax adjusts
policy=0;		% 1 -- transfers adjust, 
				% 2 -- during the first ndebt periods, debt adjusts. Afterwards, labor income tax rate
				% 3 -- during the first ndebt periods, debt adjusts. Afterwards, transfers
				
nt=200;			% number of transition periods starting in t=2010,..
update=1;		% 1 -- linear update
ndebt=8;		% number of periods during the transition where extra expenditures are financed by debt
%ndebt=10;		% does not converge for policy=2
tolt= 1e-6;     % tolerance for transition


if periodlength==1
	nage=75;
	Rage=46;            % first period of retirement 
	nr=nage-Rage+1;     % number of retirement years 
	nw=Rage-1;          % number of working years 	
elseif periodlength==5
	nage=15;
	Rage=10;            % first period of retirement 
	nr=nage-Rage+1;     % number of retirement years */
	nw=Rage-1;          % number of working years */	
else
	disp('wrong parameter period length');
	pause; 
end

% ----------------------------------- Step 2 -----------------------------------  
%
% reading of the efficiency-age profile and survival probabilities
% from excel data file 
% ===========

nage1=75;
age=linspace(20,94,nage1);
if run==1
    efage=xlsread('efficiency_profile.xls',1,'A1:A45');
    efage=efage/mean(efage);
    save('Ch7_US_debt','efage','-append');
else
    load('Ch7_US_debt','efage');
end
    
year0=(year1-1950)/5+1;

if year0==round(year0)
	popgrowth=nrate(year0-movavperiods+1:year0,case_UNscen+1);
	popgrowth=mean(popgrowth);
	popgrowth=popgrowth/100;        % data is expressed in percentage numbers
	sp=survivalprobs(1:nage,year0-movavperiods+1:year0);
	sp=mean(sp,2);
else
	disp('year1 must be a multiple of 5');
	pause;
end

sp0 = sp;

if periodlength==5;		% transformation of annual data to 5-year data: survival probs and efficiency age
	if case_productivity==0
		efage1=zeros(nw,1)
		for i=1:1:nw
			efage1(i)=mean(efage((i-1)*5+1:i*5))		% average productivity
        end
		efage=efage1;
	elseif case_productivity==1
		efage=ones(nw,1);
    else
		disp('wrong parameter case_productivity');
		pause;
    end
end

ef=efage;

% ----------------------------------- Step 3 -----------------------------------  
%
% calibration
% ===========
%

% Parameterization as in Trabandt/Uhlig
% final goods production 
alpha1=0.35;
delta=0.083;
rbbar=1.04;		% annual real interest rate on bonds
taun=0.28;
taunbar=taun;
taulbar=0.28;		% both taul+taup=0.28!!, see Mendoza, Razin (1994), p. 311
tauk=0.36;
taukbar=tauk;
if case_tauc==1
    tauc=0.05;
else
    tauc=0.10;
end
taucbar=tauc;
taup=0.124;
taupbar=0.124;		
if case_repl==1
    replacement_ratio=0.352;    % gross replacement ratio
else 
    replacement_ratio=0.352-0.1;
end
bybar=0.63;
%bybar=1.049;
%bybar=1;
if year1==2010
	gybar=0.18;
else
	gybar=0.239;
	gybar=0.18;
end

varphi=1;		% Frisch labor supply elasticity
varphi=0.3;
lbar=0.25;		% steady-state labor supply
eta1=2.0;		% 1/IES
kappa=21.5;

if case_growth==1
	ygrowth=1.02;		% annual growth factor
else
	ygrowth=1.00;
end

if periodlength==5      % transformation of annual paramaters to 5-year values
	delta=1-(1-delta)^5;
	rbbar=rbbar^5;
	ygrowth=ygrowth^5;
	bybar=bybar/5;
	popgrowth=(1+popgrowth)^5-1
end

debtoutputratio=bybar;
% load('Ch7_US_debt','beta1');
beta1 = 1.2740392;
beta1


% computation of cohort mass
mass=ones(nage,1);
for i=2:1:nage
	mass(i)=mass(i-1)*sp(i-1)/(1+popgrowth);
end
mass=mass/sum(mass);
massinitial=mass;

% load benchmark values: first run Ch7_US_debt.m
load('Ch7_US_debt','trbench','ybench','cbench','lbench','penbench','taxesbench','xstartss','bbench','gbench');
trbar=trbench;
ybar=ybench;
c=cbench;
labors=lbench;
pen=penbench;
taxes=taxesbench;
x0=xstartss;
bigb=bbench;
bigg=gbench;

utilityss=utility(c,labors);
utilityss


x00=x0;
if nw==10       % add guess for labor supply of the nw-old worker 
                % for the case that retirement period is one period later 
		ntemp=size(x0,1);
		x00=[x0(1:nage+9); 0.3; x0(nage+10:ntemp)];
		x0=x00;
end
	
% initialisation
taun1=taunbar;
tauk1=taukbar;


% -------------------------------------------------------------------------
%
% Step 4: Computation of the Initial Steady State in 2010
%
%-------------------------------------------------------------------------_  





% from steady-state computation
x0= [   0.017480913 
     0.052919677 
     0.098921501 
      0.14354518 
      0.18560132 
      0.22365105 
      0.25410633 
      0.27207996 
      0.27764643 
      0.23851540 
      0.19268324 
      0.14181508 
     0.088889711 
     0.039179495 
      0.31946936 
      0.33602093 
      0.34184321 
      0.33735841 
      0.33240906 
      0.32692851 
      0.31948597 
      0.30919177 
      0.29962217 
      0.11666830 
      0.24589630 
     0.023866726 
      0.10771268 
    0.0091176774];
y=ss_aggregate(x0);
disp('meanc(abs(y)): ss_aggregate'); 
mean(abs(y))
%pause;

[xagg, Fval] = fsolve(@ss_aggregate,x0);
    disp('accuracy');
    y=mean(abs(Fval));
    y

    
xagginitial=xagg;
	
kbar=xagg(nage+nw);
nbar=xagg(nage+nw+1);
taup=xagg(nage+nw+3);
	
tr=xagg(nage+nw+4);
disp('tr benchmark: ') 
tr


kbar0=kbar;
nbar0=nbar;
taup0=taup;
taul0=taunbar-taup0;
asset0=[0; xagg(1:nage-1)];
totalassets=asset0'*mass(1:nage);
kshare0=kbar0/totalassets;
wbar0=wagerate(kbar,nbar);
dbar0=interest(kbar,nbar);
labors=xagg(nage:nage+nw-1); 
labors0=labors;	
nopt1=ef.*labors(1:nw);
nnew=mass(1:nw)'*nopt1;	
mean_labor=nopt1'*mass(1:nw);
mean_labor=mean_labor/sum(mass(1:nw));
pen0=replacement_ratio*wbar0*mean_labor;
trbar0=tr;

disp('initial ss_aggregate: ');	
ss_aggregate(xagg)
%pause;
disp('mean(abs(ss))=0: '); 
mean(abs(ss_aggregate(xagg)))
%pause;
x00=xagg;
	
	
debtb=bigb;
[c, sy, pen, taxes]=ss_computec(x00);
cinitial=c;
syinitial=sy;
% computation of initial bequests	
temp=(1-sp(1:nage-1)).*mass(1:nage-1);
bequests=temp'*asset0(2:nage);
bequests=bequests*(1+(1-tauk)*(dbar0-delta));
bequestsinitial=bequests;
bequests		
    

% -------------------------------------------------------------------------
%
% Step 5: Computation of the Final Steady State in 2100 (or later)
%
%--------------------------------------------------------------------------  




year0=(finalyear-1950)/5+1;
if year0==round(year0)
	popgrowth1=nrate(year0-movavperiods+1:year0,case_UNscen+1);
	popgrowth1=mean(popgrowth1);
	popgrowth1=popgrowth1/100;		% data is expressed in percentage numbers
	sp1=survivalprobs(1:nage,year0-movavperiods+1:year0);
	sp1=mean(sp1');
else
	disp('year1 must be a multiple of 5');
	pause;
end


if periodlength==5
	popgrowth1=(1+popgrowth1)^5-1;
end

popgrowth=popgrowth1;
sp=sp1';

% computation of cohort mass: stationary population
mass=ones(nage,1);
for i=2:1:nage
	mass(i)=mass(i-1)*sp1(i-1)/(1+popgrowth1);
end
mass=mass/sum(mass);	
massfinal=mass;


% initial value new steady state
% number of working years */
if nw==10   % add guess for labor supply of the nw-old worker
	ntemp=size(x0,1);
	x00=[x0(1:nage+9); 0.1; x0(nage+10:ntemp)];
end
x0=x00;


case_tau=1;
% constant level of pensions?
if case_level_pen==1
	case_level_pen1=1;
end

% guess from the GAUSS program
% in GAUSS, the non-linear-eqs solver finds a solution
% while, for the same starting value, matlab fsolve does not
xfinal0=[ 0.012463313 
     0.046490950 
     0.093406870 
      0.13957845 
      0.18394923 
      0.22515182 
      0.25928091 
      0.28094567 
      0.29041689 
      0.24963854 
      0.20411340 
      0.15435574 
      0.10158521 
     0.048490417 
      0.33097012 
      0.35007921 
      0.35814673 
      0.35539214 
      0.35203829 
      0.34794695 
      0.34150367 
      0.33168376 
      0.32218646 
      0.13174788 
      0.22572670 
     0.023556443 
      0.18246243 
   -0.0068581798 ];

x0 = xfinal0;

y=ss_aggregate(x0);
disp('mean(abs(y)) for initial value in final steady state: '); 
mean(abs(y))


[xagg, Fval] = fsolve(@ss_aggregate,x0);
    disp('accuracy of final steady-state solution');
    y=mean(abs(Fval));
    y
% pause;


xstartfinal=xagg;

kbar=xagg(nage+nw);
nbar=xagg(nage+nw+1);
taup=xagg(nage+nw+3);

trbar1=xagg(nage+nw+4);
	
kbar1=kbar;
nbar1=nbar;
taun1=taun;
taup1=taup;
wbar1=wagerate(kbar,nbar); 
dbar1=interest(kbar,nbar);
labors=xagg(nage:nage+nw-1); 
labors1=labors;
asset1=[0; xagg(1:nage-1)];	
nopt1=ef.*labors(1:nw);
nnew=mass(1:nw)'*nopt1;	
mean_labor=nopt1'*mass(1:nw);
mean_labor=mean_labor/sum(mass(1:nw));
totalassets=asset1'*mass(1:nage);
kshare1=kbar1/totalassets;
	
if case_level_pen1==0
	pen1=replacement_ratio*wbar1*mean_labor;
elseif case_level_pen1==1
	pen1=pen0;
end

disp('initial~final steady state values: ');
	disp('K: ');
    [kbar0, kbar1]
	disp('N: ');
    [nbar0, nbar1]
	disp('tr: '); 
    [trbar, trbar1]
	disp('taup: ');
    [taup0, taup1]
	disp('taul: ');
    [taulbar-taup0, taulbar-taup1]
	disp('pen: ');
    [pen0, pen1]
%    pause;

[c, sy, pen, taxes]=ss_computec(xstartfinal);
cfinal=c;
syfinal=sy;	

%
% Figure with initial and final steady state
% if labor taxes are constant and transfers adjust to balance
% the government budget
%
    periods=linspace(20,90,15);
    periods1=linspace(20,60,9);
    figure
    subplot(2,2,1);
    plot(periods,asset0,periods,asset1);
    xlabel('Age');
    title('Wealth');
    legend('2010','2100');
    subplot(2,2,2);
    plot(periods1,labors0,periods1,labors1);
    xlabel('Age');
    title('Working hours');
    subplot(2,2,3);
    plot(periods,cinitial,periods,cfinal);
    xlabel('Age');
    title('Consumption');
    subplot(2,2,4);
    plot(periods,syinitial*100,periods,syfinal*100);
    xlabel('Age');
    title('Savings rate');
%    pause;



% -------------------------------------------------------------------------
%
% Step 6: Computation of the transition
%
%	-- projection of the path of (k,n,tau^n,tau^p,tr) over the years t=0,..,nt
%
% -------------------------------------------------------------------------

% computation of the demographic variables during the transition from yearinitial to finalyear
yearinitial0=(yearinitial-1950)/5+1;
yearfinal0=(finalyear-1950)/5+1;
nsp=yearfinal0-yearinitial0+1;
spvec=zeros(nage,nsp);
popgrowthvec=zeros(nsp,1);


iyear=yearinitial0-1;
i=0;
for iyear=yearinitial0:1:yearfinal0
	i=i+1;
	if iyear==round(iyear)
		popgrowth=nrate(iyear-movavperiods+1:iyear,case_UNscen+1);
		popgrowth=mean(popgrowth);
		popgrowth=popgrowth/100;		% data is expressed in percentage numbers
		sp=survivalprobs(1:nage,iyear-movavperiods+1:iyear);
		sp=mean(sp');
		spvec(1:nage,i)=sp;
		if periodlength==5
			popgrowth=(1+popgrowth)^5-1;
        end
		popgrowthvec(i)=popgrowth;
    end
end

massvec=zeros(nage,nt);			% mass of the generation j in period tp
mass0=massinitial;
massvec(1:nage,1)=mass0;
mass1=zeros(nage,1);
for i=2:1:nt
	i0=i;
	if i0>nsp  % constant population parameters sp and popgrowth after 2100
		i0=nsp
    end
	popgrowthm=popgrowthvec(i0);
	
	mass1(1)=mass0(1)*(1+popgrowthm);
    for j=2:1:nage
		mass1(j)=mass0(j-1)*spvec(j-1,i0);
    end
	mass1=mass1/(sum(mass1));	% normalization of mass of people to one
	mass0=mass1;
	massvec(1:nage,i)=mass1;
end

	
disp('computation of massvec and spvec complete');


% initial guesses for K_t,k^s_t, N_t,n^s_t, tau^p_t
% linear interpolation between old and new steady state
%
kst=zeros(nt,nage);
ktold=zeros(nt,1);
ktnew=zeros(nt,1);
ctnew=zeros(nt,1);
cst=zeros(nt,nage);
nst=zeros(nt,nw);
ntold=zeros(nt,1);
ntnew=zeros(nt,1);
ktold=linspace(kbar0,kbar1,nt);
ktold = ktold';
ntold=linspace(nbar0,nbar1,nt);
ntold = ntold';
tauntold=linspace(taunbar,taunbar,nt);
tauntold = tauntold';
tauntnew=zeros(nt,1);
tauptold=linspace(taup0,taup1,nt);
tauptold = tauptold';
tauptnew=zeros(nt,1);
pentold=linspace(pen0,pen1,nt);
pentold = pentold';
pentnew=zeros(nt,1);
trtold=ones(nt,1)*trbar0;
trtnew=zeros(nt,1);
btold=ones(nt,1)*bigb;
btnew=zeros(nt,1);

ksharetold=linspace(kshare0,kshare1,nt);
ksharetold=ksharetold';     % column vector
ksharetnew=zeros(nt,1);
xold=[ktold; ntold; tauntold; tauptold; pentold; ksharetold];
utilityt=zeros(nt+nage-1,3);
consumprofilet=zeros(nage,nt+nage-1);
assetprofilet=zeros(nage,nt+nage-1);
incomeprofilet=assetprofilet;
income1profilet=assetprofilet;
laborprofilet=zeros(nw,nt+nage-1);
% life-time utility in old steady state at tt=-nage
utilityt(1,1)=utility(cinitial,labors0);
utilityt(1,2)=-nage;	% last cohort that has consumption profile equal to initial steady state

% initialization of factor prices vector for individual decision
wseq=zeros(nage,1);
dseq=zeros(nage,1);
penseq=zeros(nage,1);
taunseq=zeros(nage,1);
taulseq=zeros(nage,1);
taupseq=zeros(nage,1);
trseq=zeros(nage,1);


ktold1=ktold;
ntold1=ntold;
tauntold1=tauntold;
tauptold1=tauptold;
ksharetold1=ksharetold;

if policy>=2;
	ndebt0=ndebt;
	ndebt=5;		% only converges for ndebt=5
end



q=0;
kritt=1+tolt;
while (q<=maxit) || (q>1 && kritt<tolt)
    q=q+1; 
    clc;

    disp('iteration over sequence capital stock/employment: '); 
    q
    kritt
    wt=wagerate(ktold,ntold);
	dbart=interest(ktold,ntold); 
	rbbart=1+(1-tauk)*(dbart-delta);
	taunt=tauntold;
	taupt=tauptold;
	pent=pentold;
	
    tt=nt;
    if update==1
        [ktnew,ntnew,tauntnew,tauptnew,pentnew,ksharetnew,ctnew,trtnew,btnew]=getkn();
        kritt=mean(abs(ktnew-ktold));
        ktold=phi*ktold+(1-phi)*ktnew;
        ntold=phi*ntold+(1-phi)*ntnew;
        tauntold=phi*tauntold+(1-phi)*tauntnew;
        tauptold=phi*tauptold+(1-phi)*tauptnew;
        pentold=phi*pentold+(1-phi)*pentnew;
        ksharetold=phi*ksharetold+(1-phi)*ksharetnew;
        trtold=phi*trtold+(1-phi)*trtnew;
        btold=phi*btold+(1-phi)*btnew;
    else
        disp('not implemented yet');
        pause;
    end
	
end


if policy>=2    % iteration over final debt

ndebt=0;    
while ndebt<ndebt0
	ndebt=ndebt+1;

    q=0;
    kritt=1+tolt;
    while (q<=maxit) || (q>1 && kritt<tolt)
        q=q+1; 
        clc;

        disp('iteration over sequence capital stock/employment: '); 
        q
        kritt
        wt=wagerate(ktold,ntold);
        dbart=interest(ktold,ntold); 
        rbbart=1+(1-tauk)*(dbart-delta);
        taunt=tauntold;
        taupt=tauptold;
        pent=pentold;
	
        tt=nt;
        if update==1
            [ktnew,ntnew,tauntnew,tauptnew,pentnew,ksharetnew,ctnew,trtnew,btnew]=getkn();
            kritt=mean(abs(ktnew-ktold));
            ktold=phi*ktold+(1-phi)*ktnew;
            ntold=phi*ntold+(1-phi)*ntnew;
            tauntold=phi*tauntold+(1-phi)*tauntnew;
            tauptold=phi*tauptold+(1-phi)*tauptnew;
            pentold=phi*pentold+(1-phi)*pentnew;
            ksharetold=phi*ksharetold+(1-phi)*ksharetnew;
            trtold=phi*trtold+(1-phi)*trtnew;
            btold=phi*btold+(1-phi)*btnew;
        else
            disp('not implemented yet');
            pause;
        end
	
    end
end     % loop: increasing the number of debt periods 
end     % policy >= 2 with debt increase
  
t=toc;
disp('Computation of transition complete');
disp(['Elapsed time: ',num2str(t/60),' minutes']);
pause;


if policy==1
	kt_policy1=ktold;
	nt_policy1=ntold;
	taunt_policy1=tauntold;
	taupt_policy1=tauptold;
	utilityt_policy1=utilityt;
	pent_policy1=pentold;
	ct_policy1=ctnew;
	trt_policy1=trtnew;
	bt_policy1=btnew;
    save('Ch7_US_debt_transition','kt_policy1','nt_policy1','taunt_policy1','-append');
    save('Ch7_US_debt_transition','taupt_policy1','utilityt_policy1','pent_policy1','-append');
    save('Ch7_US_debt_transition','ct_policy1','trt_policy1','bt_policy1','-append');
elseif policy==2;	
	kt_policy2=ktold;
	nt_policy2=ntold;
	taunt_policy2=tauntold;
	taupt_policy2=tauptold;
	pent_policy2=pentold;
	utilityt_policy2=utilityt;
	ct_policy2=ctnew;
	trt_policy2=trtnew;
	bt_policy2=btnew;
    save('Ch7_US_debt_transition','kt_policy2','nt_policy2','taunt_policy2','-append');
    save('Ch7_US_debt_transition','taupt_policy2','utilityt_policy2','pent_policy2','-append');
    save('Ch7_US_debt_transition','ct_policy2','trt_policy2','bt_policy2','-append');
elseif policy==3;	
	kt_policy3=ktold;
	nt_policy3=ntold;
	taunt_policy3=tauntold;
	taupt_policy3=tauptold;
	utilityt_policy3=utilityt;
	pent_policy3=pentold;
	ct_policy3=ctnew;
	trt_policy3=trtnew;
	bt_policy3=btnew;
    save('Ch7_US_debt_transition','kt_policy3','nt_policy3','taunt_policy3','-append');
    save('Ch7_US_debt_transition','taupt_policy3','utilityt_policy3','pent_policy3','-append');
    save('Ch7_US_debt_transition','ct_policy3','trt_policy3','bt_policy3','-append');
elseif policy==0;	
	kt_policy0=ktold;
	nt_policy0=ntold;
	taunt_policy0=tauntold;
	taupt_policy0=tauptold;
	utilityt_policy0=utilityt;
	pent_policy0=pentold;
	ct_policy0=ctnew;
	trt_policy0=trtnew;
	bt_policy0=btnew;
    save('Ch7_US_debt_transition','massvec','kt_policy0','nt_policy0','taunt_policy0','-append');
    save('Ch7_US_debt_transition','taupt_policy0','utilityt_policy0','pent_policy0','-append');
    save('Ch7_US_debt_transition','ct_policy0','trt_policy0','bt_policy0','-append');
end

    welfaret=utilityt(:,1)/utilityt(1,1);
    welfaret=(welfaret.^(1/(1-eta1)))-1;
    utilityt(:,3)=welfaret;

    nt1=30;
    periodst = linspace(2010,2010+5*nt1-5,nt1);
    periodst = periodst';
    
    figure
    subplot(3,2,1);
    plot(periodst,ktnew(1:nt1));
    xlabel('year');
	title('capital stock');
    subplot(3,2,2);
    plot(periodst,ntnew(1:nt1));
	title('labor');
    xlabel('year');
    subplot(3,2,3);
    plot(periodst,tauntnew(1:nt1)-tauptnew(1:nt1));
	title('\tau^n');
    xlabel('year');
    subplot(3,2,4);
    plot(periodst,tauptnew(1:nt1));
	title('\tau^p');
    xlabel('year');
    subplot(3,2,5);
    plot(periodst,trtnew(1:nt1));
	title('government transfers');
    xlabel('year');
    subplot(3,2,6);
    plot(periodst,btnew(1:nt1));
	title('government debt');
    xlabel('year');
    pause;

close all;
