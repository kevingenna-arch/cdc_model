%   Ch6_AK_value.m
%
%   VALUE FUNCTION ITERATION
%
%   Auerbach-Kotlikoff Model from Chapter 9.1, Heer/Maussner, 2009,
%   'Dynamic General Equilibrium Modeling: Computation Methods and
%   Applications'
%
%   author: Burkhard Heer
%   date: May 13, 2020
% 
%	Algorithm: See Algorithm 9.2.1 in Heer/Mauﬂner DGE Modeling
%

clear all     % clear variables and functions from memory
close all     % closes all figure windows
clc           % clears the command window and homes the cursor.


disp('This program computes the solution to the Auerbach-Kotlikoff model');
disp(' ' );
disp('in the Chapter 9.1 of Heer/Maussner, Dynamic General Equilibrium Modeling.');
disp(' ' );
disp('Hit any key when ready.....');
%pause;
tic;        % computational time

% define global variables
def_global_AK_value

VI_method='linear';           % 'linear' --- linear interpolation, 'cubic' --- cubic spline 
% Vq = interp2(,method) X,Y,V,Xq,Yq) returns interpolated values of a function of two 
% variables at specific query points using linear/cubic interpolation. 
% The results always pass through the original sampling of the function. 
% X and Y contain the coordinates of the sample points. 
% V contains the corresponding function values at each sample point. 
% Xq and Yq contain the coordinates of the query points.


%
%   STEP 1: Parameterization and calibration
%
% Parameterization
%fhandle_function1 = @test;
%fhandle_function2 = @wvalue;

beta0=0.96;         % discount factor 
r=0.045;            % initial value of the interest rate 
sigma=2;            % coefficient of relative risk aversion 
alpha0=0.36;        % production elasticity of capital 
rep=0.3;            % replacement ratio 
delta=0.1;          % rate of depreciation 
tr=20;              % retired
%tr=2;
t=40;               % working time 
%t=2;
tau=rep/(2+rep);    % income tax rate 
gam=2;              % disutility from working 

psi0=0.001;         % parameter of utility function 
phi=0.8;            % update parameter of aggregates
tol=0.001;          % percentage deviation of final solution 
tolk=0.001;         % percentage deviation of final solution for a^1
tol_golden=1e-8;    % tolerance for Golden Section Search
neg=-1e10;          % large negative value for initialization of
                    % value function
nq1=30;             % maximum number of iterations over a^60
nq=30;              %  maximum number of iterations over K,L

% 
% asset grid: individual wealth
assetmax=10;    % upper limit of asset grid on a 
na=200;          % number of grid points on assets
assetmin=0;
agrid=linspace(0,assetmax,na);   % asset grid	
agrid = agrid';     % column vector
aeps=(agrid(2)-agrid(1))/10;    % test for corner solution

% ---------------------------------
%
%   Step 2: initialization of K, L, tau
%
% --------------------------------- 

nbar=0.2;
kbar=(alpha0/(r+delta))^(1/(1-alpha0))*nbar;
kold=100;
nold=2;

% initialization of global variables 
k0=0;   % current period assets
k1=0;   % next-period assets

% ---------------------------------------------
%
% Step 3: iteration of policy function, wealth distribution,..
%
% ---------------------------------------------- 

q=0;
krit=1+tol;
while (q<30) | (krit>tol)
    q=q+1
    krit=abs((kbar-kold)/kbar);
    krit0=abs((nbar-nold)/nbar);
    
    w=wage(kbar,nbar)
    r=interest(kbar,nbar);
    pen=rep*(1-tau)*w*nbar*3/2;
    kold=kbar;
    nold=nbar;

    % retired agents' value function  
    vr=zeros(na,tr);    % value function 
    aropt=ones(na,tr);  % optimal asset 
    cropt=zeros(na,tr); % optimal consumption 

    % value function at at 60 = utility at age 60
    for i=1:1:na
        vr(i,tr)=u(agrid(i)*(1+r)+pen,0);
        cropt(i,tr)=agrid(i)*(1+r)+pen;
    end
    
    
    % workers' value function 
    vw=zeros(na,t);
    awopt=ones(na,t);
    cwopt=zeros(na,t);
    nwopt=zeros(na,t);

    % computation of the decision rules for the retired 
    for i=tr-1:-1:1
        clc;
        period = i;
        disp('q, i, K');
        disp([q,i+t,kbar]);
        m0=0;
        for ia=1:1:na  % loop over individual state space
            k0=agrid(ia);
            % finding the interval [ax,cx] that brackets the maximum
            ax=0;
            bx=-1;
            cx=-2;
            v0=neg;     % initial value of the value function
                        % of the t+i-year old with wealth a[ia]
            y=[0,m0-3];
            m=max(y);     % exploiting the concavity of the value function
                                % next-period wealth a'(a) is an increasing
                                % function of a
            while (ax>bx) | (bx>cx)
                m=m+1;
                k1=agrid(m);        % next-period asset
                v1=value1(k1);
                if v1>v0
                    if m==1
                        ax=agrid(m); 
                        bx=agrid(m);
                    else
                        bx=agrid(m);
                        ax=agrid(m-1);
                    end
               
                    v0=v1;
                    m0=m;   % monotonocity of the value function 
                else
                    cx=agrid(m);
                end     % if v1>v0

                
                if m==na 
                    ax=agid(m-1); 
                    bx=a(m); 
                    cx=a(m); 
                end    
            end     % while value function increases with a'    
            
            % maximum in [ax,bx]

            if ax==bx   % maximum at a'=assetmin?
                if value1(assetmin)>value1(aeps)
                    aropt(ia,i)=assetmin;
                else
                    y=golden(@value1,assetmin,aeps,agrid(2),tol_golden);
                    aropt(ia,i)=y;
                end
                    
            elseif bx==cx   % maximum at a'=assetmax?
                if value1(assetmax)>value1(assetmax-aeps)
                    aropt(ia,i)=assetmax;
                else
                    y=golden(@value1,agrid(na-1),agrid(na)-aeps,agrid(na),tol_golden);
                    aropt(ia,i)=y;
                end
            else        % inner solution
                
                y=golden(@value1,ax,bx,cx,tol_golden);
                aropt(ia,i) = y;
            end

            vr(ia,i)=value1(aropt(ia,i));
            cropt(ia,i)=(1+r)*agrid(ia)+pen-aropt(ia,i);
            
            
        end     % ia=1,..,na asset in value function of retired        
    end     % loop over the value function of the retired, i=tr,..,1
   
if q==1  
figure
plot(agrid,vr(:,1));
ylabel('Value function of the 41-year old');
xlabel('Individual wealth');
pause;

figure
plot(agrid,cropt(:,1));
ylabel('Consumption function of the 41-year old');
xlabel('Individual wealth');
pause;
    
figure
plot(agrid,aropt(:,1));
ylabel('Next-period asset function of the 41-year old');
xlabel('Individual wealth');
pause;
end

    % computation of the decision rules for the retired 
    for i=t:-1:1
         period = i;
        disp('q, i, K');
        disp([q,i,kbar]);
        m0=0;
        for ia=1:1:na  % loop over individual state space
            k0=agrid(ia);
            % finding the interval [ax,cx] that brackets the maximum
            ax=0;
            bx=-1;
            cx=-2;
            v0=neg;     % initial value of the value function
                        % of the t+i-year old with wealth a[ia]
            y=[0,m0-3];
            m=max(y);     % exploiting the concavity of the value function
                                % next-period wealth a'(a) is an increasing
                                % function of a
            while (ax>bx) | (bx>cx)
   
                m=m+1;
                k1=agrid(m);        % next-period asset
                v1=value2(k1);
                  if v1>v0
                    if m==1
                        ax=agrid(m); 
                        bx=agrid(m);
                    else
                        bx=agrid(m);
                        ax=agrid(m-1);
                    end
               
                    v0=v1;
                    m0=m;   % monotonocity of the value function 
                else
                    cx=agrid(m);
                end     % if v1>v0

                
                if m==na 
                    ax=agrid(m-1); 
                    bx=agrid(m); 
                    cx=agrid(m); 
                end    
                
            end     % while value function increases with a'    
            
            % maximum in [ax,bx]

            if ax==bx   % maximum at a'=assetmin? 
                if value2(assetmin)>value2(aeps)
                    awopt(ia,i)=assetmin;
                else
                    y=golden(@value2,assetmin,aeps,agrid(2),tol_golden);
                    awopt(ia,i)=y;
                end
                
            elseif bx==cx   % maximum at a'=assetmax?
                if value2(assetmax)>value2(assetmax-aeps)
                    awopt(ia,i)=assetmax;
                else
                    y=golden(@value2,agrid(na-1),assetmax-aeps,assetmax,tol_golden);
                    awopt(ia,i)=y;
                end
                    
            else        % inner solution                
                y=golden(@value2,ax,bx,cx,tol_golden);
                awopt(ia,i) = y;
            end

            
            
            k0=agrid(ia);
            k1=awopt(ia,i);
            
            n= 1/(1+gam)*(1-gam/((1-tau)*w)*(psi0+(1+r)*k0-k1));
            if n<0
                n=0;
            elseif n>1
                n=1;
            end
            
            c=(1-tau)*w*n+(1+r)*k0-k1;
            
            cwopt(ia,i)=c;
            nwopt(ia,i)=n;
            vw(ia,i)=value2(awopt(ia,i));
            
            
        end     % ia=1,..,na asset in value function of retired        
    end     % loop over the value function of the retired, i=tr,..,1
   
    
if q==1    
figure
plot(agrid,vw(:,1));
ylabel('Value function of the 41-year old');
xlabel('Individual wealth');
pause;

figure
plot(agrid,cwopt(:,1));
ylabel('Consumption function of the 41-year old');
xlabel('Individual wealth');
pause;
    
figure
plot(agrid,awopt(:,1));
ylabel('Next-period asset function of the 41-year old');
xlabel('Individual wealth');
pause;
    

figure
plot(agrid,nwopt(:,1));
ylabel('Labor supply function of the 41-year old');
xlabel('Individual wealth');
pause;
end

% Step 4: update of aggregate variables
%
% computation of the aggregate capital stock and employment nbar 
    kgen=zeros(t+tr,1);
    ngen=zeros(t,1);
    cgen=zeros(t+tr,1);
    kgen(1)=0;
    vexact=0;
    for j=1:1:t+tr-1
        
        if j<t+1    % worker 
            
            a1=awopt(:,j);
            kgen(j+1) = interp1(agrid,a1,kgen(j),VI_method);
            k1=kgen(j+1);
            k0=kgen(j);
            
            n= 1/(1+gam)*(1-gam/((1-tau)*w)*(psi0+(1+r)*k0-k1));
            if n<0
                n=0;
            elseif n>1
                n=1;
            end
            ngen(j)=n;
            c=(1-tau)*w*n+(1+r)*k0-k1;
            cgen(j)=c;
            vexact=vexact+beta0^(j-1)*u(cgen(j),ngen(j));
        else
            
            a1=aropt(:,j-t);
            kgen(j+1) = interp1(agrid,a1,kgen(j),VI_method);
           
            cgen(j)=(1+r)*kgen(j)+pen-kgen(j+1);
            vexact=vexact+beta0^(j-1)*u(cgen(j),0);
        end
    end
    
    % age 60 prior to death, no next-period wealth
    cgen(t+tr)=(1+r)*kgen(t+tr)+pen;
    vexact=vexact+beta0^(t+tr-1)*u(cgen(t+tr),0);

    knew=mean(kgen);
    kbar=phi*kold+(1-phi)*knew;
    nnew=mean(ngen)*2/3;
    nbar=phi*nold+(1-phi)*nnew;
    
    
end

save CH6_AK_value;

ttime=toc;
disp('Computation of transition complete');
disp(['Elapsed time: ',num2str(ttime/60),' minutes']);
pause;


    
figure
plot(agrid,vw(:,1));
ylabel('Value function of the 41-year old');
xlabel('Individual wealth');
pause;

figure
plot(agrid,cwopt(:,1));
ylabel('Consumption function of the 41-year old');
xlabel('Individual wealth');
pause;
    
figure
plot(agrid,awopt(:,1));
ylabel('Next-period asset function of the 41-year old');
xlabel('Individual wealth');
pause;
    

figure
plot(agrid,nwopt(:,1));
ylabel('Labor supply function of the 41-year old');
xlabel('Individual wealth');
pause;


periods=linspace(20,79,60);

figure
plot(periods,kgen);
ylabel('Individual wealth');
xlabel('Real-life age');
pause;

periods1=linspace(20,59,40);
figure
plot(periods1,ngen);
ylabel('Individual labor supply');
xlabel('Real-life age');
pause;

figure
plot(periods,cgen);
ylabel('Individual consumption');
xlabel('Real-life age');
















