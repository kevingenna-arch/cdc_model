function [Ex, Ey]=Firs_Moments_Sim(Model,PF,nobs)

% Computes unconditional first moments from the policy functions of DSGE models via Simulation

%{ 
    Copyright: Alfred Mauﬂner
    Revsions:  27 June 2017, first version

    Usage: [Ex, Ey]=Moments_TD(H, Rho, Omega, Mu)

    Input: Model, an instance of the DSGE object
 
            PF, structure with elements Hx_w, Hy_w, Hx_ww, Hy_ww, Hw_ss, Hy_ss which
                hold, respectively, the linear part of the policy functions for
                the states x and the controls y wrt the vector w (endogenous states x
                stacked over the shocks z), the quadratic part of the policy functions
                for the states and controls wrt to w, and the quadratic part of
                wrt to the perturbation parameter sigma.

            nobs, integer, number of simulations

    Output: Ex, nx by 1 vector of unconditional expected first moments of the states
            Ey, ny by 1 vector of unconditional expected first moments of the controls y

    
%}

% get stationary values of the states and costates
xstar=zeros(Model.nx,1);
ystar=zeros(Model.ny,1);
for i=1:Model.nx+Model.nz+Model.ny
    if strcmp(Model.Var(i).Type,'x')
        xstar(Model.Var(i).Pos)=Model.Var(i).Star;
    end
    if strcmp(Model.Var(i).Type,'y')
        ystar(Model.Var(i).Pos)=Model.Var(i).Star;
    end
end

% draw random numbers
eps=random('Normal',0,1,[nobs+1,Model.nz]);

% set starting values
xsum=zeros(Model.nx,1);
ysum=zeros(Model.ny,1);
x1=xstar;
z1=Model.Omega*eps(1,:)';
v1=[x1-xstar;z1];
xsum=x1;
% iterate
for t=1:nobs
    x2=PFX(v1);
    y1=PFY(v1);
    ysum=ysum+y1;
    xsum=xsum+x2;
    z1=Model.Rho*z1 + 0.5*Model.Muss+Model.Omega*eps(t+1,:)';
    v1=[x2-xstar;z1];
end
Ex=xsum/nobs;
Ey=ysum/nobs;

return;
 

function x=PFX(v) % the policy function
  
    x=xstar+PF.Hx_w*v;    
    if Model.order==1; return; end;
    x=x+0.5*PF.Hx_ss+0.5*kron(eye(Model.nx),v')*PF.Hx_ww*v;    
    
    return;
end

% PFY: the policy function for the costates: y=PFY(v)
%
% Input: 
%        . v:    nx+nz+1 by 1 matrix with elements [x(t)-xstar;z(t);1]
%         
% Output:
%        . y: ny by 1 matrix

function y=PFY(v) % the policy function
  
    y=ystar+PF.Hy_w*v;    
    if Model.order==1; return; end;
    y=y+0.5*PF.Hy_ss+0.5*kron(eye(Model.ny),v')*PF.Hy_ww*v;    
    
    return;
end


end