function [Ew, Ey, Varw, Vary]=Moments_TD(PF,Rho,Omega,Mu)

% Computes unconditional first and second moments from the policy functions of DSGE models

%{ 
    Copyright: Alfred Mauﬂner
    Revsions:  19 September 2016, first version
               12 May       2017, adapted to CoRRAM-M
               25 May       2017, mean of innovations added

    Usage: [Ew, Ey, Varw, Vary]=Moments_TD(H, Rho, Omega, Mu)

    Input:  PF, structure with elements Hx_w, Hy_w, Hx_ww, Hy_ww, Hw_ss, Hy_ss which
                hold, respectively, the linear part of the policy functions for
                the states x and the controls y wrt the vector w (endogenous states x
                stacked over the shocks z), the quadratic part of the policy functions
                for the states and controls wrt to w, and the quadratic part of
                wrt to the perturbation parameter sigma.

            Rho, nz by nz matrix from the VAR(1) process for the shocks z(t+1)=Rho*z(t)+sigma*Omega*eps(t+1)
          Omega, nz by nz matrix (see above equation)
             Mu, nz by 1 vector, the mean of the innovations

    Output: Ew, nw by 1 vector of unconditional expected first moments of the states w=[x;z]
            Ey, ny by 1 vector of unconditional expected first moments of the controls y
          Varw, nw by nw covariance matrix of the states (unconditional second moments)
          Vary, ny by ny covariance matrix of the controls (unconditional second moments)

    Remarks: The computations employ the formulas derived in the most recent verions of
             my lecture note "Computational Macroeconomics"
    
%}
  

%get dimensions
nx=size(PF.Hx_w,1);
nz=size(Rho,1);
ny=size(PF.Hy_w,1);
nw=nx+nz;
    
% compute Hw, Gw, and sigma_u
Hw=[PF.Hx_w; zeros(nz,nx) Rho];   
Gw=PF.Hy_w;
sigma_u=[zeros(nx,nz) zeros(nx,nz); zeros(nz,nx) Omega*Omega'];
Hww=[PF.Hx_ww;zeros(nz*nw,nw)];
Gww=PF.Hy_ww;
        
% solve the discrete Lyapunov equation for Varw
Varw=Lyapunov(Hw,sigma_u);

%solve for Vary
Vary=Gw*Varw*Gw';
    
% compute the unconditional moments from the second-order approximation
amat=eye(nw)-Hw;
bvec=zeros(nw,1);
for i =1:nw
    bvec(i)=sum(diag(Hww((i-1)*nw+1:i*nw,:)*Varw));
end
bvec=0.5*bvec+0.5*[PF.Hx_ss;zeros(nz,1)]+[zeros(nx,1);Mu];
Ew=linsolve(amat,bvec);
bvec=zeros(ny,1);
for i=1:ny
    bvec(i)=sum(diag(Gww((i-1)*nw+1:i*nw,:)*Varw))+PF.Hy_ss(i);
end
Ey=Gw*Ew+0.5*bvec;
  
end