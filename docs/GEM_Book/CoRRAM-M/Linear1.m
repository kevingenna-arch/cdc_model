function [Hx_w,Hy_w,rc] = Linear1(Model,D)
% Obtain linear part of perturbation solution of DSGE models via Schur decomposition

%{
    Copyright:  Alfred Maußner
    Date:     12 March 2015
    Revision: 17 January 2017

    Purpose:  obtain the coefficients of the linear solution of a DSGE model via the Schur decomposition

    Input:  Model, an instance of the DSGE class
            D, (nx+ny) by 2*(nx+nz+ny) matrix, the Jacobian of the model

    Output:  Hx_w, nx by (nx+nz) matrix, the linear part of the solution for the model's states
             Hy_w, ny by (nx+nz) matrix, the linear part of the solution for not predetermined variables
             rc,   integer, the index of the error messages, see the code in SolveModel for these messages

    Remark:  the program assumes the canonical model Eg(x(t+1),z(t+1),y(t+1),x(t),z(t),y(t))=0
             z(t+1)=Rho*z(t)+Omega*epsilon(t+1)

%}

% initialize
i_tol=Model.itol;
lfh=Model.lfh;
nx=Model.nx;
nz=Model.nz;
ny=Model.ny;
nu=Model.nu;
Rho=Model.Rho;

warning('off','MATLAB:singularMatrix');
lastwarn('');
Hx_w=0;
Hy_w=0;
rc=0;
Hw_w=1;
Hy_w=1;
nw=nx+nz;
nxy=nx+ny;
na=nx+ny+nz;
nv=ny-nu;

% Set up the linearized model
Cu=D(1:nu,1+na+nw:na+nw+nu);
Cwv=[-D(1:nu,1+na:na+nw), - D(1:nu,1+na+nw+nu:2*na)];
Dwv=[D(1+nu:nxy,1:nw), D(1+nu:nxy,1+nw+nu:na)];
Fwv=[D(1+nu:nxy,1+na:na+nw), D(1+nu:nxy,1+na+nw+nu:2*na)];
Du=-D(1+nu:nxy,1+nw:nw+nu);
Fu=-D(1+nu:nxy,1+na+nw:na+nw+nu);
CuiCwv=linsolve(Cu,Cwv);
if  MatrixSingular;
    fprintf(lfh,'Cu matrix is singular\n');
    save('Cu.mat','Cu');
    rc=4;
    return;
end;
Amat=[Dwv-Du*CuiCwv;[zeros(nz,nx), eye(nz), zeros(nz,nv)]];
Bmat=[-Fwv+Fu*CuiCwv;[zeros(nz,nx), Rho, zeros(nz,nv)]];
Wmat=linsolve(Amat,Bmat); 
if MatrixSingular;
    fprintf(lfh,'A matrix is singular\n');
    rc=5;
    return;
end;	
% obtain the complex, ordered Schur decomposition of Wmat
[Tmat,Smat]=schur(Wmat,'complex');
[Tmat,Smat]=ordschur(Tmat,Smat,'udi');
lambda=ordeig(Smat);
out=[1:length(Wmat);lambda'];
fprintf(lfh,'Eigenvalue no %3i = %10.6f\n',out);
% test for stability
test1=sum(abs(lambda(1:nx+nz))>1)>0;
if test1;
    fprintf(lfh,'Not all of the %4i eigenvalues are within the unit circle',nx+nz);
    rc=7;
    return;
end;
% test for determinancy
test1=sum(abs(lambda(nx+nz+1:nx+nz+nv))<1)>0;
if test1;
    fprintf(lfh,'Not all of the %4i eigenvalues are outside the unit circle',nv);
    rc=8;
    return;
end;
Hw_w=linsolve(Tmat(1:nx+nz,1:nx+nz)',(Tmat(1:nx+nz,1:nx+nz)*Smat(1:nx+nz,1:nx+nz))');
Lvw=linsolve(Tmat(1:nx+nz,1:nx+nz)',Tmat(1+nx+nz:nx+nz+nv,1:nx+nz)');	
Hw_w=Hw_w';
Lvw=Lvw';
test1=max(abs(imag(Hw_w)));
if test1>i_tol;
    fprintf(lfh,'There are complex coefficients ind Hw_w: %10.6f',itest);
end;
Hw_w=real(Hw_w);
test1=max(abs(imag(Lvw)));
if test1>i_tol;
    fprintf(lfh,'There are complex coefficients ind Lvw: %10.6f',itest);
end;
Lvw=real(Lvw);
Lvx=Lvw(1:nv,1:nx);
Lvz=Lvw(1:nv,nx+1:nx+nz);
Luw=CuiCwv*[eye(nx+nz);[Lvx,Lvz]];
Lux=Luw(:,1:nx);
Luz=Luw(:,nx+1:nx+nz);
Lyx=[Lux;Lvx];
Lyz=[Luz;Lvz];
Hy_w=[Lyx, Lyz];
Hx_w=Hw_w(1:nx,:);

    function rc=MatrixSingular()
        rc=0;
        [~,msgid]=lastwarn;
        if ~isempty(msgid);
            if strcmp(msgid,'MATLAB:singularMatrix');
                rc=true;
            end;
        end;
    end
end

