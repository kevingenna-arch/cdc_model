function [Hx_w,Hy_w,rc] = Linear2(Model,D) 
% Obtain linear part of perturbation solution of DSGE models via QZ decomposition

%{
    Copyright:  Alfred Maußner
    Date:     10 March 2015
    Revision: 17 January 2017

    Purpose:  obtain the coefficients of the linear solution of a DSGE model via the QZ decomposition

    Input:  Model, an instance of the DSGE class
            D, (nx+ny) by 2*(nx+nz+ny) matrix, the Jacobian of the model

    Output:  Hx_w, nx by (nx+nz) matrix, the linear part of the solution for the model's states
             Hy_w, ny by (nx+nz) matrix, the linear part of the solution for not predetermined variables
             rc,   integer, the index of the error messages, see the code in SolveModel for these messages

    Remarks:  the program assumes the canonical model Eg(x(t+1),z(t+1),y(t+1),x(t),z(t),y(t))=0
              z(t+1)=Rho*z(t)+Omega*epsilon(t+1)

              The linearized model is

              B*Et(w(t+1); y(t+1))=A*(w(t); y(t)), where w(t)=(x(t); z(t))

              The solution of the model are the matrices in

              x(t+1)=Hx_w*(x(t);z(t))
              y(t)  =Hy_w*(x(t);z(t))


%}


% Initialize
i_tol=Model.itol;
lfh=Model.lfh;
nx=Model.nx;
nz=Model.nz;
ny=Model.ny;
Rho=Model.Rho;

warning('off','MATLAB:singularMatrix');
lastwarn('');

rc=0;
Hw_w=1;
Hy_w=1;
na=nx+ny+nz;

% Set up the linearized model
Amat=[-D(:,1+na:2*na);[zeros(nz,nx), Rho, zeros(nz,ny)]];
Bmat=[D(:,1:na);[zeros(nz,nx),eye(nz),zeros(nz,ny)]];
save('testbaleig','Amat','Bmat');
% obtain the complex, ordered generalized Schur decomposition of the pencil (Amat, Bmat)
[T,S,Q,Z] = qz(Amat,Bmat,'complex');
[T,S,~,Z] = ordqz(T,S,Q,Z,'udi');
lambda=ordeig(T,S);
out=[1:length(Amat);lambda'];
fprintf(lfh,'Eigenvalue no %3i = %10.6f\n',out);
% test for stability
test1=sum(abs(lambda(1:nx+nz))>1)>0;
if test1;
    fprintf(lfh,'Not all of the %4i eigenvalues are within the unit circle\n',nx+nz);
    rc=7;
    return;
end;
% test for determinancy
test1=sum(abs(lambda(nx+nz+1:nx+nz+ny))<1)>0;
if test1;
    fprintf(lfh,'Not all of the %4i eigenvalues are outside the unit circle\n',nv);
    rc=8;
    return;
end;
Temp1=linsolve(S(1:nx+nz,1:nx+nz)',Z(1:nx+nz,1:nx+nz)');
if MatrixSingular; rc=6; return; end;
Temp2=linsolve(Z(1:nx+nz,1:nx+nz)',T(1:nx+nz,1:nx+nz)');
if MatrixSingular; rc=6; return; end;
Hw_w=(Temp1')*(Temp2');
Hy_w=linsolve(Z(1:nx+nz,1:nx+nz)',Z(nx+nz+1:nx+nz+ny,1:nx+nz)');
Hy_w=Hy_w';
test1=max(abs(imag(Hw_w)));
if test1>i_tol;
    fprintf(lfh,'There are complex coefficients ind Hw_w: %10.6f\n',test1);
end;
Hw_w=real(Hw_w);
Hx_w=Hw_w(1:nx,:);

test1=max(abs(imag(Hy_w)));
if test1>i_tol;
    fprintf(lfh,'There are complex coefficients ind Hy_w: %10.6f\n',test1);
end;
Hy_w=real(Hy_w);

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

