% % SolveLA.m
% Burkhard Heer, 12 August 2006
%
% main parts of this code are taken from Burnside, 1999, Real Business
% Cycles: Linear Approximation and GMM estimation,
% all remaining errors are mine
%
% solves a system of linear stochastic difference equations
% as specified in equation (2.36) in Heer/Maussner
%
%
% nx -- numbers of rows from x
%
function [Lxx,Lxz,Llx,Llz,Lux,Luz] = SolveLA(Cu,Cxl,Cz,Dxl,Fxl,Du,Fu,Dz,Fz,rho,nx)


    toli=1e-8; % Tolerance in case of numerical imprecision producing complex policy matrices
    
    % get dimension of vectors x, u, lambda
    nl=size(Dxl',1)-nx;
    nz=size(Cz',1);
    
    % Computaton of W and R
    
    if (length(Cu)==1)
        if (Cu==0)
            Cui=0;
        else
            Cui=1/Cu;
        end
    else
        Cui=inv(Cu);
    end
    
    disp('Cu inverted');
    
    temp=inv(Dxl-Du*Cui*Cxl);
    disp('Dxl-Du*Cui*Cxl inverted');
    W=-temp*(Fxl-Fu*Cui*Cxl);
    disp('W computed');
    R=temp*(Dz+Du*Cui*Cz);
    Q=temp* (Fz+Fu*Cui*Cz);
    disp('R computed');
    R*rho+Q
    
    [pr,lambr]=eig(W);
    alamb=abs(diag(lambr)) ;

    [lambs,lambz]=sort(alamb) ;


    disp('lambda: ');
    lambs
    
    if (sum(lambs(1:nx)>1.0)>0)
        disp('Not all of nx eigenvalues are within the unit circle');
        disp('no saddlepoint stabilitity');
        disp('the program terminates');
        pause
        return
    end

    if (sum(lambs(nx+1:nx+nl)<1)>0)
        disp('not all of the nl eigenvalues are outside the unit circle');
        disp('sun spots');
        disp('the program terminates');
        pause
        return
    end

    lambda=lambr(lambz,lambz) ;
    p=pr(:,lambz) ;

    lamb1=lambda(1:nx,1:nx) ;
    lamb2=lambda(nx+1:nx+nl,nx+1:nx+nl) ;

    p11=p(1:nx,1:nx) ;
    p12=p(1:nx,nx+1:nx+nl) ;
    p21=p(nx+1:nx+nl,1:nx) ;
    p22=p(nx+1:nx+nl,nx+1:nx+nl) ;

    ps=inv(p) ;
    ps11=ps(1:nx,1:nx) ;
    ps12=ps(1:nx,nx+1:nx+nl) ;
    ps21=ps(nx+1:nx+nl,1:nx) ;
    ps22=ps(nx+1:nx+nl,nx+1:nx+nl) ;

    rxe=R(1:nx,1:nz) ;
    rle=R(nx+1:nx+nl,1:nz) ;
    qxe=Q(1:nx,1:nz) ;
    qle=Q(nx+1:nx+nl,1:nz) ;

    phi0=ps21*rxe+ps22*rle ;
    phi1=ps21*qxe+ps22*qle ;

    psi=zeros(nl,nz) ;

    for i=1:nl
     psi(i,:)=-(phi0(i,:)*rho+phi1(i,:))*inv(eye(nz)-rho/lamb2(i,i))/lamb2(i,i) ;
    end

    Lxx=p11*lamb1*inv(p11) ;
    Lxz=(p11*lamb1*ps12+p12*lamb2*ps22)*inv(ps22)*psi+rxe*rho+qxe ;

    Llx=-inv(ps22)*ps21 ;
    Llz=inv(ps22)*psi ;

    
    cxl=inv(Cu)*Cxl ;
    ce=inv(Cu)*Cz ;
    Lux=cxl*[ eye(nx) ; Llx ];
    Luz=cxl*[ zeros(nx,nz) ; Llz ]+ce ;


end        