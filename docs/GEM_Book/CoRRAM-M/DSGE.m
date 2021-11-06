classdef DSGE
    % Definition of the canonical dynamic stochastic general equilibrium model
    %{
        Copyright: Alfred Mauﬂner
        Revisions: 10 January 2017, first version
                   11 May     2017, non zero means added
        This class summarizes both the properties of the model as well as
        the options to solve and simulate the model.
    %}
    
    properties
        nx;             % number of endogenous state variables
        ny;             % number of not predetermined variables
        nz;             % number of pure shocks
        nu;             % number of static equations
        Rho;            % autocorrelation matrix of shocks        
        Omega;          % standard deviations of innovations
        Skew;           % nz by nz^2 skewness matrix
        Muss;           % nz by one vector of second order derivatives of mean vector of innovations
        log=false;      % if true, the model is formulated in the natural logarithm of the variables
        Var=struct('Name','','Symbol','','Type','','Pos',0,'Print',0,'Corr',0,'Rel',0,'Star',0,'Xi',0,'Plotno',0,'Bound',[0 0]);    
        Equations;      % string, name the file which defines the model's equations
        order=2;        % order of solution, 1,2, or 3
        reduced=true;   % if false the QZ factorization is used
        etol=1e-9;      % tolerance used to check whether the model's equations hold at the stationary point
        itol=1e-10;     % tolerance used to check for complex solutions
        inobs=10;       % number of periods for impulse responses
        nobs=80;        % number of periods in simulations of the model
        nofs=500;       % number of simulations
        hp=1600;        % Hodrick-Prescott filter weight
        ls=true;        % if true, random numbers will be loaded from file
        lm=false;       % if true, symbolic matrices will be loaded
        syl=true;       % solve the Sylvester equation instead of using the vec-operator
        table='B';      % 
        outfile='Model';% name of the file that stores the symbolic definition of the model's equations
        loadpath=cd;    % path from which Matlab looks for the model
        numeric=false;  % if true derivatives will be approximated numerically
        lfh;            % handle of the logfile
        Plot=false;     % if true, impulse respones will be plotted according to the settings in the Var structure
        Print=false;    % if true, simulation results will be printed to file according to the settings in the Var structure
        CheckBounds=false; % if true, violation of bounds during simulation will be recorded
        Trendline=false; % if true, plots of impulse response display a trend line
        LegendBox=true; % if false, no lengend will be displayed in plots
        Grid=true;      % if false, no grid will be shown
        ds=false;       % if true, the model displays difference stationary growth
        Messages=cellstr(char('Model does not exist',...
                              'Invalid stationary solution', ...
                              'Illconditioned Jacobian', ...
                              'Not able to reduce model',...
                              'Not able to solve the reduced model', ...
                              'Schur failed', ...
                              'Instable model', ...
                              'Indetermined model'));

    end
    
    methods
        function obj = set.nx(obj,nx)
            if ~isnumeric(nx); error('nx must be non-negative integer value'); else obj.nx=nx; end;
        end
        function obj = set.ny(obj,ny)
            if ~isnumeric(ny); error('ny must be non-negative integer value'); else obj.ny=ny; end;
        end;
        function obj = set.nz(obj,nz)
            if ~isnumeric(nz); error('nz must be non-negative integer value'); else obj.nz=nz; end;
        end;
        function obj = set.nu(obj,nu)
            if ~isnumeric(nu); error('nu must be non-negative integer value'); else obj.nu=nu; end;
        end;
        function y=PFY(M,s); if s=='y'; y=M.ny; end; end;
        
        function M=DSGE(nx,ny,nz,nu,w,varargin)
            p=inputParser;
            p.addRequired('nx',@isnumeric);
            p.addRequired('ny',@isnumeric);
            p.addRequired('nz',@isnumeric);
            p.addRequired('nu',@isnumeric);
            p.addRequired('w', @isnumeric);
            p.addOptional('Symbols','',@iscell);
            p.addParameter('Names','',@iscell);
            p.parse(nx,ny,nz,nu,w,varargin{:});
            Symbols=p.Results.Symbols;
            Names=p.Results.Names;
                if nargin<5;
                    error('Arguments nx, ny, nz, nu and w are required');
                elseif nargin>=5;
                    M.nx=nx;
                    M.ny=ny;
                    M.nz=nz;
                    M.nu=nu;
                    if (nu>0 && nu<ny); M.reduced=true; else M.reduced=false; end;
                    if nx>0;
                        for i=1:nx;
                            if ~isempty(Names);
                                M.Var(i).Name=char(Names(i));
                            else
                                M.Var(i).Name=strcat('X',num2str(i));
                            end;
                            if ~isempty(Symbols);
                                M.Var(i).Symbol=char(Symbols(i));
                            else
                                M.Var(i).Symbol=strcat('X',num2str(i));
                            end;
                            M.Var(i).Type='x';
                            M.Var(i).Pos=i;
                            M.Var(i).Star=w(i);
                        end;
                    end;
                    if nz>0;
                        M.Rho=zeros(nz,nz);
                        M.Omega=zeros(nz,nz);
                        M.Skew=zeros(nz,nz*nz);
                        M.Muss=zeros(nz,1);
                        for i=1:nz;
                            if ~isempty(Names);
                                M.Var(i+nx).Name=char(Names(i+nx));
                            else
                                M.Var(i+nx).Name=strcat('Z',num2str(i));
                            end;
                            if ~isempty(Symbols);
                                M.Var(i+nx).Symbol=char(Symbols(i+nx));
                            else
                                M.Var(i+nx).Symbol=strcat('Z',num2str(i));
                            end;
                            M.Var(i+nx).Type='z';
                            M.Var(i+nx).Pos=i;
                            M.Var(i+nx).Star=w(i+nx);
                            M.Var(i+nx).Plotno=1;
                        end;
                    end;
                    if ny>0;
                        for i=1:ny;
                            if ~isempty(Names);
                                M.Var(i+nx+nz).Name=char(Names(i+nx+nz));
                            else
                                M.Var(i+nx+nz).Name=strcat('Y',num2str(i));
                            end;
                            if ~isempty(Symbols);
                                M.Var(i+nx+nz).Symbol=char(Symbols(i+nx+nz));
                            else
                                M.Var(i+nx+nz).Symbol=strcat('Y',num2str(i));
                            end;
                            M.Var(i+nx+nz).Type='y';
                            M.Var(i+nx+nz).Pos=i;
                            M.Var(i+nx+nz).Star=w(i+nx+nz);
                        end;
                    end;
                    na=nx+nz+ny;
                    for i=1:na;
                        M.Var(i).Plotno=0;
                        M.Var(i).Print=0;
                        M.Var(i).Corr=0;
                        M.Var(i).Rel=0;
                        M.Var(i).Xi=0;
                        M.Var(i).Bound(1)=0; M.Var(i).Bound(2)=0;
                    end;
                end;
        end;

    end
    
end


