%{
    Copyright: Alfred Maußner

    First version: 11 January 2017 (developed from SolveModelAh)

    The script computes the Jacobian D, Hessian H, and matrix of third-order derivatives T.
    Before this is done it checks whether the system of equations holds at the stationary point.
%}

% Compute D and H via numeric differentiation
if Model.numeric;
    % Collect stationary solution in vector w
    w=zeros(na,1);
    for i=1:na; w(i)=Model.Var(i).Star; end;
    % Check for proper solution
    eq_n=0;
    Sys2=@(x)Model.Equations(x,eq_n);
    test1=Sys2([w;w]);
    if max(abs(test1))>Model.etol; % check stationary solution
        rc=2;
        return;
    end;
    % compute the jacobian matrix of the model
    [D,err]=jacobianest(Sys2,[w;w]);
    % Check for invalid Jacobian
    if all(all(isnan(D))) || all(all(isinf(D)));
        rc=3;
        return;
    end;
    if Model.order==1; return; end;    
    % compute Hessian
    H=zeros((nx+ny)*2*na,2*na);
    for i=1:1:nxy;    
        eq_n=i;
        Sys2=@(x)Model.Equations(x,eq_n);    
        [temp,err]=hessian(Sys2,[w;w]);
        H(1+(i-1)*2*na:i*2*na,:)=temp;    
    end;
else
    n_m=2*na;  % number of arguments

    % check whether the file exists
    ok=exist(Model.Equations,'file');
    if ok~=2; rc=1; return; end;
    Sys=str2func(Model.Equations);
    [f, xn]=Sys();

    if ~Model.lm;    
        DS=jacobian(f,xn);
        DS=simplify(DS);
        save(Model.outfile,'DS');
    else
        load(Model.outfile,'DS');
    end;
    if Model.order>1;
        if ~Model.lm;
            HS1=jacobian(DS(:),xn);
            HS1=simplify(HS1);
            imat=zeros(n_m,nxy);
            for i=1:n_m;
                imat(i,:)=[(i-1)*nxy+1:i*nxy];
            end;
            HS=HS1(imat(:),:);
            save(Model.outfile,'HS','-append');
        else
            load(Model.outfile,'HS');
        end;
    else
        HS=0;
    end;
    if Model.order>2;
        if ~Model.lm;
            TS=jacobian(HS1(:),xn);
            TS=simplify(TS);
            m2=n_m^2;
            imat=zeros(m2,nxy);
            %clear imat;
            for i=1:m2;
                imat(i,:)=[(i-1)*nxy+1:i*nxy];
            end;
            TS=TS(imat(:),:);
            save(Model.outfile,'TS','-append');
        else
            load(Model.outfile,'TS');
        end;    
    end;     

%{
 compute the numeric Jacobian matrix, the numeric Hessian matrix, and the numeric matrix
 of third-order derivatives. In this way the base workspace will not be cluttered with
 variables which serve just one purpose,  namely to evaluate these matrices.
 %}
    if present;
        [~,nc]=size(Par);
        for i=1:1:nc;
            assignin('caller',Par(i).Symbol,Par(i).Value);
        end;
    end;
    [~,nc]=size(Model.Var);
    for i=1:1:nc;
        str1=strcat(Model.Var(i).Symbol,'1');
        str2=strcat(Model.Var(i).Symbol,'2');
        assignin('caller',[Model.Var(i).Symbol '1'],Model.Var(i).Star);
        %assignin('caller',str1,Model.Var(i).Star);
        assignin('caller',[Model.Var(i).Symbol '2'],Model.Var(i).Star);
        %assignin('caller',str2,Model.Var(i).Star);
    end;

    test=eval(f);
    if max(abs(test))>Model.etol;
        fprintf(lfh,'Invalid stationary solution:\n');
        [nt,tmp]=size(test);
        out=[(1:nt);test'];
        fprintf(lfh,'Equation no %3u = %16.8e\n',out);
        rc=2;
        return;
    end;

    D=eval(DS);
    save(Model.outfile,'D','-append');
    if Model.order>1;
        H=eval(HS); save(Model.outfile,'H','-append'); else H=0; end;
    if Model.order>2;    
        T=eval(TS);save(Model.outfile,'T','-append');
    else
        T=0;
    end;
end;
