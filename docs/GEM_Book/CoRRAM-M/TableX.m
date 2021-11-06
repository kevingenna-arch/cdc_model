function TableX(sx,rx,valid,Model)
%Write simulation results to Excel table

%{
    Copyright Alfred Maußner

    First version: 30 Januare 2017 (from previous versions of TableX and PrintX)

    Purpose: Write second moments to Excel table. The table displays the
             the variable names in column 1, the standard deviation in column 2,
             relative standard deviations wrt to nr variables in columns 2+nr,
             crosscorrelations wrt to nc variables in columns 2+nr+nc,
             and in column 2+nr+nc+1 the first-order autocorrelation.

    Input:   . sx vector with nx+nz+ny elements, the standard deviations of the variables
             . rx 2*(nx+nz+ny) by 2*(nx+nz+ny) matrix. Rows 1:(nx+nz+ny) and columns
                  1:(nx+nz+ny) hold the the contemporaneous correlations
                  between the variables (the odering follows the ordering in the Var structure of Model),
                  columns (nx+nz+ny+1):2*(nx+nz+ny) hold the correlations between the current
                  and the one-period lagged variables.
             . Model, an instance of the DSGE structure
             
   %}

% gather information
nvar=length(Model.Var); % number of variables
np=0;                   % number of variables for which moments are requested
nr=0;                   % number of variables for which relative standard deviations are requested
ncr=0;                  % number of variables for which cross correlations are requested
maxstrl=0;              % maximal string size, required for proper formatting 

for i=1:1:nvar;
    if Model.Var(i).Print>0;
        np=np+1;
        if maxstrl<length(Model.Var(i).Name); maxstrl=length(Model.Var(i).Name); end;
        if Model.Var(i).Rel>0;  nr=nr+1; end;
        if Model.Var(i).Corr>0; ncr=ncr+1; end;
    end;
end;

% prepare table entries
maxcol=2+nr+ncr; % standard deviation, rel. standard deviations, cross correlations, first-order autocorrelation
ind=zeros(np,1);
Sx=zeros(np,maxcol); % matrix that stores second moments to be written to file
Names=cell(np,1);    % cell with names of variables to be written to file
fstring='%s %6.2f';  % format string 
ColDes=cell(maxcol+1,2); % cell with information about each column of the table
ColDes{1,1}='Column A:';
ColDes{1,2}='Variable';
ColDes{2,1}='Column B:';
ColDes{2,2}='Standard Deviation';

j=0;
for i=1:1:nvar; % names of variables in the table and standard deviations
    if Model.Var(i).Print>0;
        j=j+1;        
        ind(j)=i;    
        Sx(j,1)=sx(i);
        Names{j,1}=[Model.Var(i).Name,blanks(maxstrl-length(Model.Var(i).Name))];        
    end;
end;

string1='Standard deviation relative to variable ';
nc=1;

for i=1:1:nvar; % relative deviations   
    
    if Model.Var(i).Print>0;
        if Model.Var(i).Rel>0;
            nc=nc+1;
            Sx(:,nc)=Sx(:,1)/sx(i);
            fstring=[fstring,' %6.2f'];
            ColDes{nc+1,1}=['Column ' char(65+nc) ':'];
            ColDes{nc+1,2}=[string1,Model.Var(i).Name];            
        end;
    end;
end;

% Cross correlations
string1='Cross correlation with variable ';
for i=1:1:nvar;  
    
    if Model.Var(i).Print>0;
        if Model.Var(i).Corr>0;
            nc=nc+1;
            Sx(:,nc)=rx(ind,i);
            fstring=[fstring,' %6.2f'];
            ColDes{nc+1,1}=['Column ' char(65+nc) ':'];
            ColDes{nc+1,2}=[string1, Model.Var(i).Name];
        end;
    end;
end;

% First order autocorrelation
nc=nc+1;
tmp=diag(rx(1:nvar,nvar+1:2*nvar));
Sx(:,nc)=tmp(ind);
fstring=[fstring,' %6.2f\n'];
ColDes{nc+1,1}=['Column ',char(65+nc),':'];
ColDes{nc+1,2}='First order autocorrelation';
%Names=cell2mat(Names);

% write table to file
filen=[Model.outfile,'.xlsx'];
tab1=cell2table(Names);
tab2=array2table(Sx);
mytable=[tab1,tab2];
h=height(mytable);
nr=length(ColDes);
r=['A',num2str(h+3),':B',num2str(h+2+nr)];
dtable=table(ColDes);

writetable(mytable,filen,'WriteVariableNames',false);
writetable(dtable,filen,'Range',r,'WriteVariableNames',false);

return;

end

