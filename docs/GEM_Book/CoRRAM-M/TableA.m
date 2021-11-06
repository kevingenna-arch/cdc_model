function TableA(sx,rx,valid,Model)
%Write simulation results to ASCII table

%{
    Copyright Alfred Maußner

    First version: 30 Januare 2017 (from previous versions of TableA and PrintA)

    Purpose: Write second moments to ASCII table. The table displays the
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

ColDes{1,1}='Column 1:';
ColDes{1,2}='Variable';
ColDes{2,1}='Column 2:';
ColDes{2,2}='Standard Deviation';
string1='Column ';

j=0;
for i=1:1:nvar; % names of variables in the table and standard deviations
    if Model.Var(i).Print>0;
        j=j+1;        
        ind(j)=i;    
        Sx(j,1)=sx(i);
        Names{j,1}=[Model.Var(i).Name,blanks(maxstrl-length(Model.Var(i).Name))];        
    end;
end;

string2='Standard deviation relative to variable ';
nc=1;

for i=1:1:nvar; % relative deviations   
    
    if Model.Var(i).Print>0;
        if Model.Var(i).Rel>0;
            nc=nc+1;
            Sx(:,nc)=Sx(:,1)/sx(i);
            fstring=[fstring,' %6.2f'];           
            string3=[string1,num2str(nc+1),':'];
            string4=[string2,Model.Var(i).Name];            
            ColDes{nc+1,1}=string3;
            ColDes{nc+1,2}=string4;
        end;
    end;
end;

% Cross correlations
string2='Cross correlation with variable ';
for i=1:1:nvar;  
    
    if Model.Var(i).Print>0;
        if Model.Var(i).Corr>0;
            nc=nc+1;
            Sx(:,nc)=rx(ind,i);
            string3=[string1,num2str(nc+1),':'];
            string4=[string2,Model.Var(i).Name];
            fstring=[fstring,' %6.2f'];
            ColDes{nc+1,1}=string3;
            ColDes{nc+1,2}=string4;
        end;
    end;
end;

% First order autocorrelation
nc=nc+1;
tmp=diag(rx(1:nvar,nvar+1:2*nvar));
Sx(:,nc)=tmp(ind);
fstring=[fstring,' %6.2f\n'];
string3=[string1,num2str(nc+1),':'];
ColDes{nc+1,1}=string3;
ColDes{nc+1,2}='First order autocorrelation';
Names=cell2mat(Names);

% write table to file
filen=strcat(Model.outfile,'.txt');
id=fopen(filen,'w');
fprintf(id,strcat(char(datetime('now')),'\n'));
fprintf(id,'%s\n','');

if Model.CheckBounds
     fprintf(id,'Second moments from %5i valid simulations. HP-Filter=%.0f\n',valid,Model.hp);
else
    fprintf(id,'Second moments from %5i simulations. HP-Filter=%.0f. Bounds were not checked.\n',Model.nofs,Model.hp);
end;

switch Model.order
    case 1
        fprintf(id,'Linear policy functions were used.\n');
    case 2
        fprintf(id,'Quadratic policy functions were used.\n');
    case 3
        fprintf(id,'Cubic policy functions were used.\n');
    otherwise;
end;
fprintf(id,'%s\n','');
for i=1:1:np;
    fprintf(id,fstring,Names(i,:),Sx(i,:));
end;
ntab=length(ColDes);
fprintf(id,'%s\n','');
for i=1:1:ntab;
    fprintf(id,'%s %s\n',ColDes{i,:});
end;
fclose(id);
return;

end

