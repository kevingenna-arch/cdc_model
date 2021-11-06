function maxl = GetMaxPlotNumber(Var,iz)
% Gather information about the number of panels and the number of lines per panel

%{
 
    Copyright Alfred Maußner

    First version:  06 March 2015
    Revision:       02 February 2017 (remove panels with zero lines)

    Purpose : For each of 8 potential panels find the maximum number panels to plot to
              as well as the number of curves to be plotted to each panel. Check
              for consequtive panels with at least one line
           
    Input   : Var a structure which stores the relevant information
             iz, scalar the index of the current shock

    Output  : maxl, 8 by 1 vector, the number on curves to be plotted into each panel

%}

n=length(Var);
maxl=zeros(8,1);
for i=1:1:n;
    if Var(i).Plotno>0 && Var(i).Plotno<=8;
        if (Var(i).Type=='z')&&(Var(i).Pos==iz);
            maxl(Var(i).Plotno)=maxl(Var(i).Plotno)+1;
        end;
        if (Var(i).Type=='z')&&~(Var(i).Pos==iz); continue; end;
        if (Var(i).Type=='x')||(Var(i).Type=='y');
            maxl(Var(i).Plotno)=maxl(Var(i).Plotno)+1;
        end;
    elseif Var(i).Plotno>0;
        message=strcat('CoRRAM supports only 8 panels but Var(',num2str(i),').Plotno=',num2str(Var(i).Plotno));
        error(message);
    end;
end;
% remove non consequtive panels with zeros line
%maxl=nonzeros(maxl);
return;

end

