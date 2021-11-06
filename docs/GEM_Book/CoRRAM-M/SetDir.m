function ok=SetDir(Pfad);

% Alfred Mauﬂner
% Revision history
% 14 August 2015, first version
%
% Purpose: set Matlab path to the path specified in the string variable Pfad.
%          If this directory does not exist, Matlab will create it.
%
% Input: Pfad, string variable that holds the name of the path
% Output:  ok, character, ok='T' signals successful operation.

% check if the directory exists, if it exists, set Matlab directory to this
% directory, else, create it and set Matlab directory to this directory.
    
    if exist(Pfad)==7;
        cd(Pfad);
        ok='T';
    else;
        status=mkdir(Pfad); if status; cd(Pfad); ok='T'; else;  ok='F'; return; end;
        fprintf(lfh,strcat('ok=',ok,'\n'));
    end;
    return
end

