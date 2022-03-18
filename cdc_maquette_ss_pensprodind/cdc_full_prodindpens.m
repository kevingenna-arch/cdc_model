%%%%% SS comparison on full CDC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PENSIONS

% valeurs de default (full CDC):
% work = 9, retraite à 60 ans, vie active depuis les 20 ans
% pens = 8, 40 ans de retraite jusqu'à 100 ans
% tot = work+pens = 17

% baseline: retraite à 60 ans
% calcul ES
dynare sscomp_pensprod -Dwork=9 -Dpens=8 nolog nopreprocessoroutput;
base = simuls(1, :);
base.age = 60;

% genereux: retraite à 55 ans
% calcul ES
dynare sscomp_pensprod -Dwork=8 -Dpens=9 nolog nopreprocessoroutput;
gen = simuls(1, :);
gen.age = 55;

% Macron: retraite à 65 ans
% calcul ES
dynare sscomp_pensprod -Dwork=10 -Dpens=7 nolog nopreprocessoroutput;
mac = simuls(1, :);
mac.age = 65;

% Friedman: retraite à 70 ans
% calcul ES
dynare sscomp_pensprod -Dwork=11 -Dpens=6 nolog nopreprocessoroutput;
fried = simuls(1, :);
fried.age = 70;

