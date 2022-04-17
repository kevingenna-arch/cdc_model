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

% handle results
base = stack(base, 1:width(base));
base.Properties.VariableNames = {'value', 'variable'};
base.Properties.VariableNames = {'variable', 'value'};
base.model = repmat({'base'}, height(base), 1);

gen = stack(gen, 1:width(gen));
gen.Properties.VariableNames = {'value', 'variable'};
gen.Properties.VariableNames = {'variable', 'value'};
gen.model = repmat({'55'}, height(gen), 1);

mac = stack(mac, 1:width(mac));
mac.Properties.VariableNames = {'value', 'variable'};
mac.Properties.VariableNames = {'variable', 'value'};
mac.model = repmat({'65'}, height(mac), 1);

fried = stack(fried, 1:width(fried));
fried.Properties.VariableNames = {'value', 'variable'};
fried.Properties.VariableNames = {'variable', 'value'};
fried.model = repmat({'70'}, height(fried), 1);

out = vertcat(base, gen, mac, fried);

writetable(out, "./output/pens_ss.csv");


%% Individual productivity

% baseline profile
dynare sscomp_pensprod -Dwork=9 -Dpens=8 -Dprodind=1 nolog nopreprocessoroutput;
prindbase = simuls(1, :);

% formation continue
dynare sscomp_pensprod -Dwork=9 -Dpens=8 -Dprodind=2 nolog nopreprocessoroutput;
prindfc = simuls(1, :);

% activités socialisées
dynare sscomp_pensprod -Dwork=9 -Dpens=8 -Dprodind=3 nolog nopreprocessoroutput;
prindas = simuls(1, :);

% classify models
prindbase.mod = {'base'};
prindfc.mod = {'FC'};
prindas.mod = {'AS'};

% put all together
prodind = vertcat(prindbase, prindfc, prindas);
% flip to column
prodind_mid = stack(prodind, 1:(width(prodind)-1));
prodind_mid.Properties.VariableNames = {'mod', 'variable', 'value'};
prodind_out = unstack(prodind_mid, 'value', 'mod');

writetable(prodind_out, "./output/prodind_ss.csv");