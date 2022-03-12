%%%% Simulations d'age de depart à la retraite
% noyau demographie, santé, TFP
% mars 2022
% Auteur EF

% valeurs de default (full CDC):
% work = 9, retraite à 60 ans, vie active depuis les 20 ans
% pens = 8, 40 ans de retraite jusqu'à 100 ans
% tot = work+pens = 17

% baseline: retraite à 60 ans
% calcul ES
dynare cdc_pha_stable -Dwork=9 -Dpens=8 nolog nopreprocessoroutput;
% simulations avec chocs
dynare cdc_pha_cal -DTFP=100 -Dwork=9 -Dpens=8 nolog nopreprocessoroutput;

base = simuls;


% genereux: retraite à 55 ans
% calcul ES
dynare cdc_pha_stable -Dwork=8 -Dpens=9 nolog nopreprocessoroutput;
% simulations avec chocs
dynare cdc_pha_cal -DTFP=100 -Dwork=8 -Dpens=9 nolog nopreprocessoroutput;

gen = simuls;

% Macron: retraite à 65 ans
% calcul ES
dynare cdc_pha_stable -Dwork=10 -Dpens=7 nolog nopreprocessoroutput;
% simulations avec chocs
dynare cdc_pha_cal -DTFP=100 -Dwork=10 -Dpens=7 nolog nopreprocessoroutput;

mac = simuls;

% Friedman: retraite à 70 ans
% calcul ES
dynare cdc_pha_stable -Dwork=11 -Dpens=6 nolog nopreprocessoroutput;
% simulations avec chocs
dynare cdc_pha_cal -DTFP=100 -Dwork=11 -Dpens=6 nolog nopreprocessoroutput;

fried = simuls;

yrs = [1900:5:2095];

pib = figure(1);
plot(yrs(25:end), base.y(265:280)/base.y(265), ...
	 yrs(25:end), gen.y(265:280)/gen.y(265), ...
	 yrs(25:end), mac.y(265:280)/mac.y(265), ...
	 yrs(25:end), fried.y(265:280)/fried.y(265));
legend('60', '55', '65', '70')
title('PIB normalisé au niveau 2020')


pib_ratio = figure(2);
plot(yrs(25:end), base.y(265:280)/base.y(265)./base.y(265:280), ...
	 yrs(25:end), gen.y(265:280)/gen.y(265)./base.y(265:280), ...
	 yrs(25:end), mac.y(265:280)/mac.y(265)./base.y(265:280), ...
	 yrs(25:end), fried.y(265:280)/fried.y(265)./base.y(265:280));
legend('60', '55', '65', '70');
title('Rapport PIB vs baseline')