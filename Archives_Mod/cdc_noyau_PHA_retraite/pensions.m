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
legend('60', '55', '65', '70', ...
	'location', 'best');
title('PIB normalisé au niveau 2020')
saveas(pib, '../plots/pens_pha_pib_normalise.eps', 'epsc');

pib_pct = figure(2);
plot(yrs(25:end), 100*((gen.y(265:280)./base.y(265:280))/(gen.y(265)/base.y(265)) - 1), ...
	yrs(25:end), 100*((mac.y(265:280)./base.y(265:280))/(mac.y(265)/base.y(265)) - 1), ...
	yrs(25:end), 100*((fried.y(265:280)./base.y(265:280))/(fried.y(265)/base.y(265)) - 1))
yline(0);
legend('55', '65', '70', ...
	'location', 'best');
title('PIB: ecarts en pourcentage par rapport au scenario actuel - 60 ans')
saveas(pib_pct, '../plots/pens_pha_pib_pct.eps', 'epsc');

pib_niv = figure(3);
plot(yrs, base.y(241:280), ...
	 yrs, gen.y(241:280), ...
	 yrs, mac.y(241:280), ...
	 yrs, fried.y(241:280));
legend('60', '55', '65', '70',...
 'Location', 'best')
title('PIB niveaux')
saveas(pib_niv, '../plots/pens_pha_pib_niveau_no.eps', 'epsc');



%%%%%%%%%%%%% Hypotheses croissance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% baseline: retraite à 60 ans
% calcul ES
dynare cdc_pha_stable -Dwork=9 -Dpens=8 nolog nopreprocessoroutput;
% simulations avec chocs
dynare cdc_pha_cal -DTFP=3 -Dwork=9 -Dpens=8 nolog nopreprocessoroutput;

low_base = simuls;

% genereux: retraite à 55 ans
% calcul ES
dynare cdc_pha_stable -Dwork=8 -Dpens=9 nolog nopreprocessoroutput;
% simulations avec chocs
dynare cdc_pha_cal -DTFP=3 -Dwork=8 -Dpens=9 nolog nopreprocessoroutput;

low_gen = simuls;

% Macron: retraite à 65 ans
% calcul ES
dynare cdc_pha_stable -Dwork=10 -Dpens=7 nolog nopreprocessoroutput;
% simulations avec chocs
dynare cdc_pha_cal -DTFP=3 -Dwork=10 -Dpens=7 nolog nopreprocessoroutput;

low_mac = simuls;

% Friedman: retraite à 70 ans
% calcul ES
dynare cdc_pha_stable -Dwork=11 -Dpens=6 nolog nopreprocessoroutput;
% simulations avec chocs
dynare cdc_pha_cal -DTFP=3 -Dwork=11 -Dpens=6 nolog nopreprocessoroutput;

low_fried = simuls;

yrs = [1900:5:2095];

pib = figure(11);
plot(yrs(25:end), low_base.y(265:280)/low_base.y(265), ...
	 yrs(25:end), low_gen.y(265:280)/low_gen.y(265), ...
	 yrs(25:end), low_mac.y(265:280)/low_mac.y(265), ...
	 yrs(25:end), low_fried.y(265:280)/low_fried.y(265));
legend('60', '55', '65', '70', ...
	'location', 'best');
title('PIB normalisé au niveau 2020')
saveas(pib, '../plots/pens_pha_pib_normalise_noA.eps', 'epsc');

pib_pct = figure(12);
plot(yrs(25:end), 100*((low_gen.y(265:280)./low_base.y(265:280))/(low_gen.y(265)/low_base.y(265)) - 1), ...
	yrs(25:end), 100*((low_mac.y(265:280)./low_base.y(265:280))/(low_mac.y(265)/low_base.y(265)) - 1), ...
	yrs(25:end), 100*((low_fried.y(265:280)./low_base.y(265:280))/(low_fried.y(265)/low_base.y(265)) - 1))
yline(0);
legend('55', '65', '70', ...
	'location', 'best');
title('PIB: ecarts en pourcentage par rapport au scenario actuel - 60 ans')
saveas(pib_pct, '../plots/pens_pha_pib_pct_noA.eps', 'epsc');

pib_niv = figure(13);
plot(yrs, low_base.y(241:280), ...
	 yrs, low_gen.y(241:280), ...
	 yrs, low_mac.y(241:280), ...
	 yrs, low_fried.y(241:280));
legend('60', '55', '65', '70',...
 'Location', 'best')
title('PIB niveaux')
saveas(pib_niv, '../plots/pens_pha_pib_niveau_noA.eps', 'epsc');

%%% individual prods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% baseline
dynare cdc_pha_stable -Dwork=9 -Dpens=8 -Dprodind=1 nolog nopreprocessoroutput;
dynare cdc_pha_cal -Dwork=9 -Dpens=8 -Dprodind=1 -DTFP=3 nolog nopreprocessoroutput;

prodin_base = simuls;
prodin_base.mod = repmat({'base'}, height(prodin_base), 1);
prodin_base.t = (1:height(prodin_base))';
prodin_basef = stack(prodin_base, 1:(width(prodin_base)-2));
prodin_basef.Properties.VariableNames = {'mod', 't', 'variable', 'value'};

% formation continue
dynare cdc_pha_stable -Dwork=9 -Dpens=8 -Dprodind=2 nolog nopreprocessoroutput;
dynare cdc_pha_cal -Dwork=9 -Dpens=8 -Dprodind=2 -DTFP=3 nolog nopreprocessoroutput;

prodin_fc = simuls;
prodin_fc.mod = repmat({'FC'}, height(prodin_fc), 1);
prodin_fc.t = (1:height(prodin_fc))';
prodin_fcf = stack(prodin_fc, 1:(width(prodin_fc)-2));
prodin_fcf.Properties.VariableNames = {'mod', 't', 'variable', 'value'};


% activités socialisées
dynare cdc_pha_stable -Dwork=9 -Dpens=8 -Dprodind=3 nolog nopreprocessoroutput;
dynare cdc_pha_cal -Dwork=9 -Dpens=8 -Dprodind=3 -DTFP=3 nolog nopreprocessoroutput;

prodin_as = simuls;
prodin_as.mod = repmat({'AS'}, height(prodin_as), 1);
prodin_as.t = (1:height(prodin_as))';
prodin_asf = stack(prodin_as, 1:(width(prodin_as)-2));
prodin_asf.Properties.VariableNames = {'mod', 't', 'variable', 'value'};


%% plots
pi_pib = figure(104);
plot(yrs, prodin_base.y(241:280), '-', ...
	 yrs, prodin_fc.y(241:280), '-*', ...
	 yrs, prodin_as.y(241:280), '-+');
legend('Base', 'FC', 'AS', ...
	'location', 'best');
title('PIB en niveaux')
saveas(pi_pib, '../plots/prodind_pha_pib.eps', 'epsc');

pi_nbar = figure(105);
plot(yrs, prodin_base.nbar(241:280), '-', ...
	 yrs, prodin_fc.nbar(241:280), '-*', ...
	 yrs, prodin_as.nbar(241:280), '-+');
legend('Base', 'FC', 'AS', ...
	'location', 'best');
title('Travail effectif')
saveas(pi_nbar, '../plots/prodind_pha_nbar.eps', 'epsc');


prodindsims = vertcat(prodin_basef, prodin_fcf, prodin_asf);
prodindout = unstack(prodindsims, 'value', 'mod');

writetable(prodindout, "./output/prodind_sims_pha.csv");


%%%%% augmentation partie des qualifiés %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simulate the baseline: gauss prod ind, central tfp, retirement at 60
dynare cdc_pha_stable -Dwork=9 -Dpens=8 -Dprodind=1 nolog nopreprocessoroutput;
% simulations avec chocs
dynare cdc_pha_cal -DTFP=100 -Dwork=9 -Dpens=8 -Dprodind=1 nolog nopreprocessoroutput;

q_share30 = simuls(241:280, :);


% up the skilled share to .5
dynare cdc_pha_cal -DTFP=100 -Dwork=9 -Dpens=8 -Dprodind=1 -Dchoc_q=100 nolog nopreprocessoroutput;
q_share50 = simuls(241:280, :);


plot_pib = figure(201);
plot(yrs, q_share30.y, '-', ...
	yrs, q_share50.y, '--');
legend('Base', 'Q = +5%', ...
	'Location', 'southoutside')
title('PIB');
saveas(plot_pib, '../plots/qshare_pha_pib.eps', 'epsc');

plot_l = figure(202);
subplot(1 , 2, 1);
plot(yrs, q_share30.L_NQ, '-', ...
	yrs, q_share50.L_NQ, '--', ...
	yrs, q_share30.L_Q, '-', ...
	yrs, q_share50.L_Q, '--');
legend('Base NQ', 'Q = +5%', 'Base Q', 'Q = +5%', ...
	'Location', 'southoutside')
title('Force Travail');

subplot(1, 2, 2);
plot(yrs, q_share30.nbar, '-', ...
	yrs, q_share50.nbar, '--');
legend('Base', 'Q = +5%', ...
	'Location', 'southoutside')
title('Travail Effectif');
saveas(plot_l, '../plots/qshare_pha_labour.eps', 'epsc');

plot_demo = figure(203);
subplot(1,2,1);
plot(yrs, q_share30.Ptot, '-', ...
	yrs, q_share50.Ptot, '--');
legend('Base', 'Q = +5%', ...
	'Location', 'southoutside')
title('Population');
subplot(1,2,2);
plot(yrs, q_share30.H, '-', ...
	yrs, q_share50.H, '--');
legend('Base', 'Q = +5%', ...
	'Location', 'southoutside')
title('Santé Totale');


saveas(plot_demo, '../plots/qshare_pha_demoh.eps', 'epsc');


