%%% Script pour noyau avec
% - bloc Population
% - bloc Health
% - effets de productivité A et
% - interaction entre stock de santé H et productivité/pib Y
% - effets de deux chocs:
%	+ santé +30% dés 2030 pour les 45-60
%	+ migration +10% dés 2030 pour les 20 
%
% mars 2022, EF

% clean SS
dynare cdc_pha_stable nolog nopreprocessoroutput;

%% Demo
% baseline
dynare cdc_pha_modchocs_cal -DTFP=100 nolog nopreprocessoroutput;
base = simuls(241:280, :);

% pop +10%
dynare cdc_pha_modchocs_cal -Dmig10=1 -DTFP=100 nolog nopreprocessoroutput;
pop10 = simuls(241:280, :);

% health +30%
dynare cdc_pha_modchocs_cal -Dmed30=1 -DTFP=100 nolog nopreprocessoroutput;
health30 = simuls(241:280, :);


% data work
yrs = [1900:5:2095];
% baseline vs pop
plot_pop = figure(1);
subplot(1, 2, 1);
plot(yrs, base.y, '-', ...
	yrs, pop10.y, '--');
legend('Base', 'Reduction chomage', 'Location', 'best');
title('PIB');

subplot(1, 2, 2);
plot(yrs, base.L_NQ, '-', ...
	yrs, pop10.L_NQ, '--', ...
	yrs, base.L_Q, '-', ...
	yrs, pop10.L_Q, '--');
legend('Base NQ', 'Reduction chomage NQ', 'Base Q', 'Reduction chomage Q', 'Location', 'best');
title('Offre de travail');


% baseline vs health
plot_health = figure(2);
subplot(1, 3, 1);
plot(yrs, base.y, '-', ...
	yrs, health30.y, '--');
legend('Base', 'Reduction mortalité', 'Location', 'best');
title('PIB');

subplot(1, 3, 2);
plot(yrs, base.L_NQ, '-', ...
	yrs, health30.L_NQ, '--', ...
	yrs, base.L_Q, '-', ...
	yrs, health30.L_Q, '--');
legend('Base NQ', 'Reduction mortalité NQ', 'Base Q', 'Reduction mortalité Q', 'Location', 'best');
title('Offre de travail');

subplot(1, 3, 3);
plot(yrs, base.nbar, '-', ...
	yrs, health30.nbar, '--');
legend('Base', "Reduction mortalité", ...
	'Location', "best");
title('Travail effectif');

plot_healthstock = figure(3);
subplot(1, 2, 1)
plot(yrs, base.H, '-', ...
	yrs, health30.H, '--');
legend('Base', 'Reduction mortalité', "Location", "best");
title("Stock santé des actifs");

subplot(1, 2, 2);
plot(yrs, base.Htot, '-', ...
	yrs, health30.Htot, '--');
legend('Base', 'Reduction mortalité', "Location", "best");
title('Stock santé totale');

% santé totale VS santé actifs
plot_totsant = figure(4);
plot(yrs, base.y, '-', ...
	yrs, base.y_h, '--');
legend('Base', 'Santé totale', 'Location', 'best');
title('PIB');