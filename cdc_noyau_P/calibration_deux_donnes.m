%%%% M file pour simuler et comparer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% simulations sur l'ancienne calibration
dynare cdc_prev nolog nopreprocessoroutput 

simuls_old = simuls;
[ro, ~] = size(simuls_old);
simuls_old.vers = repmat('old', ro, 1);
simuls_old.index = ((1:ro)-2)';


% simulations sur la nouvelle calibration
dynare cdc_ph nolog nopreprocessoroutput 

simuls_new = simuls;
simuls_new.vers = repmat('new', ro, 1);
[ro, ~] = size(simuls_new);
simuls_new.index = ((1:ro)-2)';

simuls = vertcat(simuls_old, simuls_new);
clearvars -except simuls_new simuls_old simuls s_*;

writetable(simuls, './output/simulations.csv');

yrs = [1900:5:2095];

% valeurs interessants
% PIB
pib_n = simuls_new.y(242:281, :);
pibn_norm = normalize(pib_n, 'range');

pib_o = simuls_old.y(242:281, :);
pibo_norm = normalize(pib_o, 'range');

% Population
pop_n = simuls_new.Ptot(242:281, :) + s_xpop;
poptot_norm = normalize(pop_n, 'range');
popproj = normalize(s_Ptot - s_xpop, 'range');


% santé
stoh_n = simuls_new.H(242:281, :);
stohn_norm = normalize(stoh_n, 'range');

stoh_o = simuls_old.H(242:281, :);
stoho_norm = normalize(stoh_o, 'range');

fig_pop = figure(1);
plot(yrs, poptot_norm, '+-', ...
    yrs, popproj, 'o-');
title('Demographie: projections INSEE et modèle calibré')
legend('Modèle','Données','location','best');


fig_pibf = figure(2);
plot(yrs, pibn_norm, '+-', ...
    yrs, pibo_norm, '*-');
title('Dynamique du PIB: deux calibrations');
legend('Actuelle','Precedente','location','best');

fig_sant = figure(3);
plot(yrs, stohn_norm, '+-', ...
    yrs, stoho_norm, '*-');
title('Dynamique du stock de santé: deux calibrations');
legend('Actuelle','Precedente','location','best');

saveas(fig_pop, './plots/pop_comp.eps', 'epsc');
saveas(fig_pibf, './plots/pibfit_comp.eps', 'epsc');
saveas(fig_sant, './plots/sante_comp.eps', 'epsc');

