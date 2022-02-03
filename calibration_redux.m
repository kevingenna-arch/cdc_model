%% Only population model
% compares the path for output embedding 
% a skeletal production function
% no health effects in this case

dynare cdcredux_p_cal -Dold=1 nolog nopreprocessoroutput;

simuls_old_demo = simuls;
[ro, ~] = size(simuls_old_demo);
simuls_old_demo.vers = repmat('old', ro, 1);
simuls_old_demo.index = ((1:ro)-2)';

dynare cdcredux_p_cal -Dold=0 nolog nopreprocessoroutput;

simuls_new_demo = simuls;
[ro, ~] = size(simuls_new_demo);
simuls_new_demo.vers = repmat('new', ro, 1);
simuls_new_demo.index = ((1:ro)-2)';

yrs = [1900:5:2095];

%% PIB
pib_n_demo = simuls_new_demo.y(241:280, :);
pib_o_demo = simuls_old_demo.y(241:280, :);


fig_pibf_demo = figure(1);
plot(yrs, pib_n_demo, '+-', ...
     yrs, pib_o_demo, '*-');
title('Dynamique du PIB: modele demo');
legend('Actuelle','Precedente','location','best');
% saveas(fig_pibf_demo, './plots/pib_popmodel.eps', 'epsc')

pop_pib_demo = array2table([yrs', pib_n_demo, pib_o_demo], "VariableNames",{'time', 'new', 'old'});
% writetable(pop_pib_demo, ['./output/cdcredux_pop_pib.csv']);

%% Health and demographics

dynare cdcredux_cal -Dold=1 nolog nopreprocessoroutput;

simuls_old = simuls;
[ro, ~] = size(simuls_old);
simuls_old.vers = repmat('old', ro, 1);
simuls_old.index = ((1:ro)-2)';

dynare cdcredux_cal -Dold=0 nolog nopreprocessoroutput;

simuls_new = simuls;
[ro, ~] = size(simuls_new);
simuls_new.vers = repmat('new', ro, 1);
simuls_new.index = ((1:ro)-2)';

yrs = [1900:5:2095];

% PIB
pib_n = simuls_new.y(241:280, :);
pib_o = simuls_old.y(241:280, :);

% Pop retraite
popret_n = simuls_new.Pret(241:280, :);
popret_o = simuls_old.Pret(241:280, :);

% Stock santé
sante_n = simuls_new.H(241:280, :);
sante_o = simuls_old.H(241:280, :);


% Total population
poptot_n = simuls_new.Ptot(241:280, :);
poptot_o = simuls_old.Ptot(241:280, :);


%% Figures
fig_redux_pib = figure(2);
plot(yrs, pib_n, '+-', ...
     yrs, pib_o, '*-');
title('PIB: modele demosante');
legend('Actuelle','Precedente','location','best');
% saveas(fig_redux_pib, './plots/redux_pib.eps', 'epsc');

fig_redux_popret = figure(3);
plot(yrs, popret_n, '+-', ...
     yrs, popret_o, '*-');
title('Pop. à la retraite: modele demosante');
legend('Actuelle','Precedente','location','best');
% saveas(fig_redux_popret, './plots/redux_popret.eps', 'epsc');


fig_redux_totpop = figure(4);
plot(yrs, poptot_n, '+-', ...
     yrs, poptot_o, '*-');
title('Population totale: modele demosante');
legend('Actuelle','Precedente','location','best');
% saveas(fig_redux_poptot, './plots/redux_poptot.eps', 'epsc');


fig_redux_sante = figure(5);
plot(yrs, sante_n, '+-', ...
     yrs, sante_o, '*-');
title('Stock de santé: modele demosanté');
legend('Actuelle','Precedente','location','best');
% saveas(fig_redux_sante, './plots/redux_sante.eps', 'epsc');
