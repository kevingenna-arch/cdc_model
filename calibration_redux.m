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
% title('Dynamique du PIB: modele demo');
legend('Actuelle','Precedente','location','best');
saveas(fig_pibf_demo, './plots/demo_pop.eps', 'epsc')

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

% labour force
nbar_n = simuls_new.nbar(241:280, :);
nbar_o = simuls_old.nbar(241:280, :);

% labour force
lq_n = simuls_new.L_Q(241:280, :);
lq_o = simuls_old.L_Q(241:280, :);

% labour force
lnq_n = simuls_new.L_NQ(241:280, :);
lnq_o = simuls_old.L_NQ(241:280, :);

%% Figures
fig_redux_pib = figure(2);
plot(yrs, pib_n, '+-', ...
     yrs, pib_o, '*-');
% title('PIB: modele demosante');
legend('Actuelle','Precedente','location','best');
saveas(fig_redux_pib, './plots/redux_pib.eps', 'epsc');

fig_redux_popret = figure(3);
plot(yrs, popret_n, '+-', ...
     yrs, popret_o, '*-');
% title('Pop. à la retraite: modele demosante');
legend('Actuelle','Precedente','location','best');
saveas(fig_redux_popret, './plots/redux_popret.eps', 'epsc');


fig_redux_totpop = figure(4);
plot(yrs, poptot_n, '+-', ...
     yrs, poptot_o, '*-');
% title('Population totale: modele demosante');
legend('Actuelle','Precedente','location','best');
saveas(fig_redux_totpop, './plots/redux_poptot.eps', 'epsc');


fig_redux_sante = figure(5);
plot(yrs, sante_n, '+-', ...
     yrs, sante_o, '*-');
% title('Stock de santé: modele demosanté');
legend('Actuelle','Precedente','location','best');
saveas(fig_redux_sante, './plots/redux_sante.eps', 'epsc');


fig_redux_nbar = figure(6);
plot(yrs, nbar_n, '+-', ...
     yrs, nbar_o, '*-');
% title('Force travail effective: modele demosanté');
legend('Actuelle','Precedente','location','best');
saveas(fig_redux_nbar, './plots/redux_nbar.eps', 'epsc');

fig_redux_lab = figure(7);
plot(yrs, lq_n, '+-', ...
     yrs, lq_o, '*-', ...
     yrs, lnq_n, 'x-', ...
     yrs, lnq_o, 'o-');
title('Pariticipation au marché du travail');
legend('Qualifiés: Actuelle', 'Qualifiés: Precedente', ... 
    'Non Qualifiés: Actuelle', 'Non Qualifiés: Precedente', ...
    'location', 'best');
% saveas(fig_redux_lab, './plots/redux_lab.eps', 'epsc');

% Plots d'ES
% Ps
pop_es_q    = simuls_new(1, simuls_new.Properties.VariableNames(~cellfun('isempty', regexp(simuls_new.Properties.VariableNames, 'P_Q_', 'once'))));
pop_es_nq   = simuls_new(1, simuls_new.Properties.VariableNames(~cellfun('isempty', regexp(simuls_new.Properties.VariableNames, 'P_NQ_', 'once'))));

% betas
beta_es_q    = simuls_new(1, simuls_new.Properties.VariableNames(~cellfun('isempty', regexp(simuls_new.Properties.VariableNames, 'beta_Q_', 'once'))));
beta_es_nq   = simuls_new(1, simuls_new.Properties.VariableNames(~cellfun('isempty', regexp(simuls_new.Properties.VariableNames, 'beta_NQ_', 'once'))));

% santé
sante_es_q    = simuls_new(1, simuls_new.Properties.VariableNames(~cellfun('isempty', regexp(simuls_new.Properties.VariableNames, 'h_Q_', 'once'))));
sante_es_nq   = simuls_new(1, simuls_new.Properties.VariableNames(~cellfun('isempty', regexp(simuls_new.Properties.VariableNames, 'h_NQ_', 'once'))));

pop_es_q    = table2array(pop_es_q)';
pop_es_nq   = table2array(pop_es_nq)';
beta_es_q   = table2array(beta_es_q)';
beta_es_nq  = table2array(beta_es_nq)';
sante_es_q  = table2array(sante_es_q)';
sante_es_nq = table2array(sante_es_nq)';

cts = (1:17)';

fig_redux_es = figure(8);

subplot(1, 3, 1);
plot(cts, 1 + pop_es_nq, '+', ...
     cts, 1 + pop_es_q, '*');
% title('Etat Stationnaire: Repartition entre Qualifés et Non Qualifé');
legend('Pop. Qualif.', 'Pop. Non Qualif.', 'location', 'best');

subplot(1, 3, 2);
plot((2:17)', beta_es_nq, '+', ...
     (2:17)', beta_es_q, '*');
% title('Etat Stationnaire: Esperance de survie, Qualifés et Non Qualifé');
legend('$\beta$ Qualif.', '$\beta$ Non Qualif.', ...
       'location', 'best', ...
       'interpreter','latex');

subplot(1, 3, 3);
plot(cts, sante_es_nq, '+', ...
     cts, sante_es_q, '*');
% title('Etat Stationnaire: Stock de Santé Qualifés et Non Qualifé');
legend('Qualif.', 'Non Qualif.', 'location', 'best');

saveas(fig_redux_es, './plots/redux_es.eps', 'epsc');