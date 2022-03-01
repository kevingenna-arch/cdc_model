%%%%%% Script pour la note de fevrier 2022


% run the SS, the EXO, the ENDO scripts

% clean SS
dynare cdcredux_stable nolog nopreprocessoroutput;

dynare cdcredux_cal_popbetas -Dold=0 nolog nopreprocessoroutput;

simuls_exo = simuls;
[ro, ~] = size(simuls_exo);
simuls_exo.vers = repmat('exo', ro, 1);
simuls_exo.index = ((1:ro)-2)';

dynare cdcredux_cal_popbetas_shocks nolog nopreprocessoroutput;

simuls_endo = simuls;
[ro, ~] = size(simuls_endo);
simuls_endo.vers = repmat('endo', ro, 1);
simuls_endo.index = ((1:ro)-2)';

yrs = [1900:5:2095];

% PIB
pib_n = simuls_endo.y(241:280, :);
pib_x = simuls_exo.y(241:280, :);

% Pop retraite
popret_n = simuls_endo.Pret(241:280, :);
popret_x = simuls_exo.Pret(241:280, :);

% Stock santé
sante_n = simuls_endo.H(241:280, :);
sante_x = simuls_exo.H(241:280, :);


% Total population
poptot_n = simuls_endo.Ptot(241:280, :);
poptot_x = simuls_exo.Ptot(241:280, :);

% labour force
nbar_n = simuls_endo.nbar(241:280, :);
nbar_x = simuls_exo.nbar(241:280, :);

% labour force
lq_n = simuls_endo.L_Q(241:280, :);
lq_x = simuls_exo.L_Q(241:280, :);

% labour force
lnq_n = simuls_endo.L_NQ(241:280, :);
lnq_x = simuls_exo.L_NQ(241:280, :);

%% Figures
fig_redux_pib = figure(1);
plot(yrs, pib_n, '+-', ...
     yrs, pib_x, 'o-');
title('PIB: ');
legend('Endo Y','Exo Y','location','best');
% saveas(fig_redux_pib, './plots/redux_pib.eps', 'epsc');

fig_redux_popret = figure(2);
plot(yrs, popret_n, '+-', ...
     yrs, popret_x, 'o-');
title('Pop. à la retraite: modele demosante');
legend('Endo Pret','Exo Pret','location','best');
% saveas(fig_redux_popret, './plots/redux_popret.eps', 'epsc');


fig_redux_totpop = figure(3);
plot(yrs, poptot_n, '+-', ...
     yrs, poptot_x, 'o-', ...
     yrs, s_Ptot, '--');
title('Population totale');
legend('Endo Ptot','Exo Ptot', 'Data', 'location','best');
% saveas(fig_redux_totpop, './plots/redux_poptot.eps', 'epsc');


fig_redux_sante = figure(4);
plot(yrs, sante_n, '+-', ...
     yrs, sante_x, 'o-');
title('Stock de santé');
legend('Endo H','Exo H','location','best');
% saveas(fig_redux_sante, './plots/redux_sante.eps', 'epsc');


fig_redux_nbar = figure(5);
plot(yrs, nbar_n, '+-', ...
     yrs, nbar_x, 'o-');
title('Force travail effective');
legend('Endo nbar','Exo nbar','location','best');
% saveas(fig_redux_nbar, './plots/redux_nbar.eps', 'epsc');

fig_redux_lab = figure(6);
plot(yrs, lq_n, '+-', ...
     yrs, lq_x, '*-', ...
     yrs, lnq_n, 'x-', ...
     yrs, lnq_x, 'o-');
title('Participation au marché du travail');
legend('Qualifiés: endo L_Q', 'Qualifiés: exo L_Q', ... 
    'Non Qualifiés: endo L_NQ', 'Non Qualifiés: exo L_NQ', ...
    'location', 'best');
% saveas(fig_redux_lab, './plots/redux_lab.eps', 'epsc');

% Plots d'ES
% Ps
pop_es_q    = simuls_endo(1, simuls_endo.Properties.VariableNames(~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'P_Q_', 'once'))));
pop_es_nq   = simuls_endo(1, simuls_endo.Properties.VariableNames(~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'P_NQ_', 'once'))));

% betas
beta_es_q    = simuls_endo(1, simuls_endo.Properties.VariableNames(~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'beta_Q_', 'once'))));
beta_es_nq   = simuls_endo(1, simuls_endo.Properties.VariableNames(~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'beta_NQ_', 'once'))));

% santé
sante_es_q    = simuls_endo(1, simuls_endo.Properties.VariableNames(~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'h_Q_', 'once'))));
sante_es_nq   = simuls_endo(1, simuls_endo.Properties.VariableNames(~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'h_NQ_', 'once'))));

pop_es_q    = table2array(pop_es_q)';
pop_es_nq   = table2array(pop_es_nq)';
beta_es_q   = table2array(beta_es_q)';
beta_es_nq  = table2array(beta_es_nq)';
sante_es_q  = table2array(sante_es_q)';
sante_es_nq = table2array(sante_es_nq)';

cts = (1:17)';

fig_redux_es = figure(7);

subplot(1, 3, 1);
plot(cts, pop_es_nq, '+', ...
     cts, pop_es_q, 'o');
% title('Etat Stationnaire: Repartition entre Qualifés et Non Qualifé');
legend('Pop. Non Qualif.', 'Pop. Qualif.', 'location', 'best');

subplot(1, 3, 2);
plot((2:17)', beta_es_nq, '+', ...
     (2:17)', beta_es_q, 'o');
% title('Etat Stationnaire: Esperance de survie, Qualifés et Non Qualifé');
legend('$\beta$ Non Qualif.', '$\beta$ Qualif.', ...
       'location', 'best', ...
       'interpreter','latex');

subplot(1, 3, 3);
plot(cts, sante_es_nq, '+', ...
     cts, sante_es_q, 'o');
% title('Etat Stationnaire: Stock de Santé Qualifés et Non Qualifé');
legend('Non Qualif.', 'Qualif.', 'location', 'best');
% saveas(fig_redux_es, './plots/redux_es.eps', 'epsc');

fit = figure(8);
subplot(1, 1, 1);
plot(yrs, s_Ptot-.3*s_xpop, '-', ...
     yrs, poptot_n, ':', ...
     yrs, poptot_x, '.');
legend('Data', 'Model Endo', 'Model Exo', 'location', 'best');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Endogeneous/Exogeneous comparisons

hs_nq10 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'h_NQ_[0-9]{2}', 'once'));
hs_nq1 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'h_NQ_\d(?!\S)', 'once'));
hnq=simuls_endo(241:280, ...
    [sort(simuls_endo.Properties.VariableNames(hs_nq1)), ...
    sort(simuls_endo.Properties.VariableNames(hs_nq10))]);


hs_q10 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'h_Q_[0-9]{2}', 'once'));
hs_q1 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'h_Q_\d(?!\S)', 'once'));
hq=simuls_endo(241:280, ...
    [sort(simuls_endo.Properties.VariableNames(hs_q1)), ...
    sort(simuls_endo.Properties.VariableNames(hs_q10))]);

b_nq10 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'beta_NQ_[0-9]{2}', 'once'));
b_nq1 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'beta_NQ_\d(?!\S)', 'once'));
bnq=simuls_endo(241:280, ...
    [sort(simuls_endo.Properties.VariableNames(b_nq1)), ...
    sort(simuls_endo.Properties.VariableNames(b_nq10))]);


b_q10 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'beta_Q_[0-9]{2}', 'once'));
b_q1 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'beta_Q_\d(?!\S)', 'once'));
bq=simuls_endo(241:280, ...
    [sort(simuls_endo.Properties.VariableNames(b_q1)), ...
    sort(simuls_endo.Properties.VariableNames(b_q10))]);


sante = figure(9);
subplot(2, 2, 1)
plot(yrs, table2array(hnq))
legend(hnq.Properties.VariableNames, 'Location', 'best', 'Interpreter', 'none')
subplot(2, 2, 2)
plot(yrs, table2array(hq))
legend(hq.Properties.VariableNames,'Location', 'best', 'Interpreter', 'none')
subplot(2, 2, 3)
plot(yrs, table2array(bnq))
legend(bnq.Properties.VariableNames, 'Location', 'best', 'Interpreter', 'none')
subplot(2, 2, 4)
plot(yrs, table2array(bq))
legend(bq.Properties.VariableNames,'Location', 'best', 'Interpreter', 'none')

reste = figure(10);
subplot(1, 4, 1);
plot(yrs, simuls_endo.y(241:280))
title('PIB');
subplot(1, 4, 2);
plot(yrs, simuls_endo.H(241:280));
title('Santé')
subplot(1, 4, 3);
plot(yrs, simuls_endo.nbar(241:280));
title('pop active')
subplot(1, 4, 4)
plot(yrs, (simuls_endo.Ptot(241:280)), ...
    yrs, (s_Ptot));
legend('Model', 'Data', 'Location','best');
title('P Tot');


% Shocks

med_q1 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'med_Q_\d(?!\S)', 'once'));
med_q10 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'med_Q_[0-9]{2}', 'once'));
med_q = simuls_endo(241:280, ...
                [sort(simuls_endo.Properties.VariableNames(med_q1)), ...
                 sort(simuls_endo.Properties.VariableNames(med_q10))]);


med_nq1 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'med_NQ_\d(?!\S)', 'once'));
med_nq10 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'med_NQ_[0-9]{2}', 'once'));
med_nq = simuls_endo(241:280, ...
                [sort(simuls_endo.Properties.VariableNames(med_nq1)), ...
                 sort(simuls_endo.Properties.VariableNames(med_nq10))]);



mig_q1 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'mig_Q_\d(?!\S)', 'once'));
mig_q10 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'mig_Q_[0-9]{2}', 'once'));
mig_q = simuls_endo(241:280, ...
                [sort(simuls_endo.Properties.VariableNames(mig_q1)), ...
                 sort(simuls_endo.Properties.VariableNames(mig_q10))]);


mig_nq1 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'mig_NQ_\d(?!\S)', 'once'));
mig_nq10 = ~cellfun('isempty', regexp(simuls_endo.Properties.VariableNames, 'mig_NQ_[0-9]{2}', 'once'));
mig_nq = simuls_endo(241:280, ...
                [sort(simuls_endo.Properties.VariableNames(mig_nq1)), ...
                 sort(simuls_endo.Properties.VariableNames(mig_nq10))]);


shocks = figure(11);
subplot(2, 2, 1)
plot(yrs, table2array(med_q))
title('Med Q');
legend(med_q.Properties.VariableNames, 'Location', 'best', 'Interpreter', 'none')
subplot(2, 2, 2)
plot(yrs, table2array(med_nq))
title('Med NQ');
legend(med_nq.Properties.VariableNames,'Location', 'best', 'Interpreter', 'none')
subplot(2, 2, 3)
plot(yrs, table2array(mig_q))
title('Mig Q');
legend(mig_q.Properties.VariableNames, 'Location', 'best', 'Interpreter', 'none')
subplot(2, 2, 4)
plot(yrs, table2array(mig_nq))
title('Mig NQ');
legend(mig_nq.Properties.VariableNames,'Location', 'best', 'Interpreter', 'none')
