%%%%%% Script pour la note de fevrier 2022


% run the SS, the EXO, the ENDO scripts

% clean SS
dynare cdcredux_stable_popbetaA nolog nopreprocessoroutput;

% scenario baseline: productivité fixée à A=1
dynare cdcredux_cal_popbetasA_shocks nolog nopreprocessoroutput;

sims_baseline = simuls;

% scenario pessimiste: dA = 1% dés 2045
dynare cdcredux_cal_popbetasA_shocks -DTFP=0 nolog nopreprocessoroutput;

sims_pess = simuls;

% scenario optimiste: dA = 1.8% dés 2045
dynare cdcredux_cal_popbetasA_shocks -DTFP=1 nolog nopreprocessoroutput;

sims_opt = simuls;

% scenario central: dA = 1.4% dés 2045
dynare cdcredux_cal_popbetasA_shocks -DTFP=100 nolog nopreprocessoroutput;

sims_ctr = simuls;

% scenario sans previsions: dA = 0% dés 2030
dynare cdcredux_cal_popbetasA_shocks -DTFP=2 nolog nopreprocessoroutput;

sims_nocro = simuls;

% scenario sans previsions: dA = 0.5% dés 2030
dynare cdcredux_cal_popbetasA_shocks -DTFP=3 nolog nopreprocessoroutput;

sims_des = simuls;

% x axis
yrs = [1900:5:2095];

% PIB
sim_y.base = sims_baseline.y(241:280);
sim_y.pess = sims_pess.y(241:280);
sim_y.opt = sims_opt.y(241:280);
sim_y.ctr = sims_ctr.y(241:280);
sim_y.nocr = sims_nocro.y(241:280);
sim_y.des = sims_des.y(241:280);

sim_A.base = sims_baseline.A(241:280);
sim_A.pess = sims_pess.A(241:280);
sim_A.opt = sims_opt.A(241:280);
sim_A.ctr = sims_ctr.A(241:280);
sim_A.nocr = sims_nocro.A(241:280);
sim_A.des = sims_des.A(241:280);

sim_aq.base = sims_baseline.A_Q(241:280);
sim_aq.pess = sims_pess.A_Q(241:280);
sim_aq.opt = sims_opt.A_Q(241:280);
sim_aq.ctr = sims_ctr.A_Q(241:280);
sim_aq.nocr = sims_nocro.A_Q(241:280);
sim_aq.des = sims_des.A_Q(241:280);


%% Figures
pibs = figure(1);
plot(yrs, sim_y.base, '-', ...
     yrs, sim_y.pess, '-*', ...
     yrs, sim_y.opt, '-+', ...
     yrs, sim_y.ctr, '--', ...
     yrs, sim_y.nocr, '-.', ...
     yrs, sim_y.des, '--s');
legend('Sans PT', 'dA=1%', 'dA=1.8%', 'dA=1.4%', 'dA=0%', 'dA=.25%', ...
     'Location', 'best');
title('PIB');
saveas(pibs, './plots/redux_PHA_pib.eps', 'epsc');
% exportgraphics(pibs,'./plots/redux_PHA_pibs.jpg','Resolution',600);

pts_macro = figure(2);
subplot(1, 2, 1);
plot(yrs, sim_A.base, '-', ...
     yrs, sim_A.pess, '-*', ...
     yrs, sim_A.opt, '-+', ...
     yrs, sim_A.ctr, '--', ...
     yrs, sim_A.nocr, '-.', ...
     yrs, sim_A.des, '--s');
legend('Sans PT', 'dA=1%', 'dA=1.8%', 'dA=1.4%', 'dA=0%', 'dA=.25%', ...
     'Location', 'best');
title('Productivité agregée')
subplot(1, 2, 2);
plot(yrs, sim_aq.base, '-', ...
     yrs, sim_aq.pess, '-*', ...
     yrs, sim_aq.opt, '-+', ...
     yrs, sim_aq.ctr, '--', ...
     yrs, sim_aq.nocr, '-.', ...
     yrs, sim_aq.des, '--s');
legend('Sans PT', 'dA=1%', 'dA=1.8%', 'dA=1.4%', 'dA=0%', 'dA=.25%', ...
     'Location', 'best');
title('Productivité agregée des qualifiés')
saveas(pts_macro, './plots/redux_PHA_TFPs.eps', 'epsc');
% exportgraphics(pts_macro,'./plots/redux_PHA_TFPs.jpg','Resolution',600);

pts_micro = figure(3);
plot(1:9, .7 + log(1:9)/10, '-', ...
     1:9, 1 + log(1:9)/10, '-');
legend('NQ', 'Q', 'location', 'best');
title('Profils de productivités individuelles');
saveas(pts_micro, './plots/redux_PHA_prodmicro.eps', 'epsc');
% exportgraphics(pts_micro,'./plots/redux_PHA_prodmicro.jpg','Resolution',600)