%% Only population model
% compares the path for output embedding 
% a skeletal production function
% no health effects in this case

dynare cdcredux_p_cal -Dold=1 nolog nopreprocessoroutput;

simuls_old = simuls;
[ro, ~] = size(simuls_old);
simuls_old.vers = repmat('old', ro, 1);
simuls_old.index = ((1:ro)-2)';

dynare cdcredux_p_cal -Dold=0 nolog nopreprocessoroutput;

simuls_new = simuls;
[ro, ~] = size(simuls_new);
simuls_new.vers = repmat('new', ro, 1);
simuls_new.index = ((1:ro)-2)';

yrs = [1900:5:2095];

%% PIB
pib_n = simuls_new.y(241:280, :);
pib_o = simuls_old.y(241:280, :);


fig_pibf = figure(2);
plot(yrs, pib_n, '+-', ...
    yrs, pib_o, '*-');
title('Dynamique du PIB: deux calibrations');
legend('Actuelle','Precedente','location','best');
saveas(fig_pibf, './plots/pib_popmodel.eps', 'epsc')

pop_pib = array2table([yrs', pib_n, pib_o], "VariableNames",{'time', 'new', 'old'});
writetable(pop_pib, ['./output/cdcredux_pop_pib.csv']);
