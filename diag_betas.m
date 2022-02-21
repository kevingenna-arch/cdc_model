hs_nq = ~cellfun('isempty', regexp(simuls.Properties.VariableNames, 'h_NQ', 'once'));
hnq=simuls(2:41, sort(simuls.Properties.VariableNames(hs_nq)));
hs_q = ~cellfun('isempty', regexp(simuls.Properties.VariableNames, 'h_Q', 'once'));
hq=simuls(2:41, sort(simuls.Properties.VariableNames(hs_q)));

b_nq = ~cellfun('isempty', regexp(simuls.Properties.VariableNames, 'beta_NQ', 'once'));
bnq=simuls(2:41, sort(simuls.Properties.VariableNames(b_nq)));
b_q = ~cellfun('isempty', regexp(simuls.Properties.VariableNames, 'beta_Q', 'once'));
bq=simuls(2:41, sort(simuls.Properties.VariableNames(b_q)));

subplot(2, 2, 1)
plot(1900:5:2095, table2array(hnq))
legend(hnq.Properties.VariableNames, 'Location', 'best', 'Interpreter', 'none')
subplot(2, 2, 2)
plot(1900:5:2095, table2array(hq))
legend(hq.Properties.VariableNames,'Location', 'best', 'Interpreter', 'none')
subplot(2, 2, 3)
plot(1900:5:2095, table2array(bnq))
legend(bnq.Properties.VariableNames, 'Location', 'best', 'Interpreter', 'none')
subplot(2, 2, 4)
plot(1900:5:2095, table2array(bq))
legend(bq.Properties.VariableNames,'Location', 'best', 'Interpreter', 'none')
