%% Pre-setup

addpath C:\dynare\dynare-5.0\matlab

% run two clean models for SS
%dynare cdc_full_new_ss nolog nopreprocessoroutput;
dynare cdc_pha_stable nolog nopreprocessoroutput;

clear all

% Intro: effet demographique
% old data
dynare cdc_pha -Dold_pop=1 nolog nopreprocessoroutput;
demo_old = simuls(241:280, :);
demo_old.mod = repmat({'INSEE20'}, height(demo_old), 1);

% new data
dynare cdc_pha nolog nopreprocessoroutput;
demo_new = simuls(241:280, :);
demo_new.mod = repmat({'INSEE21'}, height(demo_new), 1);

demograph = vertcat(demo_old, demo_new);
writetable(demograph, '../output/tr0_demo.csv');
allsims.demograph = demograph;

% intro: effet de la technologie
% baseline: dA = 1.3 %
dynare cdc_pha -DTFP=100 nolog nopreprocessoroutput;
tech_base = simuls(241:280, :);
tech_base.mod = repmat({'base'}, height(tech_base), 1);
tech_base.t = (1:height(tech_base))';

% pessimist: dA = 1 %
dynare cdc_pha -DTFP=0 nolog nopreprocessoroutput;
tech_pess = simuls(241:280, :);
tech_pess.mod = repmat({'pess'}, height(tech_pess), 1);
tech_pess.t = (1:height(tech_pess))';


% optimist: 1.8 %
dynare cdc_pha -DTFP=1 nolog nopreprocessoroutput;
tech_opt = simuls(241:280, :);
tech_opt.mod = repmat({'opt'}, height(tech_opt), 1);
tech_opt.t = (1:height(tech_opt))';


% disaster: .25%
dynare cdc_pha -DTFP=3 nolog nopreprocessoroutput;
tech_dis = simuls(241:280, :);
tech_dis.mod = repmat({'low'}, height(tech_dis), 1);
tech_dis.t = (1:height(tech_dis))';

% no forecasts: from 2030 onwards no growth at all
dynare cdc_pha -DTFP=2 nolog nopreprocessoroutput;
tech_hist = simuls(241:280, :);
tech_hist.mod = repmat({'hist'}, height(tech_hist), 1);
tech_hist.t = (1:height(tech_hist))';


tech = vertcat(tech_base, tech_opt, tech_hist, tech_dis, tech_pess);
tech = stack(tech, 1:(width(tech)-2));
tech.Properties.VariableNames = {'mod', 'period', 'variable', 'value'};
allsims.tech = tech;
writetable(tech, '../output/tr0_tech.csv');


%% TR1: quantité travail
% Reforme des retraites: maquette complete

% base
dynare cdc_full_new_ss -Dwork=9 -Dpens=8 nolog nopreprocessoroutput;
pens_base = simuls(2,:);
pens_base.mod = {'age_60'};
pens_base = stack(pens_base, 1:(width(pens_base)-1));

% 55
dynare cdc_full_new_ss -Dwork=8 -Dpens=9 nolog nopreprocessoroutput;
pens_gen = simuls(2,:);
pens_gen.mod = {'age_55'};
pens_gen = stack(pens_gen, 1:(width(pens_gen)-1));

% 65
dynare cdc_full_new_ss -Dwork=10 -Dpens=7 nolog nopreprocessoroutput;
pens_mac = simuls(2,:);
pens_mac.mod = {'age_65'};
pens_mac = stack(pens_mac, 1:(width(pens_mac)-1));

% 70
dynare cdc_full_new_ss -Dwork=11 -Dpens=6 nolog nopreprocessoroutput;
pens_fried = simuls(2,:);
pens_fried.mod = {'age_70'};
pens_fried = stack(pens_fried, 1:(width(pens_fried)-1));

pensions = vertcat(pens_base, pens_gen,...
	pens_mac, pens_fried);
pensions.Properties.VariableNames = {'mod', 'variable', 'value'};
allsims.pensions = pensions;
writetable(pensions, '../output/tr1_pensions_ss.csv')



% Reforme des retraites: maquette complete et deficits

% base
dynare cdc_full_new_ss -Dretrat=.14 -Dwork=9 -Dpens=8 nolog nopreprocessoroutput;
pens_base_def = simuls(2,:);
pens_base_def.mod = {'age_60'};
pens_base_def.retg = .14;
pens_base_def = stack(pens_base_def, 1:(width(pens_base_def)-2));

% 55
dynare cdc_full_new_ss -Dretrat=.17 -Dwork=8 -Dpens=9 -Dnolog nopreprocessoroutput;
pens_gen_def = simuls(2,:);
pens_gen_def.mod = {'age_55'};
pens_gen_def.retg = .17;
pens_gen_def = stack(pens_gen_def, 1:(width(pens_gen_def)-2));

% 65
dynare cdc_full_new_ss -Dretrat=.12 -Dwork=10 -Dpens=7 nolog nopreprocessoroutput;
pens_mac_def = simuls(2,:);
pens_mac_def.mod = {'age_65'};
pens_mac_def.retg = .12;
pens_mac_def = stack(pens_mac_def, 1:(width(pens_mac_def)-2));

% 70
dynare cdc_full_new_ss -Dretrat=.1 -Dwork=11 -Dpens=6 nolog nopreprocessoroutput;
pens_fried_def = simuls(2,:);
pens_fried_def.mod = {'age_70'};
pens_fried_def.retg = .1;
pens_fried_def = stack(pens_fried_def, 1:(width(pens_fried_def)-2));

pensions_def = vertcat(pens_base_def, pens_gen_def, ...
	pens_mac_def, pens_fried_def);
pensions_def.Properties.VariableNames = {'mod', 'def', 'variable', 'value'};
allsims.pensions_def = pensions_def;
writetable(pensions_def, '../output/tr1_pensions_dep_ss.csv');


% reforme retraites dynamique
% base
dynare cdc_pha_stable -Dwork=9 -Dpens=8 nolog nopreprocessoroutput;
dynare cdc_pha -Dwork=9 -Dpens=8 -DTFP=100 nolog nopreprocessoroutput;
pens_dyn_base = simuls(241:280, :);
pens_dyn_base.mod = repmat({'age_60'}, height(pens_dyn_base), 1);
pens_dyn_base.t = (1:height(pens_dyn_base))';
pens_dyn_base = stack(pens_dyn_base, 1:(width(pens_dyn_base)-2));
% 55
dynare cdc_pha_stable -Dwork=8 -Dpens=9 nolog nopreprocessoroutput;
dynare cdc_pha -Dwork=8 -Dpens=9 -DTFP=100 nolog nopreprocessoroutput;
pens_dyn_gen = simuls(241:280, :);
pens_dyn_gen.mod = repmat({'age_55'}, height(pens_dyn_gen), 1);
pens_dyn_gen.t = (1:height(pens_dyn_gen))';
pens_dyn_gen = stack(pens_dyn_gen, 1:(width(pens_dyn_gen)-2));
% 65
dynare cdc_pha_stable -Dwork=10 -Dpens=7 nolog nopreprocessoroutput;
dynare cdc_pha -Dwork=10 -Dpens=7 -DTFP=100 nolog nopreprocessoroutput;
pens_dyn_mac = simuls(241:280, :);
pens_dyn_mac.mod = repmat({'age_65'}, height(pens_dyn_mac), 1);
pens_dyn_mac.t = (1:height(pens_dyn_mac))';
pens_dyn_mac = stack(pens_dyn_mac, 1:(width(pens_dyn_mac)-2));
% 70
dynare cdc_pha_stable -Dwork=11 -Dpens=6 nolog nopreprocessoroutput;
dynare cdc_pha -Dwork=11 -Dpens=6 -DTFP=100 nolog nopreprocessoroutput;
pens_dyn_fried = simuls(241:280, :);
pens_dyn_fried.mod = repmat({'age_70'}, height(pens_dyn_fried), 1);
pens_dyn_fried.t = (1:height(pens_dyn_fried))';
pens_dyn_fried = stack(pens_dyn_fried, 1:(width(pens_dyn_fried)-2));


pensions_dyna = vertcat(pens_dyn_base, pens_dyn_gen, ...
	pens_dyn_mac, pens_dyn_fried);
pensions_dyna.Properties.VariableNames = {'mod', 'period', 'variable', 'value'};
allsims.pensions_dyna = pensions_dyna;
writetable(pensions_dyna, '../output/tr1_pensions_dyna.csv');

% reforme retraites dynamique - low productivity
% base
dynare cdc_pha -Dwork=9 -Dpens=8 -DTFP=3 nolog nopreprocessoroutput;
pens_dyn_base_lp = simuls(241:280, :);
pens_dyn_base_lp.mod = repmat({'age_60'}, height(pens_dyn_base_lp), 1);
pens_dyn_base_lp.t = (1:height(pens_dyn_base_lp))';
pens_dyn_base_lp = stack(pens_dyn_base_lp, 1:(width(pens_dyn_base_lp)-2));
% 55
dynare cdc_pha -Dwork=8 -Dpens=9 -DTFP=3 nolog nopreprocessoroutput;
pens_dyn_gen_lp = simuls(241:280, :);
pens_dyn_gen_lp.mod = repmat({'age_55'}, height(pens_dyn_gen_lp), 1);
pens_dyn_gen_lp.t = (1:height(pens_dyn_gen_lp))';
pens_dyn_gen_lp = stack(pens_dyn_gen_lp, 1:(width(pens_dyn_gen_lp)-2));
% 65
dynare cdc_pha -Dwork=10 -Dpens=7 -DTFP=3 nolog nopreprocessoroutput;
pens_dyn_mac_lp = simuls(241:280, :);
pens_dyn_mac_lp.mod = repmat({'age_65'}, height(pens_dyn_mac_lp), 1);
pens_dyn_mac_lp.t = (1:height(pens_dyn_mac_lp))';
pens_dyn_mac_lp = stack(pens_dyn_mac_lp, 1:(width(pens_dyn_mac_lp)-2));
% 70
dynare cdc_pha -Dwork=11 -Dpens=6 -DTFP=3 nolog nopreprocessoroutput;
pens_dyn_fried_lp = simuls(241:280, :);
pens_dyn_fried_lp.mod = repmat({'age_70'}, height(pens_dyn_fried_lp), 1);
pens_dyn_fried_lp.t = (1:height(pens_dyn_fried_lp))';
pens_dyn_fried_lp = stack(pens_dyn_fried_lp, 1:(width(pens_dyn_fried_lp)-2));


pensions_dyna_lp = vertcat(pens_dyn_base_lp, pens_dyn_gen_lp, ...
	pens_dyn_mac_lp, pens_dyn_fried_lp);
pensions_dyna_lp.Properties.VariableNames = {'mod', 'period', 'variable', 'value'};
allsims.pensions_dyna_lp = pensions_dyna_lp;
writetable(pensions_dyna_lp, '../output/tr1_pensions_dyna_lowprod.csv');

% +10% of migration as U reduction

% take previous baseline
dynare cdc_pha -DTFP=100 nolog nopreprocessoroutput;
chom_base = simuls(241:280, :);
chom_base.mod = repmat({'base'}, height(chom_base), 1);
chom_base.t = (1:height(chom_base))';

% +10% in working pop
dynare cdc_pha -DTFP=100 -Dmig10=0 nolog nopreprocessoroutput;
chom_red = simuls(241:280, :);
chom_red.mod = repmat({'shock'}, height(chom_red), 1);
chom_red.t = (1:height(chom_red))';

chom = vertcat(chom_base, chom_red);
chom = stack(chom, 1:(width(chom)-2));
chom.Properties.VariableNames = {'mod', 'period', 'variable', 'value'};
allsims.chom = chom;
writetable(chom, '../output/tr1_labour_mig.csv');


%% TR 2: individual productivity, skill shock

% SS comparison
% baseline
dynare cdc_full_base_prodind nolog nopreprocessoroutput;
prodind_full_base = simuls(2,:);
prodind_full_base.mod = {'base'};

% formation continue
dynare cdc_full_base_prodind -Dprodind=2 nolog nopreprocessoroutput;
prodind_full_fc = simuls(2,:);
prodind_full_fc.mod = {'FC'};

% activités socialisées
dynare cdc_full_base_prodind -Dprodind=3 nolog nopreprocessoroutput;
prodind_full_as = simuls(2,:);
prodind_full_as.mod = {'AS'};

% FC+AS
dynare cdc_full_base_prodind -Dprodind=4 nolog nopreprocessoroutput;
prodind_full_bo = simuls(2,:);
prodind_full_bo.mod = {'bo'};

prodind_full = vertcat(prodind_full_base, prodind_full_fc, prodind_full_as, prodind_full_bo);
prodind_full = stack(prodind_full, 1:(width(prodind_full)-1));
prodind_full.Properties.VariableNames = {'mod', 'variable', 'value'};
allsims.prodind_full = prodind_full;
writetable(prodind_full, '../output/tr2_skill_ss.csv');

% dynamic effects
% baseline
dynare cdc_pha_stable nolog nopreprocessoroutput;
dynare cdc_pha -DTFP=100 nolog nopreprocessoroutput;
prodind_dyn_base = simuls(241:280, :);
prodind_dyn_base.mod = repmat({'base'}, height(prodind_dyn_base), 1);
prodind_dyn_base.t = (1:height(prodind_dyn_base))';
prodind_dyn_base = stack(prodind_dyn_base, 1:(width(prodind_dyn_base)-2));

% fc
dynare cdc_pha_stable -Dprodind=2 nolog nopreprocessoroutput;
dynare cdc_pha -DTFP=100 -Dprodind=2 nolog nopreprocessoroutput;
prodind_dyn_fc = simuls(241:280, :);
prodind_dyn_fc.mod = repmat({'fc'}, height(prodind_dyn_fc), 1);
prodind_dyn_fc.t = (1:height(prodind_dyn_fc))';
prodind_dyn_fc = stack(prodind_dyn_fc, 1:(width(prodind_dyn_fc)-2));

% as
dynare cdc_pha_stable -Dprodind=3 nolog nopreprocessoroutput;
dynare cdc_pha -DTFP=100 -Dprodind=3 nolog nopreprocessoroutput;
prodind_dyn_as = simuls(241:280, :);
prodind_dyn_as.mod = repmat({'as'}, height(prodind_dyn_as), 1);
prodind_dyn_as.t = (1:height(prodind_dyn_as))';
prodind_dyn_as = stack(prodind_dyn_as, 1:(width(prodind_dyn_as)-2));

% as+fc
dynare cdc_pha_stable -Dprodind=4 nolog nopreprocessoroutput;
dynare cdc_pha -DTFP=100 -Dprodind=4 nolog nopreprocessoroutput;
prodind_dyn_bo = simuls(241:280, :);
prodind_dyn_bo.mod = repmat({'bo'}, height(prodind_dyn_bo), 1);
prodind_dyn_bo.t = (1:height(prodind_dyn_bo))';
prodind_dyn_bo = stack(prodind_dyn_bo, 1:(width(prodind_dyn_bo)-2));


% join all
prodind_dyn = vertcat(prodind_dyn_base, prodind_dyn_fc, prodind_dyn_as,prodind_dyn_bo);
prodind_dyn.Properties.VariableNames = {'mod', 'period', 'variable', 'value'};
allsims.prodind_dyn = prodind_dyn;
writetable(prodind_dyn, '../output/tr2_skill_dyn.csv');


% tick up share of skilled
% baseline
dynare cdc_pha_stable nolog nopreprocessoroutput;
dynare cdc_pha -DTFP=100 nolog nopreprocessoroutput;
share_base = simuls(241:280, :);
share_base.mod = repmat({'base'}, height(share_base), 1);
share_base.t = (1:height(share_base))';

% share tick
dynare cdc_pha -DTFP=100 -Dq_shock=3 nolog nopreprocessoroutput;
share_up = simuls(241:280, :);
share_up.mod = repmat({'shock'}, height(share_up), 1);
share_up.t = (1:height(share_up))';

share = vertcat(share_base, share_up);
share = stack(share, 1:(width(share)-2));
share.Properties.VariableNames = {'mod', 'period', 'variable', 'value'};
allsims.share = share;
writetable(share, '../output/tr2_skilledshock_dyn.csv');


%% Table Ronde 3: santé
% baseline
dynare cdc_pha -DTFP=100 nolog nopreprocessoroutput;
health_base = simuls(241:280, :);
health_base.mod = repmat({'base'}, height(health_base), 1);
health_base.t = (1:height(health_base))';


% 50% more heath for mid age 45-65
dynare cdc_pha -DTFP=100 -Dmed30=1 nolog nopreprocessoroutput;
health_mid = simuls(241:280, :);
health_mid.mod = repmat({'45-65'}, height(health_mid), 1);
health_mid.t = (1:height(health_mid))';

% 50% more health for early age 25-35
dynare cdc_pha -DTFP=100 -Dmed30=2 nolog nopreprocessoroutput;
health_earl = simuls(241:280, :);
health_earl.mod = repmat({'25-35'}, height(health_earl), 1);
health_earl.t = (1:height(health_earl))';

health = vertcat(health_base, health_mid, health_earl);
health = stack(health, 1:(width(health)-2));
health.Properties.VariableNames = {'mod', 'period', 'variable', 'value'};
allsims.health = health;
writetable(health, '../output/tr3_health_dyn.csv');