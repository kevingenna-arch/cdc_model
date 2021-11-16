dynare international_BC nolog

load('international_BC_results.mat', 'oo_');
c1=oo_.endo_simul(1,:);
k2=oo_.endo_simul(4,:);
save simuldata_international c1 k2;