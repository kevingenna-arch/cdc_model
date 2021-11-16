dynare international_BC nolog

load('international_BC_results.mat', 'oo_');
c=oo_.endo_simul(1,:);
save simuldata c;