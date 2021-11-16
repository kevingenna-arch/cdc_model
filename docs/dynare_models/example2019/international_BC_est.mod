var c1 c2 k1 k2 a1 a2;
varexo e1 e2;

parameters gamma delta alpha beta rho;

gamma=2;
delta=.05;
alpha=.4;
beta=.98;
rho=.85;

model;
c1=c2; exp(c1)^(-gamma) =beta*exp(c1(+1))^(-gamma)*(alpha*exp(a1(+1))*exp(k1)^(alpha-1)+1-delta);
exp(c2)^(-gamma) =beta*exp(c2(+1))^(-gamma)*(alpha*exp(a2(+1))*exp(k2)^(alpha-1)+1-delta);
exp(c1)+exp(c2)+exp(k1)-exp(k1(-1))*(1-delta)+exp(k2)-exp(k2(-1))*(1-delta)= exp(a1)*exp(k1(-1))^alpha+exp(a2)*exp(k2(-1))^alpha;
a1=rho*a1(-1)+e1; 
a2=rho*a2(-1)+e2; 
end;

initval;
k1=2.8;
k2=2.8;
c1=.8;
c2=.8;
a1=0;
a2=0;
e1=0;
e2=0;
end;

shocks;
var e1;
stderr .08;
var e2;
stderr .08;
end;

steady;

estimated_params;
rho, normal_pdf, .84,.05; 
alpha, normal_pdf, .38, .03;
stderr e1, inv_gamma_pdf, .078, inf;
stderr e2,inv_gamma_pdf,.082, inf;
end;

varobs c1 k2;

estimation(datafile=simuldata_international,mode_compute=6,mh_replic=1000,mh_jscale=.5);