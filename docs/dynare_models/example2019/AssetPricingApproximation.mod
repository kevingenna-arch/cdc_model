
var  dc, dd, v_c, v_d, x;
varexo e_c, e_x, e_d;

parameters DELTA THETA PSI MU_C MU_D RHO_X LAMBDA_DX;

DELTA=.99;
%PSI=0.133;
%THETA=1;
PSI=1.5;
THETA=(1-7.5)/(1-1/PSI); 
MU_C=0.0015;
MU_D=0.0015;
RHO_X=.979;
LAMBDA_DX=3;

model;
v_c       = DELTA^THETA * exp((-THETA/PSI)*dc(+1) + (THETA-1)*log((1+v_c(+1))*exp(dc(+1))/v_c) ) * (1+v_c(+1))*exp(dc(+1));
v_d       = DELTA^THETA * exp((-THETA/PSI)*dc(+1) + (THETA-1)*log((1+v_c(+1))*exp(dc(+1))/v_c) ) * (1+v_d(+1))*exp(dd(+1));
dc        = MU_C  + x(-1) + e_c;
dd        = MU_D + LAMBDA_DX*x(-1) + e_d;
x         = RHO_X * x(-1) + e_x;
end;

initval;
dc=MU_C;
dd=MU_D;
v_c=DELTA*exp((1-1/PSI)*dc)/(1-DELTA*exp((1-1/PSI)*dc));
v_d=(DELTA^THETA*exp(-THETA/PSI*dc+dc*(THETA-1)+dd)*(1+v_c)^(THETA-1))^(1/THETA);
x=0;
e_c=0;
e_x=0;
e_d=0;
end;

shocks;
var e_c;
stderr .0078;
var e_x;
stderr .0078*.044;
var e_d;
stderr .0078*4.5;
end;

steady;
check;

stoch_simul(order=1, periods=1000, irf=30);

save AssetPricing_simuldata v_d dd dc;