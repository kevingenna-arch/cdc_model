% Estimates the Hansen Sargent and Tallarini model by maximum likelihood.

var s c h k i d dhat dbar mus muc muh gamma R;
varexo e_dhat e_dbar;

parameters lambda deltah deltak mud b bet phi1 phi2 cdbar alpha1 alpha2 cdhat;
bet=0.9971;
deltah=0.682;
lambda=2.443;
alpha1=0.813;
alpha2=0.189;
phi1=0.998;
phi2=0.704;
mud=13.710;
cdhat=0.155;
cdbar=0.108;
b=32;
deltak=0.975;

model(linear);
R=deltak+gamma;
R*bet=1;
s=(1+lambda)*c-lambda*h(-1);
h=deltah*h(-1)+(1-deltah)*c;
k=deltak*k(-1)+i;
c+i=(1/bet-deltak)*k(-1)+d;
mus=b-s;
muc=(1+lambda)*mus+(1-deltah)*muh;
muh=bet*(deltah*muh(+1)-lambda*mus(+1));
muc=muc(+1);
d=mud+dbar+dhat;
dbar=(phi1+phi2)*dbar(-1) - phi1*phi2*dbar(-2) + cdbar*e_dbar;
dhat=(alpha1+alpha2)*dhat(-1) - alpha1*alpha2*dhat(-2) + cdhat*e_dhat;
end;

initval;
R=1/bet;
gamma=R-deltak;
d=mud;
k=100;
i=(1-deltak)*k;
c=gamma*k+d-i;
h=c;
s=c;
mus=b-s;
muh=lambda/(deltah-1/bet)*mus;
muc=(1+lambda)*mus+(1-deltah)*muh;
end;

steady(nocheck);

shocks;
var e_dhat;
stderr 1;
var e_dbar;
stderr 1;
end;

% stoch_simul(irf=0, periods=500);
% save dataHST c i;

estimated_params;
bet,normal_pdf, .9971, 0.05;
deltah,normal_pdf, 0.682,0.05;
lambda,normal_pdf, 2.443, 0.05;
alpha1,beta_pdf, 0.813, 0.2;
alpha2,beta_pdf, 0.189, 0.2;
phi1,beta_pdf, 0.998, 0.01;
phi2,beta_pdf, 0.704, 0.01;
mud,normal_pdf, 13.71, 0.05;
cdhat,normal_pdf, 0.155,0.05;
cdbar,normal_pdf, 0.108, 0.05;

end;

varobs c i;
% estimation(datafile=dataHST,first_obs=1,nobs=500,mode_compute=4,MH_jscale=2);
estimation(datafile=dataHST,first_obs=1,nobs=500,mode_compute=4,mode_check,mh_replic=5000,mh_nblocks=1,mh_jscale=0.3,diffuse_filter);
