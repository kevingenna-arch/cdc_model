// Basic OLG with schooling


periods 20;

var pim pih tauw tauc y l c r k R w WEL_M WEL_H lambc g c11 c21 c31 c41 s11 s21 s31 c12 c22 c32 c42 s12 s22 s32 c13 c23 c33 c43 s13 s23 s33 pen1 pen2 pen3;

varexo A tauk coefg coefed pil;

parameters beta alpha delta a1 a2 a3 h1 h2 h3 v lambb lambbb deped; 
beta = 0.98^(55/4); // 0.97^30
alpha = 1/3; delta = 1-(1-0.08)^(55/4); a1=1; a2=1.5; a3=1.6; h1=.5; h2=1; h3=3; v=.5; lambb=50; lambbb=-50; deped=.2;

model;

% Low-skilled
(1/(1+tauc))*1/c11 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c21(+1);
(1/(1+tauc))*1/c21 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c31(+1);
(1/(1+tauc))*1/c31 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c41(+1);

(1+tauc)*c11 + s11 = (1-tauw)*a1*h1*w;
(1+tauc)*c21 + s21 = (1-tauw)*a2*h1*w + (1+(1-tauk)*r)*s11(-1);
(1+tauc)*c31 + s31 = (1-tauw)*a3*h1*w + (1+(1-tauk)*r)*s21(-1);
(1+tauc)*c41 = (1+(1-tauk)*r)*s31(-1)+pen1;

% Middle-skilled
(1/(1+tauc))*1/c12 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c22(+1);
(1/(1+tauc))*1/c22 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c32(+1);
(1/(1+tauc))*1/c32 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c42(+1);

(1+tauc)*c12 + s12 + deped*coefed = (1-tauw)*a1*h2*(1-v)*w;
(1+tauc)*c22 + s22 = (1-tauw)*a2*h2*w + (1+(1-tauk)*r)*s12(-1);
(1+tauc)*c32 + s32 = (1-tauw)*a3*h2*w + (1+(1-tauk)*r)*s32(-1);
(1+tauc)*c42 = (1+(1-tauk)*r)*s32(-1)+pen2;

% High-skilled
(1/(1+tauc))*1/c13 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c23(+1);
(1/(1+tauc))*1/c23 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c33(+1);
(1/(1+tauc))*1/c33 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c43(+1);

(1+tauc)*c13 + s13 + deped*coefed = (1-tauw)*a1*h3*(1-v)*w;
(1+tauc)*c23 + s23 + deped*coefed = (1-tauw)*a2*h3*(1-v)*w + (1+(1-tauk)*r)*s13(-1);
(1+tauc)*c33 + s33 = (1-tauw)*a3*h3*w + (1+(1-tauk)*r)*s23(-1);
(1+tauc)*c43 = (1+(1-tauk)*r)*s33(-1)+pen3;




WEL_M=log(c12)-lambc+beta*log(c22(+1))+beta^2*log(c32(+2))+beta^3*log(c42(+3));
WEL_H=log(c13)-lambc+beta*log(c23(+1))-beta*lambc+beta^2*log(c33(+2))+beta^3*log(c43(+3));


pim=(1-pil)*min(1,max(0,(lambb-lambc)/(lambb-lambbb)));
%pil=min(1,max(0,(lambc-lambbb)/(lambb-lambbb)));
pih=(1-pim-pil);


0=WEL_H-WEL_M;

pen1=.7*(a1*h1*w(-3)+a2*h1*w(-2)+a3*h1*w(-1))/3;
pen2=.7*(a1*h2*(1-v)*w(-3)+a2*h2*w(-2)+a3*h2*w(-1))/3;
pen3=.7*(a1*h3*(1-v)*w(-3)+a2*h3*(1-v)*w(-2)+a3*h3*w(-1))/3;

pil(-3)*pen1+pim(-3)*pen2+pih(-3)*pen3=tauw*(pim*a1*h2*(1-v)*w+pim(-1)*a2*h2*w+pim(-2)*a3*h2*w + pih*a1*h3*(1-v)*w+pih(-1)*a2*h3*(1-v)*w+pih(-2)*a3*h3*w+ pil*a1*h1*w+pil(-1)*a2*h1*w+pil(-2)*a3*h1*w);

g+pih*deped*(1-coefed)+pih(-1)*deped*(1-coefed)+pim*deped*(1-coefed)=tauc*(pim*c12+pim(-1)*c22+pim(-2)*c32+pim(-3)*c42+pih*c13+pih(-1)*c23+pih(-2)*c33+pih(-3)*c43+pil*c11+pil(-1)*c21+pil(-2)*c31+pil(-3)*c41);

g=coefg*y;

l=h1*(pil*a1*h1+pil(-1)*a2*h1+pil(-2)*a3*h1)+h2*(pim*a1*h2*(1-v)+pim(-1)*a2*h2+pim(-2)*a3*h2)+h3*(pih*a1*h3*(1-v)+pih(-1)*a2*h3*(1-v)+pih(-2)*a3*h3);

y = A*k(-1)^alpha*l^(1-alpha); 
w = A*(1-alpha) * (k(-1)^alpha) * (l^-alpha);         //wage foc from firm problem
R =  A* alpha * k(-1)^(alpha-1) * l^(1-alpha);  //interest rate foc from firm problem
R = (1-tauk)*r + delta;

c = pil*c11 + pil(-1)*c21 + pil(-2)*c31 + pil(-3)*c41 + pim*c12 + pim(-1)*c22 + pim(-2)*c32 + pim(-3)*c42 +pih*c13 + pih(-1)*c23 + pih(-2)*c33+ pih(-3)*c43;

y = k - (1-delta)*k(-1) + c + g + pih*deped*(1-coefed) + pih(-1)*deped*(1-coefed) + pim*deped*(1-coefed);

//check1 = y -(k - (1-delta)*k(-1) + c + g);
//check2 = k - (pil*s11 + pil(-1)*s21) + (pim*s12 + pim(-1)*s22) + (pih*s13 + pih(-1)*s23);

end;

initval;
c  	=	 0.578033;
c11 	=	 0.184384;
c21 	=	 0.192557;
c31 	=	 0.201092;
c41 	=	 0.201092;

c12 	=	 0.184384;
c22 	=	 0.192557;
c32 	=	 0.201092;
c42 	=	 0.201092;
c13 	=	 0.184384;
c23 	=	 0.192557;
c33 	=	 0.201092;
c43 	=	 0.201092;
k  	=	 1;
w  	=	 .5;
r  	=	 1.5;
l=.5;
pen1=.4;
pen2=.5;
pen3=.6;
s11 	=	 0.0101944;
s21 	=	 0.0772153;
s12 	=	 0.0101944;
s22 	=	 0.0772153;
s13 	=	 0.0101944;
s23 	=	 0.0772153;
WEL_M=.5;
WEL_H=.8;
pil=1/3;
pim=.5;
pih=.5;
g   =    .2;
y  	=	 1;
tauw= .1;
tauk   = .1;
tauc=.1;
coefg=.03;
coefed=.2;
A=1;
end;
resid;
steady;
check;



endval;
coefed=.5;
%A=1.1;
end;
resid;
steady;


%shocks;
%var h2;
%periods 1 2 3 4 5 6 7 8 9 10;
%values 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1;
%end;

 
perfect_foresight_setup;
perfect_foresight_solver;

run plotsim.m

