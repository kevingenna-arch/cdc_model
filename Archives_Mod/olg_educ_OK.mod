// Basic OLG with schooling


periods 50;

var tauw pil pih y l c r k R w AIML AIMH lambc g c11 c21 c31 s11 s21 c12 c22 c32 s12 s22 pie pen1 pen2 check1 check2;

varexo tauk tauc;

parameters beta alpha delta a1 a2 h1 h2 v lambb lambbb; 
beta = 0.98^(55/4); // 0.97^30
alpha = 1/3; delta = 1-(1-0.08)^(55/4); a1=1; a2=1; h1=.8; h2=1.5; v=.2; lambb=50; lambbb=-50;

model;

% Low-skilled
(1/(1+tauc))*1/c11 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c21(+1);
(1/(1+tauc))*1/c21 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c31(+1);

(1+tauc)*c11 + s11 = (1-tauw)*a1*h1*w;
(1+tauc)*c21 + s21 = (1-tauw)*a2*h1*w + (1+(1-tauk)*r)*s11(-1);
(1+tauc)*c31 = (1+(1-tauk)*r)*s21(-1)+pen1;

% High-skilled
(1/(1+tauc))*1/c12 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c22(+1);
(1/(1+tauc))*1/c22 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c32(+1);

(1+tauc)*c12 + s12 = (1-tauw)*a1*h2*(1-v)*w;
(1+tauc)*c22 + s22 = (1-tauw)*a2*h2*w + (1+(1-tauk)*r)*s12(-1);
(1+tauc)*c32 = (1+(1-tauk)*r)*s22(-1)+pen2;


AIML=(1-tauc)*a1*h1*w+((1-tauc(+1))*a2*h1*w(+1))/(1+(1-tauk(+1))*r(+1))+pen1(+2)/((1+(1-tauk(+1))*r(+1))*(1+(1-tauk(+2))*r(+2)));
AIMH=(1-tauc)*a1*h2*(1-v)*w+((1-tauc(+1))*a2*h2*w(+1))/(1+(1-tauk(+1))*r(+1))+pen2(+2)/((1+(1-tauk(+1))*r(+1))*(1+(1-tauk(+2))*r(+2)));


pil=min(1,max(0,(lambb-lambc)/(lambb-lambbb)));
%pih=min(1,max(0,(lambc-lambbb)/(lambb-lambbb)));
pih=1-pil;

%pil=.7;
%pih=.3;

lambc=(AIMH-AIML)/(a1*h1*w*v);

pen1=.7*(pil(-2)*a1*h1*w(-2)+pil(-2)*a2*h1*w(-1));
pen2=.7*(pih(-2)*a1*h2*(1-v)*w(-2)+pih(-2)*a2*h2*w(-1));

pil(-2)*pen1+pih(-2)*pen2=tauw*(pil*a1*h1*w+pil(-1)*a2*h1*w+pih*a1*h2*(1-v)*w+pih(-1)*a2*h2*w);

g=tauc*(pil*c11+pil(-1)*c21+pil(-2)*c31+pih*c12+pih(-1)*c22+pih(-2)*c32);


l=h1*(pil*a1*h1+pil(-1)*a2*h1)+h2*(pih*a1*h2*(1-v)+pih(-1)*a2*h2);

y = k(-1)^alpha*l^(1-alpha); 
w = (1-alpha)*y/l;
R = alpha*y/k(-1);
R = (1-tauk)*r +delta;

c = pil*c11 + pil(-1)*c21 + pil(-2)*c31 + pih*c12 + pih(-1)*c22 + pih(-2)*c32;
pie = y - w*l - R*k(-1);

y = k - (1-delta)*k(-1) + c + g;

check1 = y -(k - (1-delta)*k(-1) + c + g);
check2 = k - (pil*s11 + pil(-1)*s21) + (pih*s12 + pih(-1)*s22);

end;

initval;
c  	=	 0.578033;
c11 	=	 0.184384;
c21 	=	 0.192557;
c31 	=	 0.201092;
c12 	=	 0.184384;
c22 	=	 0.192557;
c32 	=	 0.201092;
k  	=	 1;
w  	=	 .5;
r  	=	 1.5;
l=.5;
pen1=.5;
pen2=.5;
pie=0;
s11 	=	 0.0101944;
s21 	=	 0.0772153;
s12 	=	 0.0101944;
s22 	=	 0.0772153;
AIML=.5;
AIMH=.8;
%lambc   =   .5;
pil=.5;
pih=.5;
g   =    .1;
y  	=	 1;
tauw= .1;
tauk   = .1;
tauc=.1;
end;
resid;
steady;


endval;
tauc=.2;
end;
resid;
steady;


%shocks;
%var tauw;
%periods 8;
%values .1;
%end;

 
perfect_foresight_setup;

perfect_foresight_solver;

run plotsim.m