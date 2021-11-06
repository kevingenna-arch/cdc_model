// Basic OLG with population growth and technological growth


periods 25;

var c c1 c2 c3 c4 k w R r s1 s2 s3 y pie pen g tauw; % check1 check2

varexo tauk tauc rhop;

parameters beta alpha delta a1 a2 a3 l;
beta = 0.98^(55/4); // 0.97^30
alpha = 1/3; delta = 1-(1-0.08)^(55/4); a1=1; a2=1; a3=.5; l=a1+a2+a3;

model;

(1/(1+tauc))*1/c1 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c2(+1);
(1/(1+tauc))*1/c2 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c3(+1);
(1/(1+tauc))*1/c3 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c4(+1);

(1+tauc)*c1 + s1 = (1-tauw)*a1*w;
(1+tauc)*c2 + s2 = (1+(1-tauk)*r)*s1(-1)+(1-tauw)*a2*w;
(1+tauc)*c3 + s3 = (1+(1-tauk)*r)*s2(-1)+(1-tauw)*a3*w;
(1+tauc)*c4      = (1+(1-tauk)*r)*s3(-1)+pen;


pen=rhop*(a3*w(-1)+a2*w(-2)+a1*w(-3))/3;

pen=(tauw*a1*w+tauw*a2*w+tauw*a3*w);

g=tauc*c1+tauc*c2+tauc*c3+tauc*c4;

y = k(-1)^alpha*l^(1-alpha); 
w = (1-alpha)*y/l;
R = alpha*y/k(-1);
R = (1-tauk)*r +delta;

c = c1 + c2 + c3 + c4;
pie = y - w*l - R*k(-1);

y = k - (1-delta)*k(-1) + c + g;


% check1 = y -(k - (1-delta)*k(-1) + c);
% check2 = k - s1 - s2;

end;

initval;
c  	=	 0.578033;
c1 	=	 0.184384;
c2 	=	 0.192557;
c3 	=	 0.201092;
c4 	=	 0.201092;
k  	=	 0.0874097;
w  	=	 0.243223;
r  	=	 1.6043;
s1 	=	 0.0101944;
s2 	=	 0.0772153;
s3 	=	 0.0772153;
g   =    .1;
y  	=	 0.656702;
tauk   =          .1;
tauc=.1;
rhop=.7;
end;
resid;
steady;


endval;
tauk   =          .1;
tauc=     .1;
rhop=.6;
end;
resid;
steady;


%shocks;
%var tauw;
%periods 8;
%values .1;
%end;

%simul(periods=25);
 
perfect_foresight_setup;

%perfect_foresight_solver(stack_solve_algo=6);
perfect_foresight_solver;

run plotsim.m