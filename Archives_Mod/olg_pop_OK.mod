// Basic OLG with population growth and technological growth

% periods 25;

var test P1 P2 P3 P4 Ptot y k l w R r c tr g beq c1 c2 c3 c4 s1 s2 s3;

varexo tauw tauk tauc beta2 beta3 beta4 mig2 mig3 mig4 xpop;

parameters beta alpha delta a1 a2 a3;
beta = 0.97; // 0.97^30 0.98^(55/4)
alpha = 1/3; 
delta = .02;  //1-(1-0.08)^(55/4)
delta = .02;  //1-(1-0.08)^(55/4)
a1=.8; 
a2=.9; 
a3=1;

model;

% pop dynamics 
P1=xpop;
P2=beta2(-1)*P1(-1)+mig2;
P3=beta3(-1)*P2(-1)+mig3;
P4=beta4(-1)*P3(-1)+mig4;
Ptot=P1+P2+P3+P4;

(1/(1+tauc))*1/c1 = (1/(1+tauc(+1)))*beta*beta2*(1+(1-tauk(+1))*r(+1))/c2(+1);
(1/(1+tauc))*1/c2 = (1/(1+tauc(+1)))*beta*beta3*(1+(1-tauk(+1))*r(+1))/c3(+1);
(1/(1+tauc))*1/c3 = (1/(1+tauc(+1)))*beta*beta4*(1+(1-tauk(+1))*r(+1))/c4(+1);

(1+tauc)*c1 + s1 = (1-tauw)*a1*w + beq;
(1+tauc)*c2 + s2 = (1-tauw)*a2*w + (1+(1-tauk)*r)*s1(-1) + beq;
(1+tauc)*c3 + s3 = (1-tauw)*a3*w + (1+(1-tauk)*r)*s2(-1) + beq;
(1+tauc)*c4      = tr            + (1+(1-tauk)*r)*s3(-1) + beq;


y = k(-1)^alpha*l^(1-alpha); 
w = (1-alpha) * (k(-1)^alpha) * (l^-alpha);         //wage foc from firm problem
R =  alpha * k(-1)^(alpha-1) * l^(1-alpha);  //interest rate foc from firm problem
R = (1-tauk)*r + delta;

test= l - (P1*a1 + P2*a2 + P3*a3);

// Commenter cette equation si on ajoute l'emploi aggrégé.
beq=((1-beta2(-1))*P1(-1)*(1+(1-tauk)*r)*s1(-1)+(1-beta3(-1))*P2(-1)*(1+(1-tauk)*r)*s2(-1)+(1-beta4(-1))*P3(-1)*(1+(1-tauk)*r)*s3(-1))/Ptot;

g=P1*tauc*c1+P2*tauc*c2+P3*tauc*c3+P4*tauc*c4+P2*tauk*r*s1(-1)+P3*tauk*r*s2(-1)+P4*tauk*r*s3(-1);

tr=(tauw*P1*a1*w+tauw*P2*a2*w+tauw*P3*a3*w)/P4;


k = P1*s1 + P2*s2 + P3*s3;
c = P1*c1 + P2*c2 + P3*c3 + P4*c4;
//l = P1*a1 + P2*a2 + P3*a3;

y = k - (1-delta)*k(-1) + c + g;

end;

initval;
xpop=.1;
P1   	=	 .1;
P2   	=	 .1;
P3   	=	 .1;
P4   	=	 .1;
Ptot 	=	 .4;
c    	=	 0.698724;
c1   	=	 0.0836133;
c2   	=	 0.127202;
c3   	=	 0.193514;
c4   	=	 0.294395;
k    	=	 0.175085;
l    	=	 .04;
w    	=	 0.296016;
R    	=	 1.6907;
r    	=	 1.12049;
s1   	=	 0.011631;
s2   	=	 0.0906495;
s3   	=	 0.0728049;
y    	=	 0.888048;
tr   	=	 0.17761;
g    	=	 0.0698724;
beq  	=	 0.05;
tauw   =          .02;
tauk   =          0;
tauc   =          .02;
beta2 =.999;
beta3 =.999;
beta4=.999;
mig2=0;
mig3=0;
mig4=0;
test=0;
end;

steady;
resid;
check;



endval;
beta2 =.98;
tauw   =          .2;
%beta2 =.9;
%mig2=.1;
end;

steady;
resid;


perfect_foresight_setup(periods = 250);
perfect_foresight_solver;


%run plotsim.m


