// Basic OLG with population growth and technological growth

periods 15;

var c c1 c2 c3 c4 k l w R r s1 s2 s3 y pie tr g beq P1 P2 P3 P4 Ptot beta2 beta3 beta4 h2 h3 h4 check1; % check2

varexo tauw tauk tauc med2 med3 med4 h1;

parameters beta alpha delta a1 a2 a3 deltah;
beta = 0.98^(55/4); // 0.97^30
alpha = 1/3; delta = 1-(1-0.08)^(55/4); a1=.5; a2=1; a3=.5; deltah=.2;

model;

P1=1;
P2=beta2(-1)*P1(-1);
P3=beta3(-1)*P2(-1);
P4=beta4(-1)*P3(-1);
Ptot=P1+P2+P3+P4;

h2=1/(1-beta2);
h3=1/(1-beta3);
h4=1/(1-beta4);

h2=(1-deltah)*h1+med2;
h3=(1-deltah)*h2(-1)+med3;
h4=(1-deltah)*h3(-1)+med4;


(1/(1+tauc))*1/c1 = (1/(1+tauc(+1)))*beta*beta2*(1+(1-tauk(+1))*r(+1))/c2(+1);
(1/(1+tauc))*1/c2 = (1/(1+tauc(+1)))*beta*beta3*(1+(1-tauk(+1))*r(+1))/c3(+1);
(1/(1+tauc))*1/c3 = (1/(1+tauc(+1)))*beta*beta4*(1+(1-tauk(+1))*r(+1))/c4(+1);

(1+tauc)*c1 + s1 = (1-tauw)*a1*w + beq;
(1+tauc)*c2 + s2 = (1-tauw)*a2*w + (1+(1-tauk)*r)*s1(-1) + beq;
(1+tauc)*c3 + s3 = (1-tauw)*a3*w + (1+(1-tauk)*r)*s2(-1) + beq;
(1+tauc)*c4      = tr            + (1+(1-tauk)*r)*s3(-1) + beq;


Ptot*beq=(P1(-1)*(1-beta2(-1))*(1+(1-tauk)*r)*s1(-1)+P2(-1)*(1-beta3(-1))*(1+(1-tauk)*r)*s2(-1)+P3(-1)*(1-beta4(-1))*(1+(1-tauk)*r)*s3(-1));
tr=(P1*tauw*a1*w+P2*tauw*a2*w+P3*tauw*a3*w)/P4;

g=P1*tauc*c1+P2*tauc*c2+P3*tauc*c3+P4*tauc*c4+P1(-1)*tauk*r*s1(-1)+P2(-1)*tauk*r*s2(-1)+P3(-1)*tauk*r*s3(-1);

l =  P1*a1+P2*a2+P3*a3;
y = k(-1)^alpha*l^(1-alpha); 
w = (1-alpha)*y/l;
R = alpha*y/k(-1);
R = r + delta;

c = P1*c1 + P2*c2 + P3*c3 + P4*c4;
pie = y - w*l - R*k(-1);

k =  P1*s1 + P2*s2 + P3*s3;

check1 = y -(k - (1-delta)*k(-1) + c + g);

% check2 = k - s1 - s2- s3;
%y = k - (1-delta)*k(-1) + c + g;

end;

initval;
P1   	=	 1;
P2   	=	 .9;
P3   	=	 .8;
P4   	=	 .7;
Ptot 	=	 3;
c    	=	 0.698724;
c1   	=	 0.0836133;
c2   	=	 0.127202;
c3   	=	 0.193514;
c4   	=	 0.294395;
k    	=	 0.175085;
l    	=	 4;
w    	=	 0.296016;
R    	=	 1.6907;
r    	=	 1.12049;
s1   	=	 0.011631;
s2   	=	 0.0906495;
s3   	=	 0.0728049;
y    	=	 0.888048;
pie  	=	 0;
tr   	=	 0.17761;
g    	=	 0.0698724;
beq  	=	 0;
tauw   =          .3;
tauk   =          .1;
tauc   =          .1;
beta2 =.98067; // Survie à l'age de 25 ans en france en 1990 (source insee 2012-2016)
beta3 =.94927; // Survie à l'age de 45 ans en france en 1990 (source insee 2012-2016)
beta4=.81687;  // Survie à l'age de 65 ans en france en 1990 (source insee 2012-2016)
h1=2;
h2=2;
h3=2;
h4=2;
med2=1;
med3=1;
med4=1;
end;
resid;
steady;


endval;
med2 = 1.5;
%tauw   =          .2;
end;
resid;
steady;


%shocks;
%var xpop;
%periods 1:50;
%values .5;
%end;

perfect_foresight_setup;
perfect_foresight_solver;

run plotsim.m