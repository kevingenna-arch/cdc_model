// Basic OLG with population growth and technological growth


% periods 50;

var y c D k R rd r w g pen penb penw Tw Tc Tk tauc  c1 c2 c3 c4 s1 s2 s3 ;

varexo tauw tauk rhob rhow;

parameters beta alpha delta a1 a2 a3 l;
beta = 0.98^(55/4); // 0.97^30
alpha = 1/3; delta = 1-(1-0.08)^(55/4); a1=.8; a2=.9; a3=1; l=a1+a2+a3 ;

model;

(1/(1+tauc))*1/c1 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c2(+1);
(1/(1+tauc))*1/c2 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c3(+1);
(1/(1+tauc))*1/c3 = (1/(1+tauc(+1)))*beta*(1+(1-tauk(+1))*r(+1))/c4(+1);

(1+tauc)*c1 + s1 = (1-tauw)*a1*w;
(1+tauc)*c2 + s2 = (1+(1-tauk)*r)*s1(-1)+(1-tauw)*a2*w;
(1+tauc)*c3 + s3 = (1+(1-tauk)*r)*s2(-1)+(1-tauw)*a3*w;
(1+tauc)*c4      = (1+(1-tauk)*r)*s3(-1)+pen;

penb=((1-tauw)*a1*w+(1-tauw)*a2*w+(1-tauw)*a3*w)/3;
penw=((1-tauw)*a3*w(-1)+(1-tauw)*a2*w(-2)+(1-tauw)*a1*w(-3))/3;
pen=rhob*penb+rhow*penw;


Tw=tauw*a1*w+tauw*a2*w+tauw*a3*w;
Tc=tauc*c1+tauc*c2+tauc*c3+tauc*c4;
Tk=tauk*r*s1(-1)+tauk*r*s2(-1)+tauk*r*s3(-1);

g=.05*y;
D=.2*y;

//test=rd -(1-tauk)*r;

D-D(-1)=rd*D(-1)+g+pen-(Tw+Tc+Tk);


y = k(-1)^alpha*l^(1-alpha); 
w = (1-alpha) * (k(-1)^alpha) * (l^-alpha);         //wage foc from firm problem
R =  alpha * k(-1)^(alpha-1) * l^(1-alpha);         //interest rate foc from firm problem
R = (1-tauk)*r + delta;

c = c1 + c2 + c3 + c4;
k + D = s1 + s2 + s3;

y = k - (1-delta)*k(-1) + c + g;

end;

initval;
D   =    .5;
c  	=	 0.578033;
c1 	=	 0.184384;
c2 	=	 0.192557;
c3 	=	 0.201092;
c4 	=	 0.201092;
k  	=	 0.0874097;
w  	=	 0.243223;
r  	=	.01;
s1 	=	 0.0101944;
s2 	=	 0.0772153;
s3 	=	 0.0772153;
g   =    .2;
y  	=	 0.656702;
penw=.3;
tauk   =          .1;
tauc=   .1;
tauw=.1;
rhow=.5;
rhob=.5;
end;
resid;
steady;
check;



% endval;
% tauw=.2;
% %rhow=.7;
% end;
resid;
steady;


%shocks;
%var tauw;
%periods 1:8;
%values .1;
%end;

 
perfect_foresight_setup(periods = 50);

%perfect_foresight_solver(stack_solve_algo=6);
perfect_foresight_solver;

% run plotsim.m

