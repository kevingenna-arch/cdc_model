// Basic OLG with population growth and technological growth
var c c1 c2 k w r s y check1 check2;
varexo x;

parameters beta alpha delta n; % x;
beta = 0.4010; // 0.97^30
alpha = 1/3; 
delta = 0.9; 
n=0.1614; 
% x=0.3478;

kss = ((1+1/beta)*(1+n)/(1-alpha))^(1/(alpha-1));
wss = (1-alpha)*kss^alpha; 
rss = alpha*kss^(alpha-1)-delta;
sss = (1+n)*kss; 
c1ss = sss/beta; 
c2ss = (1+rss)*sss;
yss = kss^alpha; 
css = c1ss+c2ss/(1+n);

model;
// In effective labour form
1/c1 = beta*(1+r(+1))/c2(+1);
c1 + s = w;
c2(+1) = (1+r(+1))*s;
y = k(-1)^alpha; // Y/XL = (K/XL)^alpha
w = (1-alpha)*k(-1)^alpha;
r = alpha*k(-1)^(alpha-1) - delta;
c = c1 + c2/(1+n+x); // C = XL*c1 + X(-1)*L(-1)*C2;
// ECNOMY WIDE FEASIBILITY CONSTRAINT
y = (1+n+x)*k - (1-delta)*k(-1) + c;
check1 = y -( (1+n+x)*k - (1-delta)*k(-1) + c);
check2 = (1+n+x)*k - s;
end;

initval;
x 	=	.1;
k 	= 	kss; 
c1  = 	c1ss; 
c2 	= 	c2ss; 
r 	=	rss; 
w 	=	wss; 
s 	=	sss; 
c 	=	css; 
y 	= 	yss;
end;

shocks;
var x;
periods 50:150;
values .3;
end;

steady;


perfect_foresight_setup(periods=500);
perfect_foresight_solver(maxit = 10);