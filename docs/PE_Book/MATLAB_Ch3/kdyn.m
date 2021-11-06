function [kdiff] = kdyn(par,x)
% Dynamics of the capital stock
% Equation (.3.38) in Chapter 3.3.4 in Heer (2018)
	
    alpha1=par.alpha;
    n=par.n;
    
    kt2=zeros(par.bigt+2,1);
    kt2(1)=par.k0;
	kt2(2)=x;
	for i=3:par.bigt+2,
		kt2(i)=kt2(i-2)+kt2(i-2)^(alpha1) -(1+n)*kt2(i-1) ;
        kt2(i)=(kt2(i-1)+kt2(i-1)^(alpha1))/(1+n) - (1+alpha1*kt2(i-1)^(alpha1-1))/((1+n)^2) * kt2(i);
    end;
	kdiff=kt2(par.bigt+2)-par.kfinal; 

kdiff

end
