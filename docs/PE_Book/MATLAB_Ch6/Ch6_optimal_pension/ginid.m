function [ y ] = ginid(x,g)
% computes the gini for distribution where x has measure g(x) 
% input: x - column vector
%        g - frequency of x
    ng=size(x,1);
    if size(x,2)>1
        disp('wrong size of input x into Gini - must be column vector');
        pause;
    end
    
    if size(g,2)>1
        disp('wrong size of input g into Gini - must be column vector');
        pause;
    end
        
    temp = [x zeros(ng,1)];
    x=max(temp,[],2);
	xmean=x'*g;
    y=[x g];
    y=sort(y,1);
    x=y(:,1);
    g=y(:,2);
    f=zeros(ng,1);      % accumulated frequency 
    f(1)=g(1)*x(1)/xmean;
	gini=1-f(1)*g(1);
    for i=2:1:ng
        f(i)=f(i-1)+g(i)*x(i)/xmean;
        gini=gini-(f(i)+f(i-1))*g(i);
    end

    y= gini;

end

