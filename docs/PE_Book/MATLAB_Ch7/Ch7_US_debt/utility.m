function [ y ] = utility(c,n)
% lifetime utility 

    def_global_USdebt
	y=0;
	for i=1:1:nw
		y = y+beta1^(i-1)*prod(sp(1:i))*instantutil(c(i),n(i));
    end

    for i=nw+1:1:nw+nr
		y = y+beta1^(i-1)*prod(sp(1:i))*instantutil(c(i),0);
    end
		

end

