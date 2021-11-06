function [ y ] = getss(x)
% auxiliary function to compute steady state

    [y1,y2,y3,y4,y5,y6]=getvaluess(x(1),x(2),x(3),x(4),x(5),x(6));
    y(1) = y1-x(1);
    y(2) = y2-x(2);
    y(3) = y3-x(3);
    y(4) = y4-x(4);
	y(5) = y5-x(5);
	y(6) = y6-x(6);
    
end

