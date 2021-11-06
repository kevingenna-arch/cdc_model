function [ Jac ] = jac_cd( fun,x0 )
%Calculates the Jacobian of a function 'fun' (using centred differences)
% INPUTS: 
%   fun - string of function name with output that is column vector (m)
%   x0  - row vector of the point where function is to be evaluated (n)
%
% OUTPUT:
%  Jac  - m x n Jacobian matrix

% Sijmen Duineveld, 2-10-2015;
% code based on Charles Bos' code for gradient, who adapted from Doornik:
% http://www.matrixlab-examples.com/gradient.html).

if size(x0,1) ~= 1;
    error('Input x0 is not row vector')
end
n = size(x0,2);

y0  = feval(fun,x0);
if size(y0,2) ~= 1;
    error('Output y from function is not column vector')
end
yl  = size(y0,1);

Jac = NaN(yl,n);
for ix = 1:n;
    dx      = zeros(1,n);
    dx(ix)  = max(1e-8,abs(x0(ix)/1e8));
    
    x1 = x0 + dx;    
    y1  = feval(fun,x1);
    
    x2 = x0 - dx;    
    y2  = feval(fun,x2);
    
    Jac(:,ix)  = (y1-y2)./(x1(ix)-x2(ix));  
    
end

end

