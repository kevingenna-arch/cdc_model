function [xmax] = golden(f,ay,by,cy,tol)
% golden section search
%  computes the maximum (!) of a one-dimensional function
%  Inputs:  f  function
%           ay,by,cy grid points with ay<by<cy
%  Output:      maximum
%  Remark:      see Press et al., chapter 10.1 
    r1=0.61803399; 
    r2=1-r1;
    x0=ay;
    x3=cy;
    if abs(cy-by)<=abs(by-ay)
        x1=by; 
        x2=by+r2*(cy-by);
    else
        x2=by; 
        x1=by-r2*(by-ay);
    end
    f1=-f(x1);
    f2=-f(x2);
    while abs(x3-x0)>=tol*max(1,(abs(x1)+abs(x2)))
        if f2<f1
            x0=x1;
            x1=x2;
            x2=r1*x1+r2*x3;
            f1=f2;
            f2=-f(x2);
        else
            x3=x2;
            x2=x1;
            x1=r1*x2+r2*x0;
            f2=f1;
            f1=-f(x1);
        end
    end
    if f1<=f2
        xmax=x1;
    else
        xmax=x2;
    end    
end


