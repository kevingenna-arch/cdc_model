function vc = VCX(x,maxlag)
%Compute variance covariance matrix for the time series in the matrix x for lags 0, ..., maxlag

%{
    Author: Alfred Mauﬂner

    First version: 04 February 2017

    Purpose: compute the covariance matrix between the columns of the matrix x
             a lags 0 through maxlag.

    Input:    x   : nobs by nvar matrix
            maxlag: scalar

    Output:  vc nvar by nvar by maxlag+1 three dimensional array. vc(i,j,k) store
             the covariance between i an j at lag k. The respective formula
             (taken from the EViews user guide) is

             cv(i,j,k)=x(k:nobs,i)'*x(1:nobs-k,j)/nobs-xi*xj
 
             where xi and xi are the means of colum i and j respectively.

%}


[nobs,nvar]=size(x);
xbar=mean(x,1);
vc=zeros(nvar,nvar,maxlag+1);

for k=1:maxlag+1; % loop over lags starts here
    for i=1:nvar; % loop over variable i starts here
        for j=1:nvar; % loop over variable j starts here
            vc(i,j,k)=x(k:nobs,i)'*x(1:nobs-k+1,j)/nobs - xbar(i)*xbar(j);
        end;
    end;
end;
    
return;

end

