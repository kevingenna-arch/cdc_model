function [ y ] = rvalue(a)
% interpolation of the retired agent value
% for the next-period asset grid point
% for given ix,ieta,ieps

    def_global
%    y = interp1(agrid,vr11,a,VI_method,'extrap');
% does linear interpolation with interp1 take too long?

    if a==assetmin 
        y= vr11(1); 
    
    elseif a>=assetmax  
        y=vr11(na);
    else
        ia0=sum(a>agrid);
        lambda=(agrid(ia0+1)-a) / (agrid(ia0+1)-agrid(ia0) );

		if ia0==na
			y=vr11(ia0);
        else
			y=lambda*vr11(ia0)+(1-lambda)*vr11(ia0+1);
        end
    end
		
end

