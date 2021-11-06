function [ lopt0, copt0, aopt0, val0 ] = getoptpolicy(wseq1,rseq1,penseq1,trseq1,tauwseq1,taubseq1)
%   gets the optimal policy for given wages, pensions, transfers and taxes
%   solves the individual optimization problem with value function
%   iteration

%	local val0, aopt0, copt0, lopt0, v0;	
%	local c,y,cmax, m0, m, vr, df, ax, bx, cx, v1, ncol, ieta1, wu0, vw01, vw1, vw0, w0, l0; 
%	local a1, p, ve1, ve2;
	
    def_global

	rseq=rseq1;		% rseq1,.. need to be stored as global variables so that
	wseq=wseq1;		% can be used in value1 and wvalue
	penseq=penseq1;
	tauwseq=tauwseq1;
	taubseq=taubseq1;
	trseq=trseq1;
	
	val0=zeros(nage,neps,neta,na);
	aopt0=zeros(nage,neps,neta,na);
	copt0=zeros(nage,neps,neta,na);
	lopt0=zeros(t,neps,neta,na);	% labor supply worker

	%
	% Part 1.: Old's problem
	%

	% value of the value function in the last period of lifetime
    
    for ia=1:1:na
        for ieta=1:1:neta
            for ieps=1:1:neps
				c=(1+(1-taur)*rseq(nage))*agrid(ia)+penseq(nage,ieps)+trseq(nage);		
				copt0(nage,ieta,ieps,ia)=c;		% all theta/epsilon types
				y=u(c,0);
				if c>0
					val0(nage,ieta,ieps,ia)=y;	% behave alike in old age
                else
					val0(nage,ieta,ieps,ia)=neg;
                end
            end
        end
    end
		

	% solving the individual value function problem
	% in old age, epsilon=epsilon1, eta=eta1

	it=tr;
	for it=tr-1:-1:1        % all periods t=T+1,T+2,..T+TR 
		for ia=1:1:na        % asset holding in period t
            for ieta = 1:1:neta
                for ieps = 1:1:neps
                    


                    % maximum consumption at age it+t
					cmax=(1+(1-taur)*rseq(it+t))*agrid(ia)+penseq(it+t,ieps)+trseq(it+t);
					
					if cmax>0
						m0=0;			% computational parameter, monotonocity
		
				% extracting the matrix of the value function array
				% that contains the value for the retired agent
				% that is one year older and has individual
				% types ieta and ieps
			
                    	vr1=val0(t+it+1,ieta,ieps,:);		% CHECK
						vr11(1:na)=vr1;			% check if this is a vector
                        vr11=vr11';              % column vector
                        
				% check monotonocity of value function
						df=vr11(2:na)-vr11(1:na-1);
						if sum(df<=0)>0
							disp('value function not monotone');
							pause;
                        end

				% golden section search to find optimum value of a
				% triple ax,bx,cx 			
						ax=0; bx=-1; cx=-2;
						v0=neg;
						m=max(0,m0-2);	 % monotonocity of the value function
						while (ax>bx) || (bx>cx);
							m=m+1;
							v1=value1(agrid(m));
							if v1>v0
                                if m==1 
                                    ax=agrid(m); 
                                    bx=agrid(m);
                                else
                                    bx=agrid(m); 
                                    ax=agrid(m-1);
                                end
                                v0=v1;
                                m0=m;   % monotonocity of the value function 
                            else
                                cx=agrid(m);
                            end
                            if m==na 
                                ax=agrid(m-1); 
                                bx=agrid(m); 
                                cx=agrid(m); 
                            end
                        end

                        if ax==bx       % corner solution: a'=0?
                            v1=value1(ax+psi1);
                            if v1>=v0	
                                bx=ax+psi1;
                                a1=golden(fhandle_function1,ax,bx,cx,tol_golden);
                            else
                                a1=ax;
                            end
				
                        elseif bx==cx;	% corner solution: a'=assetmax?
                            v0=value1(agrid(na));
                            v1=value1(agrid(na)-psi1);
                            if v0<=v1
                                bx=agrid(na)-psi1;
                                a1=golden(fhandle_function1,ax,bx,cx,tol_golden);
                            else
                                a1=agrid(na);
                            end
                        else		
                            a1=golden(fhandle_function1,ax,bx,cx,tol_golden);
                        end
						c=(1+(1-taur)*rseq(it+t))*agrid(ia)+penseq(it+t,ieps)+trseq(it+t)-a1*(1+garate);
					
                        % all iempl/eta/epsi1lon types behave alike in old age
                        aopt0(it+t,ieta,ieps,ia)=a1;
                        copt0(it+t,ieta,ieps,ia)=c;	
                        val0(it+t,ieta,ieps,ia)=value1(a1);
                    else 	% cmax>0 ?
                        aopt0(it+t,ieta,ieps,ia)=0;
                        copt0(it+t,ieta,ieps,ia)=0;	
                        val0(it+t,ieta,ieps,ia)=neg;		
                    end
			
                    if ia>1
                        if val0(it+t,ieta,ieps,ia)<val0(it+t,ieta,ieps,ia-1)
                            disp('no monotonocity in value function');
                            pause;
                        end
                    end
							
				end
            end
        end
    end	% it = tr,tr-1,..,1

    % ____________________________________________________________________________
    %
    %	Part 2: the young's problem
    %

    % compuation of the decsion rules for the worker 
    for it=t:-1:1       % all periods t=1,2,..T 
		m0=0;

        it
        itrans
        for ia = 1:1:na     % asset holding in period it 
			asset0=agrid(ia);
            for ieps=1:1:neps

				% select the part in the value-function array v(.) where
				% the value of the one-year older retired/worker (age=it+t+1) is stored
				% for the case that the worker has type ieta and iempl next period
				% likewise: optimal consumption and labor
                
                v11=val0(it+1,1,ieps,:);		% CHECK
				v11e(1:na)=v11;			% check if this is a vector
                v11e=v11e';              % column vector

	
                v12=val0(it+1,2,ieps,:);		% CHECK
				v12e(1:na)=v12;			% check if this is a vector
                v12e=v12e';              % column vector		
			
				df=v11e(2:na)-v11e(1:na-1);
				if sum(df<=0)>0
					disp('value function v11e not monotone');
                    pause;
                end
					
				df=v12e(2:na)-v12e(1:na-1);
				if sum(df<=0)>0
					disp('value function v12e not monotone');
                    pause;
                end

				% first step: l=0 optimal solution?
				% 1. computation of optimal a' for l=0
				% 2. computation of optimal a' for l=lmin, small constant
				
				for ieta = 1:1:2
					w0=wseq(it)*efage(it)*eps1(ieps)*eta1(ieta);
                   
					% triple ax,bx,cx for golden section search
					ax=0; 
                    bx=-1; 
                    cx=-2;
					vw0=neg;
					vw1=wvalue(agrid(1));
					m=0;
					% do until ax<=bx and bx<=cx;
                    while (ax>bx) || (bx>cx)
                        m=m+1;
						vw1=wvalue(agrid(m));
						if vw1>vw0
							if m==1 
                                ax=agrid(m); 
                                bx=agrid(m);
                            else	
								bx=agrid(m); 
                                ax=agrid(m-1);
                            end
							vw0=vw1;
                        else
							cx=agrid(m);
                        end
						
						if m==na 
							ax=agrid(m-1); 
                            bx=agrid(m); 
                            cx=agrid(m); 
                        end
                    end

					if ax==bx  
						a1=agrid(1);
						vw0=wvalue(a1);
						vw1=wvalue(a1+psi1);
						if vw1>vw0		% corner solution: a'=0?	
							bx=assetmin+psi1;
                            a1=golden(fhandle_function2,ax,bx,cx,tol_golden);
						end
						
					elseif bx==cx
						a1=agrid(na);
						vw0=wvalue(a1);
						vw1=wvalue(a1-psi1);
						if (vw0<=vw1) 		% corner solution: a'=assetmax?
							if vw1>neg		% both solutions for bx and cx imply c<0
								bx=a1-psi1;
							else
								bx=ax+psi1;
							end
                            a1=golden(fhandle_function2,ax,bx,cx,tol_golden);
						end
									
					else
                        a1=golden(fhandle_function2,ax,bx,cx,tol_golden);
					end
					
					aopt0(it,ieta,ieps,ia)=a1;
					w0=wseq(it)*efage(it)*eps1(ieps)*eta1(ieta);
	
					% computation of labor supply 
					if case_endogenouslabor==1;
						c=gam*( (1+(1-taur)*rseq(it))*asset0+trseq(it)+(1-tauwseq(it)-taubseq(it))*w0-a1*(1+garate) );
						if c>0 
							l0=1-c/ ( (1-tauwseq(it)-taubseq(it))*w0 ) * (1-gam)/gam;
							if l0<0
								l0=0;
								c= (1+(1-taur)*rseq(it))*asset0+trseq(it)-a1*(1+garate);
							end
                        end

					else
						l0=laborexogenous;
						c=(1+(1-taur)*rseq(it))*asset0+trseq(it)+(1-tauwseq(it)-taubseq(it))*w0*l0-a1*(1+garate);
					end

                    if c>0
                        vw01=wvalue(a1);			
						lopt0(it,ieta,ieps,ia)=l0;
						copt0(it,ieta,ieps,ia)=c;
						val0(it,ieta,ieps,ia)=vw01;
                    else
                        vw01=neg;			
						lopt0(it,ieta,ieps,ia)=0;
						copt0(it,ieta,ieps,ia)=0;
						val0(it,ieta,ieps,ia)=vw01;
                        
                    end
                
    			end     % ieta
        	end         % ieps */	
    	end             % ia
	end                 % it 

    lopt = lopt0;
    copt = copt0;
    aopt = aopt0;
    v = val0;
%	save aopt0,lopt0,copt0;	
%	if _ssfinal==1; save val0final=val0; endif;		// just a check
	
end

