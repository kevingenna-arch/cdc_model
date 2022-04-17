@#for i in 1:NT
	% productivity by age
	a_NQ_@{i}=.7+log(@{i})/10;
	a_Q_@{i}=1+log(@{i})/10;
@#endfor
