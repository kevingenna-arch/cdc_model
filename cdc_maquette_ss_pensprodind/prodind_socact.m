@#for i in 1:NT
	a_NQ_@{i} = .8 + normpdf(@{i}, 6+2, 4/1);
	a_Q_@{i} = 1.1 + normpdf(@{i}, 6+2, 3/1);
@#endfor