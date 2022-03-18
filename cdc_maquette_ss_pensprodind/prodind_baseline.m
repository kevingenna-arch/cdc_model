@#for i in 1:NT
	a_NQ_@{i} = .8 + normpdf(@{i}, 6, 4/1.2);
	a_Q_@{i} = 1.1 + normpdf(@{i}, 6, 3/1.2);
@#endfor