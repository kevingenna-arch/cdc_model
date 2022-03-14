@#for i in 1:NT
	a_NQ_@{i} = .7 + normpdf(@{i}, 6, 4);
	a_Q_@{i} = 1 + normpdf(@{i}, 6, 3);
@#endfor