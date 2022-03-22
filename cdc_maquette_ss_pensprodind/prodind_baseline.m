@#for i in 1:NT
	a_NQ_@{i} = .85 + normpdf(@{i}, 6, 4/1.5);
	a_Q_@{i} = 1.15 + normpdf(@{i}, 6, 3/1.5);
@#endfor