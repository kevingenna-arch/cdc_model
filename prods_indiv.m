x = 0:9;

microprod = figure(1);
subplot(2, 2, 1);
plot(x*5 + 20, 1 + log(x)/10, ...
     x*5 + 20, .7 + log(x)/10);
legend('Q', 'NQ', 'Location','southeast')
title('Actuelle')
ylim([.6 1.3])

subplot(2, 2, 2);
plot(x*5 + 20, 1 + normpdf(x, 6, 3/2), ...
     x*5 + 20, .7 + normpdf(x, 6, 4/2));
legend('Q', 'NQ', 'Location','southeast')
title('Baseline')
ylim([.6 1.3])

subplot(2, 2, 3);
plot(x*5 + 20, 1 + normpdf(x, 6, 3), ...
     x*5 + 20, .7 + normpdf(x, 6, 4));
legend('Q', 'NQ', 'Location','southeast')
title('Formation continue')
ylim([.6 1.3])

subplot(2, 2, 4);
plot(x*5 + 20, 1 + normpdf(x, 6+2, 3/1), ...
     x*5 + 20, .7 + normpdf(x, 6+2, 4/1));
legend('Q', 'NQ', 'Location','southeast')
title('Activités Socialisés')
ylim([.6 1.3])
