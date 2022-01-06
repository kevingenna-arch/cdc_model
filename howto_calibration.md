# Calibration


Cette note reporte les étapes suivie pour effectuer la calibration du modèle CDC sur la base des nouvelles projections de l'INSEE à 2100. 
Cela sert de guide de référence pour toute calibration future.

#### Idée générale
L'approche est de remplacer une variable endogène qu'on est intéressé à calibrer par une variable exogène qui soit assez flexible pour reproduire la dynamique d’intérêt de la variable endogène.
Prenons un modèle simple comme le suivant
$$Y = BX + \epsilon$$
où $Y$ est la variable endogène d’intérêt, $X$ une variable exogène, $B$ un paramètre autrement calibré, et $\epsilon$ un terme d'erreur.

La stratégie est de considérer temporairement $Y$ comme une variable exogène avec un sentier/évolution donné par les données, et $X$ comme une variable endogène.
Si on impose un sentier précis pour $Y$, le solver de Dynare trouvera des valeurs de $X$ qui matchent $Y$, ainsi que des chocs compatibles.

On peut donc sauvegarder ces valeurs de $X$ (et, le cas échéant, de $\epsilon$) et les utiliser comme valeurs pour le modèle une fois qu'on considère à nouveau $X$ comme variable exogène et $Y$ comme variable endogène.

#### Autres Variables

##### `xpop`
Variable mystérieuse qui donne des informations sur la population au tout debout des générations. 
Notamment, si $\pi_s$ est la proportion entre qualifiés et non qualifiés on a $P_{s,1} = \pi_s \times xpop$, puis chaque cohorte évolue selon sa propre dynamique $P_{s,j,t} = \beta_{s,j,t-1} P_{s,j-1,t-1} + mig_{s,j,t}$.


### Calibration population

Le choc à cibler est celui de migration. 
Pour ce faire on porte $mig_{s,j}$ comme variable endogène, et on transforme $P_{s,j}$ en variable exogène. 
Or, les chocs migratoires sont $(T-1)\times S$, où $T$ est le nombre total de générations, 17, et $S$ les deux niveaux de qualification. 
Cela laisse une différence de deux variables, notamment $P_{s,1}$.
Ces deux variables dépendent de `xpop`, donc on éteint temporairement cette variable aussi.
Comme `initval`s pour les variables de population on utilise des valeurs proches de celles d’état stationnaire précédemment trouvées, $.8$.
En revanche, on utilise le bloc `shocks;` pour passer les valeurs historiques des variables de population, $P_{s, j, t}$.

Afin de sauvegarder les résultats, on utilise le commande `save_params_and_steady_state` qui permet de stocker les valeurs des paramètres et des variables à l’état stationnaire.
De plus, on simule beaucoup plus de périodes que prévu, notamment 100 au lieu de 40 quinquennats.
Les chocs seront donc une séquence au milieu de ces périodes.

Le modèle sera donc résolu à partir des `initval` données, Dynare définira un ES cohérent et simulera $\pm 500$ ans.
Si les chocs frappent en $t'$ les agents commenceront à anticiper ces chocs avant, comme le modèle est déterministe.
Une fois la simulation terminée, les séries temporelles sont sauvegardées et les résultats extraits de ces séries.
Dans ce cas on sauvegarde les séries des chocs migratoires.


_Nota bene sur les données_: les données utilisées dans ce cas sont celles historiques de l'INSEE fusionnées avec les projections du scenario central jusqu'à 2121.
La période d’intérêt couvre finalement de 1900 à 2100, avec des données officielles.
On agrège les observations par an et par age unitaire à des quinquennats (de 4.99 ans), pareil pour les classes d'age sauf une agrégation complète pour les ultra-centenaires.
Pour la normalisation des niveaux des classes de population on prend comme référence la classe des $0-4$ en 2000, l'année centrale.
Comme le modèle se base sur 17 générations, on ne considère pas les données des classes $0-20$.
Les mêmes données de population sont pour l'instant utilisées pour tout niveau de qualification.





#### Ressources
+ Dynare forum
	* [Solving for SS](https://forum.dynare.org/t/olg-model-solving-for-steady-state/6327)
	* [Simulating and calibrating](https://forum.dynare.org/t/simulating-olg-models-with-real-data/443)
	* [Optimal taxation](https://forum.dynare.org/t/olg-optimal-taxation/5538)
	* [Diamond 1965 basic OLG](https://github.com/davidrpugh/pyeconomics/wiki/OLG-Models)
	* [General DF search](https://forum.dynare.org/search?q=OLG)
+ Dynare helpfile
	* [initval, endval, histval](https://www.dynare.org/manual/the-model-file.html#initial-and-terminal-conditions)
	* [SS](https://www.dynare.org/manual/the-model-file.html#steady-state)
	* [Chocs](https://www.dynare.org/manual/the-model-file.html#shocks-on-exogenous-variables)
	* [vartype flipping](https://www.dynare.org/manual/the-model-file.html#change_type)
	* [`save_params_and_steady_state` pour la calibration](https://www.dynare.org/manual/the-model-file.html#save_params_and_steady_state)

#### Autres notes
+ pull out var names from simulations and simulations series: 
```
simuls = array2table([oo_.endo_simul',oo_.exo_simul]);
simuls.Properties.VariableNames = [M_.endo_names; M_.exo_names];
writetable(simuls, 'simuls.csv', 'Delimiter', ',')
% select only matching vars
matched = ~cellfun('isempty', regexp(simuls.Properties.VariableNames, 'mig_', 'once'));
mig_shocks = simuls(:, simuls.Properties.VariableNames(matched));
mig_flipped = rows2vars(mig_shocks);
mig_flipped.Properties.RowNames = table2array(mig_flipped(:, 1));
writetable(mig_flipped(:, (start:end)+1), 'mig_shocks.xlsx', 'WriteVariableNames',false);
% la derniere ligne laisse une ligne de plus (sans nom) dans le fichier
```
+	_NB_: Dynare rajoute toujours une periode de plus à $t-1$ pour les simulations.
Cela change le timing des chocs aussi: un choc prevu à $t\in (12, 15)$ se trouve en vrai en $t\in(12, 16)$