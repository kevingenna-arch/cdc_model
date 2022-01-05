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

##### xpop
Variable mystérieuse qui donne des informations sur la population au tout debout des générations. 
Notamment, si $\pi_s$ est la proportion entre qualifiés et non qualifiés on a $P_{s,1} = \pi_s \times xpop$, puis chaque cohorte évolue selon sa propre dynamique $P_{s,j,t} = \beta_{s,j,t-1} P_{s,j-1,t-1} + mig_{s,j,t}$.


### Calibration population

Le choc à cibler est celui de migration. 
Pour ce faire on porte $mig_{s,j}$ comme variable endogène, et on transforme $P_{s,j}$ en variable exogène. 
Or, les chocs migratoires sont $(T-1)\times S$, où $T$ est le nombre total de générations, 17, et $S$ les deux niveaux de qualification. 
Cela laisse une différence de deux variables, notamment $P_{s,1}$.



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