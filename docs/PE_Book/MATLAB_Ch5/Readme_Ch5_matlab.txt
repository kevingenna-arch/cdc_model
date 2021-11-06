Burkhard Heer, Public Economics. The Macroeconomic Perspective

The files in Ch5_Matlab_programs.rar are:
========================================

Ch5_data.m:	generates Fig. 5.7 and computes empirical business-cycle statistics
		uses input files Ch5_data_matlab.xlsx and Fig_5_1_data_bh.xlsx as inputs

Ch5_welfare_taul.m:	computes the partial and general equilibrium effects of a 1 
			percentage point change in the labor income tax (that is spent on 
			government consumption)

Ch5_welfare_tauk.m: 	computes the general equilibrium effects of a change in capital
			income taxes that is financed by a change in labor income taxe

Ch5_laffer.m:	computes the Laffer curves of capital and labor income taxation

Ch5_lucas.m:	computes the Lucas growth model


Ch5rbcstochtax.mod: computes the impulse responses to a shock on the labor and capital
		income taxes in the RBC model with stochastic taxes

Ch5rbcstochtaxlin.mod:	linear version of Ch5rbcstochtax.mod, 
			impulse responses of the variables cannot be interpreted
			as percentage deviations


Readme_Ch5_matlab.txt  (this file)


Hints:
=====
In order to run Ch5rbcstochtax.mod, you need to install 
DYNARE (freely downloadable):

http://www.dynare.org/documentation-and-support/quick-start

In order to start the program, change the path in Matlab to the directory containing the file
and type:

dynare Ch5rbcstochtax.mod



All errors in the programs are mine.