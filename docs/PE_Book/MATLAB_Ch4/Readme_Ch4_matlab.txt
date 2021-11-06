Burkhard Heer, Public Economics. The Macroeconomic Perspective

The files in Ch4_Matlab_programs.rar are:
========================================


Ch4_data.m:	generates Fig. 4.1 and computes empirical business-cycle statistics
		uses input files FRED_data.txt and Fig_4_1_data.xlsx

Ch4_subs_private.m	computes Fig. 4.9 and 4.10, uses many different input files

Ch4_subs_private_pub_dyn.m	computes the transition in Fig. 4.11, using many different 
			input files (findkinitial.m, findclpk.m, findclp.m,...)
			PERMANENT INCREASE of G by 1%


Ch4_subs_private_pub_dyn1.m	computes the transition in Fig. 4.11, using many different 
			input files (findkinitial.m, findclpk.m, findclp.m,...)
			TEMPORARY INCREASE of G by 1%


Ch4rbc1.mod: 	computes the stochastic Ramsey model with goverment consumption with
		Cobb-Douglas utility function

Ch4rbc1linear.mod:	linear version of Ch4rbc1.mod, 
			impulse responses of the variables cannot be interpreted
			as percentage deviations

Ch4rbc2.mod:	computes the stochastic Ramsey model with goverment consumption with
		additively separable utility function (in consumption and labor)

Ch4rbc2linear.mod:	linear version of Ch4rbc2.mod, 
			impulse responses of the variables cannot be interpreted
			as percentage deviations


Ch4newkeynes.mod:	computes the impulse responses of the New Keynesian model in Chapter 4

Ch4newkeyneslin.mod:	linear version of Ch4newkeynes.mod, 
			impulse responses of the variables cannot be interpreted
			as percentage deviations


Readme_Ch4_matlab.txt  (this file)


Hints:
=====
In order to run Ch4rbc1.mod, Ch4rbc1.mod, and Ch4newkeynes.mod you need to install 
DYNARE (freely downloadable):

http://www.dynare.org/documentation-and-support/quick-start

In order to start the program, change the path in Matlab to the directory containing the file
and type:

dynare Ch4rbc1.mod




All errors in the programs are mine.