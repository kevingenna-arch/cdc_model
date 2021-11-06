Burkhard Heer

The files in Ch1_Matlab_programs.rar are:
========================================


Ch2_data.m	(main program file, described in Chapter 2)
Ch2_Ramsey1.m  	(main program file, described in Chapter 2)
Ch2_Ramsey2.m   (main program file, described in Chapter 2)
Ch2_rbc.m  	(main program file, described in Chapter 2)

solveLA.m 	(source file, which computes the solution to the log-linearized RBC model; the routine is taken
		from burnside. All errors are mine.)
hpfilter.m	(source file, which computes the hp-filtered series)
solab.m 	(source file, which computes the solution to linear system of difference equations; the routine is taken
		from Paul Klein. All errors are mine.)
jac_cd		(source file, which computes the Jacobian using central differences. The routine is based on Doornik's code. All errors are mine. )


FRED_data_2_24okt2015.txt	(data 1947.1-2015.2, excluding PCE index for prices)
FRED_data_2_28okt2015.txt;	(data 1959.1-2015.2, including PCE index for prices)

Readme_Ch2_matlab.txt  (this file)

Instructions:
=============

In order to run the Matlab programs, you need to store all files in the same directory and make this 
the working directory.

In order to run Ch2_rbc.m, I use the routine solveLA that is provided by Burnside. You can retrieve this file
and other Matlab code from Christian Zimmermann's Repec QM&RBC codes online page:

https://dge.repec.org/codes.html

There, you also find a link to the Matlab code for RBC models and GMM estimation by Craig Burnside.

Help by Sijmen Duineveld is gratefully acknowledged.

All errors in the programs are mine.


