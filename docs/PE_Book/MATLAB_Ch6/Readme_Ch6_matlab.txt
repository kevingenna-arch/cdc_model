Burkhard Heer, Public Economics. The Macroeconomic Perspective

The files in Ch6_matlab_programs.rar are:
========================================


Ch6_social_security1.m: computes the steady-state and transition effects of 
	an abolition of pay-as-you-go pensions, case 1: inelastic labor supply


Ch6_social_security2.m: computes the steady-state and transition effects of 
	an abolition of pay-as-you-go pensions, case 2: elastic labor supply


Ch6_social_security3.m: computes the steady-state effects 
	an abolition of pay-as-you-go pensions, case 3: contributions-based pensions

Ch6_social_security4.m: computes the steady-state effects
	an abolition of pay-as-you-go pensions, case 4: growth

Ch6_optimal_pension_ss.m: 	computes the stationary equilibrium for
	the large-scale OLG model with individual uncertainty in the Section 9.3.1
	NOTICE: The computation takes VERY LONG. On my computer, it took approx 6 weeks!! 
	The code is a one-to-one translation of the corresponding Gauss code, which
	only takes 11 minutes!!
	Before you start the program, take caution that the program can run uninterrupted
	for such a long period.


Input files:
============

survival_probs_US.xlsx:	Excel input file for Ch6_optimal_pension.g with US population
	data

efficiency_profile.xlsx: efficiency-age profile estimated by Hansen (input into Ch6_optimal_pension_ss.m)


Readme_Ch6_matlab.txt  (this file)


Instructions:
=============

In order to run the matlab programs, you need to store all files in the same directory and make this 
the working directory.

CAUTION: The computation time of ch6_optimal_pension_ss.m is extremely long (many weeks!).
All errors in the programs are mine.
