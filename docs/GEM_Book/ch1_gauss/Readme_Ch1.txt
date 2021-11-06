Gauss Computer Programs Chapter 1,
2nd edition Heer/Maussner, Dynamic General Equilibrium Modeling
===============================================================

The files

Benchmark_LL.g
Data.txt        
GetPar.g        
nettoav.txt       
SolveLA.src     
Var1.g          

store various programs referred to in Chapter 1 of the book.

Benchmark_LL.g and SolveLA.src compute impulse responses and second moments for
the benchmark model (see Figure 1.6 and Table 1.2). You must either copy the code
in SolveLA.src into Benchmark_LL.g, put an #include SolveLA.src statement at
the top of the file Benchmark_LL.g or add the files in SolveLA.src to
your user.lcg library to be able to run Benchmark_LL.g.

Var1.g estimates the impulse responses displayed in Figure 1.7. The necessary data
are stored in Data.txt.

GetPar.g computes the parameter values used in our simulation of the benchmark model
from German data. It also computes second moments from German data. Data are stored
in Data.txt and nettoav.txt.

The programs run under Gauss 7.0 and higher.

In case you have any questions please contact me at 

alfred.maussner@wiwi.uni-augsburg.de

