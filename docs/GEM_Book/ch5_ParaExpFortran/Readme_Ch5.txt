Fortran Computer Programs Chapter 5,
2nd edition Heer/Maussner, Dynamic General Equilibrium Modeling
===============================================================

The files in

LP.zip

    GaussNewton.for
    LP.dsp         
    LP.dsw         
    LP.for         
    LP.opt         
    MNR.for        
    Parameters.txt 
    Psi.for        
    QN.FOR         
    QuickWin.for   
    SolveLA.for    
    ToolBox.for                 

    Lp.exe

store the Fortran code that implements the parameterized expectations solution of the limited participation
model of Chapter 5 of our book.

The files in

SGNNI_b.zip

    Differentiation.for              
    Function.for                     
    MinPack.for                      
    MNR.for                          
    Parameters.txt                   
    PE.FOR                           
    QuickWin.for                     
    Search1.for                      
    SGNNI.dsp                        
    SGNNI.dsw                        
    SGNNI.opt
    SGNNI.exe                        
    SGNNI_b.FOR                      
    Simulations.for                  
    ToolBox.for                      

store the Fortran code that implements the parameterized expectations solution of the stochastic
growth model with a non-negativity constraint on investment as described in
Chapter 5 of our book.


The code is Fortran 95 and some procedures require not only Fortran 95 intrinsic functions
but call routines from the IMSLF and Compaq CXML libraries. If you do not have access to these
libraries you will not be able to compile and link the respective peace of program code and you must
find a work around. Also, some lines of code - input and output to screen windows - may work
only in the Windows Visual Studio environment which we used to compile and link our programs. You may
want to out-comment these lines of code or replace them with code adequate for your environment. 

You should, however, be able to run the executables Lp.exe and SGNNI_b.exe, respectively. The program's options
are supplied by the file Parameters.txt, which you can edit in order to change the
parameters of the model and in order to choose different solution methods.

Though I will not be able to give you advise on how to use our programs in any specific Fortran
environment I will answer your questions regarding those pieces of Fortran code that implements
certain numerical methods. You can contact me via my email address: 

alfred.maussner@wiwi.uni-augsburg.de

