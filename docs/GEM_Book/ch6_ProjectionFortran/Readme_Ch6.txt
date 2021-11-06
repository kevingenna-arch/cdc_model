Fortran Computer Programs Chapter 6,
2nd edition Heer/Maussner, Dynamic General Equilibrium Modeling
===============================================================

The files in

Equity.zip

    Differentiation.for  
    Equity.dsp
    Equity.dsw
    Equity.exe           
    Equity.for
    Equity.opt
    Function.for     
    Integration.for
    LA.for               
    MinPack.for                     
    MNR.for  
    Parameters.txt     
    Parameters_Explained.txt             
    QA.for                                                 
    QuickWin.for                    
    Search1.for                     
    Simulations.for 
    SolveLA.for
    SolveQA.for                
    ToolBox.for                     
    
store the Fortran code that implements the projection solution of the equity premium model
considered in Chapter 6 of our book.

    
The files in
    
SGNNI_c.zip
    
    Differentiation.for              
    Function.for                     
    MinPack.for                      
    MNR.for   
    Optimization.for                       
    Parameters.txt  
    Parameters_Explained.txt                 
    Projection.for
    QuickWin.for                     
    Search1.for                      
    SGNNI_c.dsp                        
    SGNNI_c.dsw                        
    SGNNI_c.exe
    SGNNI_b.for                      
    SGNNI_c.opt
    Simulations.for                  
    ToolBox.for                      
    
store the Fortran code that implements the projection solution of the stochastic
growth model with a non-negativity constraint on investment as described in
Chapter 5 of our book.
    
    
The code is Fortran 95 and some procedures require not only Fortran 95 intrinsic functions
but call routines from the IMSLF and Compaq CXML libraries. If you do not have access to these
libraries you will not be able to compile and link the respective peace of program code and you must
find a work around. Also, some lines of code - input and output to screen windows - may work
only in the Windows Visual Studio environment which we used to compile and link our programs. You may
want to out-comment these lines of code or replace them with code adequate for your environment. 
    
You should, however, be able to run the executables Equitz.exe and SGNNI_c.exe, respectively. The program's options
are supplied by the file Parameters.txt, which you can edit in order to change the
parameters of the model and in order to choose different solution methods.
    
Though I will not be able to give you advise on how to use our programs in any specific Fortran
environment I will answer your questions regarding those pieces of Fortran code that implements
certain numerical methods. You can contact me via my email address: 
    
alfred.maussner@wiwi.uni-augsburg.de
    
    
