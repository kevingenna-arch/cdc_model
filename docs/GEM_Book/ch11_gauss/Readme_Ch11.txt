Gauss Computer Programs Chapter 11,
2nd edition Heer/Maussner, Dynamic General Equilibrium Modeling
===============================================================

The files

Derivatives.src
Filter.src
Function.src
Integration.src
Nleq.src
Optimization.src
Search.src
Tools.src
Tools.dll

store the Gauss procedures that implement the various numerical methods
considered in Chapter 11 of the book. It should be obvious from the files' names - if not noted
in the book - which procedures are stored in which file.

In order to use these programs you either can copy them into your main file so that this file
contains all code required to solve your problem or build a Gauss library so that these
procedures are available to all of your programs. Please see the Gauss reference on how
to build Gauss libraries.

To use the procedure SolveLA2 you must copy the file Tools.dll into
the subdirectory \dlib of the directory where the Gauss executable file
resides (e.g. \program files\gauss\dlib)


The programs run under Gauss 7.0 and higher.

In case you have any questions please contact me at 

alfred.maussner@wiwi.uni-augsburg.de
