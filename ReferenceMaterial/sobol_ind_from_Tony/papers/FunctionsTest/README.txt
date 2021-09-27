
The folder FunctionsTest contains the following files:

* FunctionsExamples.tex: 
tex file containing analytical expressions of four functions test and the associated true values of first-order and total effect Sobol' indices

* bratley.m:
matlab file defining the bratley et al. function

* gfunction.m:
matlab file defining the gfunction of Sobol'

* ishigami.m:
matlab file defining the Ishigami function

* morokoff.m:
matlab file defining the Morokoff and Caflisch function

* sobol_theoretical.mw:
Maple code defining the function soboltheoric


IMPORTANT:

The file sobol_theoretical.mw is not a text file. To be read, it must be pulled from command line. 
The file was written with Maple 17 (I do not know if the code works with older versions).
The file contains the following elements:
* the function soboltheoric that calculates either first-oder indices, closed second-oder indices or total effect indices.
* the first-order and total effect indices calculated for each of the four functions test.

I will add soon comments to the code to ease the use of the function soboltheoric.

