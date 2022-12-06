# PSO_code
Particle Swarm Optimisation code to fit CEST spectra

This code can be compiled using the provided Makefile (just run make in the current directory).
Arguments needed for the code to run are the following:
First argument is the name of the file with the spectra inside: the filename should start with spect_XYZ, but only the substring XYZ is needed). 
Second argument is the B1 measured at the position of the spectra given.
Third argument is a substring given to the results output file.

The ouput file is a csv file containing all the concentration, exchange rate, T2 and fitted chemical shift of the various pools, together with the observed T1 obs and T2 of the water pool, as well as the error in the fit over all the frequencies measured. 

A lot of parameters can be set inside the PSO_v5.c file at the moment: this may be included inside a separate file subsequently.
