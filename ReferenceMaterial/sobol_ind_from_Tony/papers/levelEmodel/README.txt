# Steps to compile the fortran source
# https://gcc.gnu.org/wiki/GFortranGettingStarted

# First step:
# On all systems:

gfortran -c batman.f
gfortran -c bios.f
gfortran -c coeff.f
gfortran -c gtm1.f
gfortran -c levele_model.f
gfortran -c mgspac.f
gfortran -c mgtime.f
gfortran -c nearf.f
gfortran -c setzer.f
gfortran -c sift.f
gfortran -c tres.f


# Second step:

# On Unix systems:

gfortran levele.f  batman.o bios.o coeff.o gtm1.o levele_model.o mgspac.o mgtime.o nearf.o setzer.o sift.o tres.o -o levele.out


# On Windows:

gfortran levele.f  batman.o bios.o coeff.o gtm1.o levele_model.o mgspac.o mgtime.o nearf.o setzer.o sift.o tres.o -o levele.exe


# Once the files are compiled the model is called by:
./levele.out (or ./level.exe on Windows)


# IMPORTANT: The model require a .txt file named « try.txt » containing:
# on the first line the total number of points on which the model is to be evaluate
# then, a matrix containing the evaluation points
# the file must be place in the « /src » repository that contains the fortran files
 

