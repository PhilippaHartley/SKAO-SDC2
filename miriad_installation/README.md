This directory contains files that were used to install a 1024 antenna version of MIRIAD on a Ubuntu 22.04.2 LTS system. 

The bash script is intended for testing different installation options within the above custom files. If the script finds a directory called miriad in its own directory it will first delete it before re-untarring the source files into a new miriad. It then runs configure, before copying any necessary modified files from the designated directory to the necessary location in the new miriad directory.

Main differences in GNUmakedefs file from default settings produced by configure:

1.

-fallow-argument-mismatch flag passed to gfortran in order to solve the following error (see https://github.com/zorkzou/Molden2AIM/issues/9) :

 Type mismatch in argument ?data? at (1); passed REAL(4) to COMPLEX(4)

2.

 -mcmodel=large  flag passed to gfortran and gcc. It prevents some errors that result from memory allocation, e.g, "...relocation truncated to fit: R_X86_64_PC32..."

3.

-fPIC passed to gcc. It prevents the C compiler errors associated with memory allocation.

4.

full paths to rpfits and wcs libraries from the binary installation are added to GNUmakedefs