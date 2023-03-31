This directory contains custom files that were used to install a 1024 antenna version of MIRIAD on a Ubuntu 22.04.2 LTS system. 

WARNING: this space just collects together the settings used on a specfic machine at a particular point in time. Please refer to the developers' [pages](https://www.atnf.csiro.au/computing/software/miriad/INSTALL.html) for the correct and full instructions.

### Bash script

The bash script, `build_miriad.sh` was used for testing different installation options within the custom files. If the script finds a directory called `miriad` in its own directory it will first delete it before re-untarring the source files into a new `miriad`. It then runs `configure`, before copying any necessary modified files from the designated directory to the necessary location in the new `miriad` directory. 
 
A binary version of miriad was installed prior to running the source code installation steps. Some libraries from the binary version were used during the source code installation.

### maxdim*.h files

These `inc/linux64/` files have been edited in order to compile miriad for SKA simulations. The number of telescopes required is specfied by `MAXANT=1024`. 

### GNUmakedefs file 

This is a modified version of the `linux64/` file produced by `configure`. It was modified after running configure. The modifications solve memory allocation and other issues. Main differences between this version of the GNUmakedefs file and the default settings:

1.

-fallow-argument-mismatch flag passed to gfortran in order to solve the following error (see https://github.com/zorkzou/Molden2AIM/issues/9) :

 Type mismatch in argument ?data? at (1); passed REAL(4) to COMPLEX(4)

2.

 -mcmodel=large  flag passed to gfortran and gcc. It prevents some errors that result from memory allocation, e.g, "...relocation truncated to fit: R_X86_64_PC32..."

3.

-fPIC passed to gcc. It prevents the C compiler errors associated with memory allocation.

4.

Full paths to rpfits and wcs libraries from the binary installation are added to GNUmakedefs.