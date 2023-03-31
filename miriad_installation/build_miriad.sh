#!/bin/bash
export LD_LIBRARY_PATH=/usr/local/miriad-binary/linux64/lib
sudo rm -r miriad
sudo bzcat miriad-code.tar.bz2 | sudo tar xvf -
sudo bzcat miriad-common.tar.bz2 | sudo tar xvf -
cd miriad
sudo ./configure
sudo cp /home/p.hartley/software/miriad/maxdim.h inc/linux64/
sudo cp /home/p.hartley/software/miriad/maxdimc.h inc/linux64/
sudo cp /home/p.hartley/software/miriad/GNUmakedefs linux64/
sudo make