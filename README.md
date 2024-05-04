  mdotwd: Simulation of a double White Dwarf binary system with FLASH
  ======

  This is a FLASH4 setup, together with a Gravity Unit (see the Roche directory). It simulates the accretion from one White Dwarf (WD) to another through Roche lobe overflow.

  Some notes:
  
  1. The problem does not use self-gravity, it depends on a  gravity unit called Roche
  2. To use the gravity unit, put it in source/Gravity/GravityMain/Roche
  3. After initial calculations, the code dumps some useful parameters in the files
   
  - rho.par - Stores some results from the Lane-Emden solver.
  - initial_don_profile.dat - Saves the donor profile.
  - initial_acc_profile.dat - Saves the accretor profile.

## Requirements:

- [FLASH](https://flash.rochester.edu/site/flashcode.html) installed (publicly available, high-performance computing, multi-physics application code)
- Access to High-Performance Parallel Computing HPC; computations are numerically demanding
- [yt project](https://yt-project.org/) for plots and various figures

## Setup, compiling and running my FLASH simulation:

Assuming FLASH has been installed, follow the steps to configure, compile, and run the code. The first step below specifies that we want to use the /Roche files to model gravity, in 2 dimensions, using cartesian coordinates:
  
  -  ./setup mdotwd -2d +cartesian -auto -with-unit=physics/Gravity/GravityMain/Roche -maxblocks=3000
  - make -j 4
  -  ./FLASH
