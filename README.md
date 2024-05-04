  mdotwd: Simulation of a double White Dwarf binary system with FLASH
  ======

  This is a FLASH4 setup, together with a Gravity Unit (see the Roche directory). It simulates the accretion from one White Dwarf (WD) to another through Roche lobe overflow.

  Some notes:
  
  1. The problem does not use self-gravity, it depends on a  gravity unit called Roche
  2. To use the gravity unit, put it in source/Gravity/GravityMain/Roche
  3. After initial calculations, the code dumps some useful parameters in the files
   
 ![DWD1-ezgif com-video-to-gif-converter](https://github.com/ianpaga/mdotwd/assets/57350668/7f6e663c-e2f9-446a-a0cb-82f1e05926e0)

## Requirements:

- [FLASH](https://flash.rochester.edu/site/flashcode.html) installed (publicly available, high-performance computing, multi-physics application code)
- Access to High-Performance Parallel Computing HPC; computations are numerically demanding
- [yt project](https://yt-project.org/) for plots and various figures

## Setup, compiling and running my FLASH simulation:

Assuming FLASH has been installed, follow the steps to configure, compile, and run the code. The first step below specifies that we want to use the /Roche files to model gravity, in 2 dimensions, using cartesian coordinates:
  
  -  ./setup mdotwd -2d +cartesian -auto -with-unit=physics/Gravity/GravityMain/Roche -maxblocks=3000
  - make -j 4
  -  ./FLASH

<img width="731" alt="roche" src="https://github.com/ianpaga/mdotwd/assets/57350668/b356a297-8bda-4e07-a894-3585f308741a">
![initialconditionDDgrid_Slice_z_density](https://github.com/ianpaga/mdotwd/assets/57350668/941d789b-ae79-4417-a892-dbfcb087b06f)
![DDzoomed_Slice_z_density](https://github.com/ianpaga/mdotwd/assets/57350668/198d935e-7e8c-4a81-826f-c2e64a0c71fe)
