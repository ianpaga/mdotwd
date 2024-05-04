  mdotwd
  ======

  This is a FLASH4 setup, along with a Gravity Unit (see the Roche directory),
  it is meant to simulate accretion from one WD to another through Roche lobe overflow. 

  Some notes:
  1. Problem does not use self gravity, it depends on a  gravity unit called Roche. 
  2. To use the gravity unit, put it in source/Gravity/GravityMain/Roche. 
  3. After initial calculations, the code dumps some useful parameters in the files
  - rho.par - Stores some results from the Lane-Emden solver.
  - initial_don_profile.dat - Saves the donor profile.
  - initial_acc_profile.dat - Saves the accretor profile.