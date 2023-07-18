! Usage:
! python3 onclick.py 


! the required files (are all optional) 
! for harmonic parameterisation: filename.fch or filename.fchk (this is the vibrational analysis made by Gaussian)
! for charge fitting: filename.esp (this is the electrostatic potential calculations made be Gaussian)
! for Lennard Jones paraemeters: filename.mol2 or filename.pdb  
! for torasional parameters: filename.scn (these are the potential energy scan made by Gaussian, multiple files allowed)  
! the configuration file: cfg 


! for guidance on how to introduce your preferences in the cfg file read the comments below 
! in this file you can enter your preferences here line by line 






! ===================================================================================================
! WARNING ::: your comment lines have to be starting with "!" 
!             otherwise confusion with the keywords is possible

!             you have to make sure that all of your entries are consistent, 
!             otherwise the program will try to make them consistent his way 
!             and this can't be guaranteed each time


! for the details of the parameterization processes check the *.log file  
! the preferences in this file will be copied in the log file  
               

!____________________________________________________
!  General options ____________________________
! to disable the equivalencing scheme of the program use "nosymm"

! give your run an indicative name if you like 
runtag nma_mcsa


!____________________________________________________
!  Lennard-Jones options ____________________________

! to give an atom fixed Lenard-Jones parameters use 
! flj atom_index  epsilon Rmin epsilon_1-4 Rmin_1-4  (no zeros here)
! for example:  
! flj 3 -0.024    1.34  -0.024    1.34
! you can introduce more lines the same way 

! to force the program to assign multiple atoms the same Lenard_Jones parameters use 
! elj atom1_index  atom2_index  atom3_index  ..... 
! for example:   
! elj 2 3 4 



!____________________________________________________
!  charge options ___________________________________

! to give an atom a fixed charge use 
! fc atom_index  charge_value 
! for example:  
! fc 3 -0.148


! to force the program to assign multiple atoms the same charge use 
! ec atom1_index  atom2_index  atom3_index  ..... 
! for example:   
! ec 2 3 4 


! to force the program to assign a total charge to a group of atoms use 
! gc atom1_index  atom2_index  atom3_index ..... total_charge 
! for example:   
! gc 9 10 11 0.0 


! the charges to be kept constant after the first stage 
! in RESP these charges are the non aliphatic Hydrogens 
! of course you can extend this possibility as you wish 
! to do so use 
! kc atom1_index  atom2_index ... 
! for example: 
! kc 5 8 13 



!____________________________________________________
! torsional options _________________________________

! if you wish to only fit the torsional parameters use the keyword "qm_mm"  
! then the program will read the QM and MM conformational energies
! and the torsional angles with the desired periodicities, a full hypothetical examples is stated below   
qm_mm
qm  1 2 3 3   
mm  3 4 5 9 
dihedral 1 2 3 4 7 8 9 7  
dihedral 2 3 4 5 8 9 6 3
dihedral 3 4 5 6 0 9 8 5   



! to give a torsional angle fixed periodicities use 
! fp atom1_index  atom2_index  atom3_index  atom4_index  n1 n2 ... 
! note, the first four numbers are the atom indexes and the rest of numbers are the periodicities
! for example:  
fp 2 1 5 6 3
fp 5 7 9 10 3


! to give a torsional angle fixed parameters use 
! ft atom1_index  atom2_index  atom3_index  atom4_index  force_constant periodicity phase_angle 
! for example:  
ft 2 1 5 6 1.324 1 0
ft 2 1 5 6 0.543 2 180
ft 2 1 5 8 3.897 1 0 


! it is optional to change the fitting constraints, here are the default values  
phases 
kmax  7
maxenergy  30
!nsteps  20000
!weight b 300

! and here is the definition
! to force the phase angles to take values of either 0 or 180 use "phases 0"
! to allow the phase angle to take any value and the force constants to be positive or negative use "phases"
! otherwisw the program will assign either 0 or 180 to the phase shifts and the force constants values > 0
! kmax is the maximum value the force constants can reach 
! maxenergy is the maximum energy above which the conformers are ignored 
! tnsteps is the number of torsional fitting iterations 
! to give more or less weight for conformers use weight 
! to give more or less weight for conformers whose energy below a cutoff use 
! weight c cutoff_value weight_value
! alternatively you can give weights according to Boltzmann distribution, then use 
! weight b temperature_value




 










