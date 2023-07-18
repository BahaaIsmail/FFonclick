
# this version is documented 

from os import walk , getcwd , system , remove
import sys
from functions import *

# creating tempars file: to store the parameters calculated for harmonic, LJ and charges parameters 
tf = open('tempars','w')  
tf.close()

#_____________________________________________________
##  exploring the input files  
hessf , espf , mol2 , scnfs , cfg , srts = 0 , 0 , 0 , [] , 0 , []
path = getcwd()
for root , dirs , files in walk(path) :
    for f in files :
        if '.fch' in f or '.fchk' in f :     # formatted checkpoint file necessary for harmonic parameterization, and a candidate for topology.  
            hessf = f            
        elif '.esp' in f :                   # electrostatic potential file necessary for charge fitting.
            espf = f 
        elif '.scn' in f :                   # potential energy scan file(s)  necessary for  torsional fitting if not qm_mm is used. 
            scnfs += [f]  
        elif '.mol2' in f :                  # necessary for atom typing of LJ parameters if not hessian file is  there.  
            mol2 = f
        elif 'cfg' == f :                    # confoguration file where the parameterization options are introduced (optional).  
            cfg = f
        elif '.srt' in f :                   # the main output file(s) which containes the parameterized force field.  
            srts += [f]	


#___________________________________________________________________________________________________________________
## the name of the output files 
''' two output files would be supplied by the program 
    1- the 'srt' file which containes the parameterized force field. 
    2- the 'run' file which contains  important details on the parameterization procedures. 
- Because the user ususally needs to make several trials on the same molecules to reach the best force field he should give a name for the 
  output file at each trial. He can assign an indicative name by using the keyword 'runtag somename' in the configuration file
- Otherwise the program will assign an indexed name after the name of the directory at each trial e.g. mymolecule1.strt, mymolecule2.srt, etc
- So in the latter the program needs to check the index of the last saved srt file and makes a +1 increment 
- If the user assigns his own names he must be carefull in the context of changing the runtag everytime he makes a new trial, otherwise 
  the program will override the last saved 'srt' and 'run' files. 
- One last thing; no spaced name is recommended for the outputfiles. if spaces exist, the program 
  will replace with underscores and give a notification on that
'''

name = 0
if cfg: 
    sf = open(cfg)
    for line in sf : 
        if 'runtag' in line :             
            line = line.split()
            if len(line) > 2 :
                name = ''.join([i+'_' for i in line[1:]]).strip('_')
                print ( 'Spaces were detected in runtg, so it is modified to be ' + name)
            elif len(line) == 2 : 
                name = line[1]
            break 
    sf.close()

if not name :     
    name = path.replace('\\',' ').replace('/',' ').split()[-1]
    n = sorted(srts)
    if len(srts) > 0 : 
        n = int(n[-1].replace('.srt','').split('_')[-1])+1
    else : 
        n = 1
    n = str(n)
    name = name + '_' + n 
    
srt = name + '.srt'
run = name + '.run' 
# finally we have named the srt and run outputfiles 
#__________________________________________________________________ 




# ___________________________________________________________________ 
# openning the run file to record all what happens during the parameterization steps
# main remarks on the parameterization steps will appear on the screen and be written in the run file 
rf = open(run,'w')
pro = '''\n\n
  ___________________________________________________
 |                                                   |
 |                   Onclickff program               |
 |                                                   |
 | Bahaa Mostafa Abdelazim, Biophysics Department,   | 
 | Faculty of science, Cairo University              |
 |                                                   |
 | For guidance in the valid input and configuration | 
 |     files visit the README file                   |
 |                                                   |
 | For detailed information on the parameterization  |
 |     steps and strategies visit the run file       | 
 |___________________________________________________|
\n\n'''
print ( pro )
rf.write(pro)


#__________________________________________________
## two optional entries in the cfg file are expected here
'''
nosymmetry: to disable the equivqlencing scheme allover the parameterization steps
            the equivalencing scheme detects the atoms which are supposed to have or contribute to the same values
qmm_mm: if the user wishes to only run torsional fitting skipping the first three parts of para
        see the cfg for mor details 
'''
nosymm = 0
qm_mm = 0
if cfg :
    sf = open(cfg)
    for line in sf :
        if not '#' in line : 
            if 'qm_mm' in line :
                qm_mm = 1
            elif 'nosymmetry' in line : 
                nosymm = 1 
    sf.close()

if nosymm : 
    rf.write ('The equivalencing scheme of OnclickFF was disabled\n\n')
else : 
    rf.write ('The equivalencing scheme of OnclickFF was enabled\n\n')
    


#__________________________________________________________
## the source of the molecular topology
if hessf : 
    rf.write ('\n\nMolecular topology will be extractd from {}  \n'.format(hessf))
elif mol2 : 
    rf.write ('\n\nMolecular topology will be extractd from {}  \n'.format(mol2))



#_____________________________________________________________
## taking a copy of the configuration file
# only the valid entries will be considered   
if cfg :    
    rf.write('\n\nA copy of the configuration file \n')
    sf = open('cfg')
    keywords = ['runtag','flj','elj','fc' , 'gc' , 'ec' , 'kc' ,'fp' , 'ft' , 'phases' , 'kmax' , 'maxenergy' , 'nsteps' , 'weight']
    incomplete = []    
    for line in sf :
        if not line.startswith('#') :                
            check = line.split()
            if check and check[0] in ['nosymmetry','qm_mm'] :
                if len (line) > 1 : 
                    rf.write(line)    
                else : 
                    incomplete += [line]
            elif intercept(check,keywords) and len(line) > 1: 
                if len (line) > 1 : 
                    rf.write(line)    
                else : 
                    incomplete += [line]                
    sf.close()


if cfg and incomplete : 
    rf.write('\n\nThe following entries were not assigned values: \n')
    for i in incomplete : 
        rf.write(i + '\n')



   


#______________________________________________________________
# In case no input files are introduced the program halts 
if not (hessf or espf or scnfs or mol2 or qm_mm) : 
    noinput = '''\n\n
  _______________________________________
 |       No input files detected         |
 | (.fchk , .esp , .scn , .mol2 , .cfg)  | 
 |                                       |
 |          Program Halted               |
 |_______________________________________|
    \n\n'''
    print ( noinput )
    rf.write(noinput)
    rf.close()
    sys.exit()




#___________________________________________________________________________________________________
### extracting the atom symbols, connectivities and bond orders from fchk if exists, or from mol2 if exists
# this piece of data is necessary for LJ parameters
# if neither fchk or mol2 files detected, the PES file will be checked for atom symbols and connectivities
# if no PES file was there, the ESP file will be chckee for atom symbols


'''
'tmp' is the temporary file which will be created to store the topology data of the molecule
(e.g. atom elements, bonds and bond orders).
and these data can be taken from the following files in order (hess, mol2, PES, ESP)
a sample of that file can be
C H H H C O . . . 
1 2 1.0
1 3 1.0
.
.
'''

if hessf and not qm_mm:
    tf = open('tmp','w')
    #_______________
    # atom symbols (the full list should be supplied here) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    d = {1:'H' , 2:'M' ,  6:'C' ,  7:'N' ,  8:'O' ,  9:'F' , 10:'E' ,
                13:'A' , 15:'P' , 16:'S' , 17:'L' , 35:'B' , 53:'I' }
    sf = open(hessf)
    # extracting the number of atoms in the molecule as 'nat' during the extraction of atomic numbers
    for line in sf :
        if 'Atomic numbers' in line :
            break
    nat = 0
    for line in sf :
        if 'Nuclear charges' in  line :
            break
        else :
            line = line.split()    # this line containes the atomic numbers which will be converted by the dictionary 'd' into atom symbols
            nat += len(line) 
            for s in line :            
                tf.write(d[int(s)] + ' ')
    tf.write('\n')

    #____________________________________________________________________________
    # bonds and bond orders
    atoms = [i for i in range(1,nat+1)]
    for line in sf : 
        if 'MxBond' in line : 
            line = line.split()
            maxb = int(line[2])   # the maximum number of bound atoms detected for an atom
            break
    for line in sf : 
        if 'IBond' in line : 
            break 
    # bonds (connectivities)
    L = []       # to collect the full list of the connected atoms to each atom in atoms order taking maxb into account 
    for line in sf : 
        if 'RBond' in line : 
            break 
        else : 
            line = line.split()
            L += [int(i) for i in line]

    # bond orders  (may be needed for fixing the aromatic bonds)
    R = []       # to collect the full list of the bond orders in atoms order taking maxb into account 
    for line in sf :
        if 'Virial Ratio' in line :
            break
        else :
            line = line.split()
            R += [float(i) for i in line]

    for i in atoms :
        j = maxb*(i-1)
        s = L[j:j+maxb]
        t = R[j:j+maxb]
        for k in range(maxb) : 
            if s[k] > i:    # this condotion will remove the repetitions and the zeros as well 
                tf.write(str(i)+ ' ' + str(s[k]) +  ' ' + str(t[k]) + '\n')                
    sf.close()
    tf.close()    
elif not qm_mm :
    nohess = '''\n\n
  ____________________________________________
 |  WARNING::: No Hessian detected            |
 |             (.fchk or .fch)                |
 |                                            | 
 |             No harmonic parameters in SRT  |
 |____________________________________________|
    \n\n'''
    print ( nohess )
    rf.write(nohess)


# if fchk doesn't exist, mol2 will be used if exists 
if mol2 and not (hessf or qm_mm) :
    tf = open('tmp','w')
    sf = open(mol2)
    for line in sf :
        if line.startswith('@<TRIPOS>ATOM') :
            break
    for line in sf :
        if '@<TRIPOS>BOND' in line : 
            break
        else :
            line = line.split()
            a = line[1].strip('0123456789')
            tf.write(a + ' ')
    tf.write('\n')
    for line in sf :
        line = line.split()
        if len(line) == 4 :
            [i,j,k] = line[1:]
            tf.write(i + ' ' + j + ' ' + k + '\n')
    sf.close()
    tf.close()



# if neither hess nor mol2 exists, we use the PES files to extract the atom symbols and the connectivities but not the bond orders
# so we cannot determine LJ parameters because bond orders are essential for that
# perhaps as a future plan we arrange a code for detecting bond orders
if scnfs and not (hessf or mol2 or qm_mm): 
    noord = '''\n\n
  ___________________________________________________________
 |  WARNING::: No Bond orders detected                       |
 |             (.fchk, .fch or .mol2)                        |
 |                                                           | 
 |             Lenard-Jones parameters generation is halted  |
 |___________________________________________________________|
    \n\n'''
    print ( noord )
    rf.write(noord)

    scnf = scnfs[0]
    tf = open('tmp','w')
    sf = open(scnf)
    bonds = []
    atoms = {}
    pool = []
    
    for line in sf : 
        if '! A' in line : 
            break 
        elif '! R' in line :             
            b = line.split()[2].strip('R(').strip(')').split(',')
            [i,j] = [b[0] , b[1]]
            bonds += [[int(i),int(j)]]
            pool += [int(i),int(j)]
    nat = int(max(pool))
    sf.close()
    sf = open(scnf)
    for line in sf : 
        if 'Symbolic' in line or 'Charge =' in line : 
            break
    for line in sf : 
        if 'Symbolic' in line or 'Charge =' in line :    # required because the order of the both wordsin espfile can be different
            break 
    for k in range(1,nat+1) : 
        a = sf.readline().split()[0]
        atoms[k] = {'a':a}        
    sf.close()
    for i in atoms : 
        tf.write(atoms[i]['a'] + ' ') 
    tf.write('\n')
    for b in bonds : 
        [i,j] = b
        tf.write(str(i) + ' ' + str(j) + ' 0 \n')
    sf.close()
    tf.close()
elif not (scnfs or qm_mm) : 
    nopes = '''\n\n
  _____________________________________________
 |  WARNING::: No PES scans detected           |
 |             (.scn)                          |
 |                                             |  
 |             Torsional fitting is halted     |
 |_____________________________________________|
    \n\n'''
    print ( nopes )
    rf.write(nopes)


# if non of hess, mol2, or PES files exists the program uses the esp file if exists to extract the atom symbols
if espf and not (hessf or mol2 or scnfs or qm_mm):
    nocon = '''\n\n
  __________________________________________________________
 |  WARNING::: No Connectivities detected                   |
 |             (.fchk, fch, .mol2, or .scn)                 |
 |                                                          |
 |             The equivalencing scheme is disabled         |
 |                 in charge fitting                        |
 |             Non-methyl Hydrogens may need to be assigned | 
 |                in the configuration file                 |
 |__________________________________________________________|
    \n\n'''
    print ( nocon )
    rf.write(nocon)
    
    tf = open('tmp','w')
    sf = open(espf)
    for line in sf : 
        if 'NAtoms=' in line : 
            nat = int(line.split()[1])
            break
    sf.close()
    sf = open(espf)
    for line in sf : 
        if 'Symbolic' in line or 'Charge =' in line : 
            break
    for line in sf : 
        if 'Symbolic' in line or 'Charge =' in line :    # required because the order of both in espfile varies
            break            
    for k in range(1,nat+1) : 
        a = sf.readline().split()[0]        
        tf.write(a + ' ' )
    tf.write('\n')
    sf.close()
    tf.close()
elif not espf : 
    noesp = '''\n\n
  ______________________________________
 |  WARNING::: No ESP file detected     |
 |             (.esp)                   |
 |                                      |
 |             Charge fitting is halted |
 |______________________________________|
    \n\n'''
    print ( noesp )
    rf.write(noesp)
	


# defining atom labels from the tmp file
if (hessf or mol2 or espf or scnfs) and not qm_mm :  
    sf = open('tmp')
    elements = sf.readline().split()
    d = {i+1:elements[i]+str(i+1) for i in range(len(elements))}    
    sf.close()










#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
#_______________________________________________________________________________________________
### Running the parameterization steps in the following order: Harmonic, LJ, charges, torsional


# to copy a temporary text file in the run file, the text file will be removed after copying
def copyin(txt) :
    rf.write('\n\n')
    sf = open(txt)
    for line in sf :
        rf.write(line)
    sf.close()
    rf.write('\n\n')
    remove(txt)


# to make sure that a part in the parameterization is successfully done an empty file is created to notify the main script 
# this function removes such files 
def done(output) : 
    for root , dirs , files in walk(path) :
        for f in files :
            if  f == output:
                remove(output)
                return 1  


#_____________________________________________________________________________________________
## Harmonic parameterization 
print ( '\nHarmonic Parameterization ..  ', end=' ')
BONDS , UREYS , ANGLES , IMPROPERS  = []  , [] , [] , []  # Lists to store lables and parameters for each harmonic coordinate (from tempars file)
if hessf and not qm_mm:
    # running the script which calculates the harmonic force constants 
    system('python3 hess.py ' + hessf + ' > harmonics')   
    copyin('harmonics')
    if done('hessdone') : 
        print ( 'Done'  )   # this means that the harmonic parameters are now stored in the tempars file 
    else : 
        print ( 'Failed,   check the run file' )   # this means that no harmonic parameters are stored in the tempars file      
    sf = open('tempars')          # the parameters will be extracted from tempars and stored in lists to be finally written to the srt file. 
    for line in sf : 
        if 'bond' in line :      
            line = line.split()
            i,j = int(line[1]) , int(line[2])      # the atoms which form the bond 
            i,j = d[i] , d[j]                      # the lables of each i and j 
            f,q = float(line[3]) , float(line[4])  # the force constant and the equilibrium value of the bond  
            BONDS += [[i,j,f,q]]
        elif 'urey' in line : 
            line = line.split()
            i,j = int(line[1]) , int(line[2])
            i,j = d[i] , d[j]
            f,q = float(line[3]) , float(line[4])
            UREYS += [[i,j,f,q]]       
        elif 'angle' in line : 
            line = line.split()
            i,j,k = int(line[1]) , int(line[2]), int(line[3])
            i,j,k = d[i] , d[j] , d[k]
            f,q = float(line[4]) , float(line[5])  
            ANGLES += [[i,j,k,f,q]]             
        elif 'improper' in line : 
            line = line.split()
            i,j,k,l = int(line[1]) , int(line[2]), int(line[3]) , int(line[4])
            i,j,k,l = d[i] , d[j] , d[k] , d[l]
            f = float(line[5])
            IMPROPERS += [[i,j,k,l,f]]
    sf.close()
else : 
    print ( 'Not requested' )
    
    

#_____________________________________________________________________________________________
## van der Waal's Interactions 
print ( '\nLennard-Jones Parameterization ..' , end=' ' )
LJs = []
if (hessf or mol2) and not qm_mm :
    if cfg :
        system('python3 typer.py  cfg  > LJ') 
    else :
        system('python3 typer.py  0  > LJ')  
    copyin('LJ')
    if done('ljdone') : 
        print ( 'Done' )
    else : 
        print ( 'Failed,   check the run file' )
    
    sf = open('tempars')
    for line in sf : 
        if 'lj' in line : 
            line = line.split()
            i = d[int(line[1]) ]
            LJs += [[i]+[k for k in line[2:]]]
    sf.close()
else : 
    print ( 'Not requested' )
    

#_____________________________________________________________________________________________
## Charge fitting 
print ( '\nCharge Fitting ..  ', end=' ' )
qs = []
if espf and not qm_mm:     
    if cfg : 
        system('python3 resp.py ' + espf + ' cfg  > charges')
    else : 
        system('python3 resp.py ' + espf + ' 0    > charges')
    sf = open('tempars')
    for line in sf : 
        if 'charges' in line : 
            qs = line.split()[1:]
            qs = [float(q) for q in qs]
    sf.close()
    Q = float(round(sum(qs),0))       
    sf.close()            
    copyin('charges')
    if done('espdone') : 
        print ( 'Done' )
    else : 
        print ( 'Failed,  check the run file' )    
else : 
    print ( 'Not requested' )
    

#_____________________________________________________________________________________________
## Torsional fitting
DIHEDRALS = []
print ( '\nTorsional Fitting ..  ', end=' ' )
if qm_mm : 
    system('python3 mcsa.py  qm_mm > torsionals')   
    copyin('torsionals')
    if done('pesdone') : 
        print ( 'Done' )
    else : 
        print ( 'Failed,   check the run file' )
elif scnfs and hessf and espf: 
    scans = ''.join([f+' ' for f in scnfs])
    if cfg : 
        system('python3 mcsa.py ' + scans + ' cfg  > torsionals')
    else :
        system('python3 mcsa.py ' + scans + ' 0   > torsionals')
    copyin('torsionals')
    if done('pesdone') : 
        print ( 'Done' )
    else : 
        print ( 'Failed,   check the run file' )       
else : 
    print ( 'Not requested' )       


if qm_mm or (scnfs and hessf and espf) :    
    sf = open('dihedrals')
    for line in sf : 
        line = line.split()
        if qm_mm :
            t = [i for i in line[:4]]
        else :
            t = [d[int(i)] for i in line[:4]]
        p = [i for i in line[4:]]
        DIHEDRALS += [t+p]
    sf.close()
    remove('dihedrals')





    
#_____________________________________________________________________________________________
#_____________________________________________________________________________________________
# exporting the parameters (toppars , SRT file)
srt = srt.split('.')[0]+'.srt'

srt = open(srt,'w')
srt.write('\nread rtf card append\n\nRESI *****          ')
if espf and not qm_mm:
    srt.write(str(Q) + '\n')
srt.write('GROUP\n')
if espf and not qm_mm:
    for i in d : 
        srt.write('ATOM %5s%10s%10s\n'%(d[i],d[i],str(qs[i-1])))
    srt.write('\n')
for b in BONDS : 
    [i,j] = b[:2]
    srt.write('BOND %4s%8s\n'%(i,j))
for m in IMPROPERS : 
    [i,j,k,l] = m[:4]
    srt.write('IMPR %4s%8s%8s%8s\n'%(i,j,k,l))
srt.write('\nEND\n\n')
    
srt.write('\nread param card flex append\n\n')
    
srt.write('BONDS\n')
if hessf and not qm_mm: 
    for b in BONDS : 
        [i,j,f,q] = b 
        srt.write(i+'%5s%10s%8s\n'%(j,str(f),str(q)))
srt.write('\nANGLES\n')
if hessf and not qm_mm: 
    for a in ANGLES : 
        [i,j,k,f,q] = a 
        srt.write(i+'%5s%5s%10s%8s\n'%(j,k,str(f),str(q)))
srt.write('\nDIHEDRALS\n')
if qm_mm or (scnfs and hessf and espf) : 
    for d in DIHEDRALS : 
        [i,j,k,l,f,p,t] = d 
        srt.write(i+'%5s%5s%5s%10s%6s%8s\n'%(j,k,l,str(f),str(p),str(t)))    
srt.write('\nIMPROPERS\n')
if hessf and not qm_mm: 
    for a in IMPROPERS : 
        [i,j,k,l,f] = a 
        srt.write(i+'%5s%5s%5s%10s   0     0.00\n'%(j,k,l,str(f))) 
srt.write('\nNONBONDED\n')
if (hessf or mol2) and not qm_mm: 
    for w in LJs : 
        [i,s,r,s4,r4] = w 
        srt.write(i+'    0.000    %8s%8s    0.000    %8s%8s\n'%(s,r,s4,r4))         
srt.write('\nEND\nRETURN')        
        
srt.close()


#remove('tempars')
#remove('tmp')



ter = '''\n
    ____________________________________
   |  Successful Termination            |
   |____________________________________|
    \n'''
print ( ter )
rf.write(ter)

rf.close()
#x = raw_input()
