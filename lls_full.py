
'''
Next time: consider the maxenergy in QME , MME , Ts , confdi , etc
'''


# Linear Least Squares Method for torsional fitting 
### 12-1-2018
import sys
from math import *
from random import random , uniform
from functions import *

## entries just for trials until we integrate it with other parts
# the temporary files temp and tempars gave to be kept to run this script 

print '\n\n'
print '__________________________________________________________________________________'
print 'Torsional Fitting ________________________________________________________________'




#___________________________________________________________________________________________________________________________________________________________
#___________________________________________________________________________________________________________________________________________________________
#         ________________________________________________________________
             
#                               THE QUANTUM SCAN FILE 
#       ___________________________________________________________________
#         ___     ___    __                        ____  ___           ___
#        /       /      /  \     |\   |           |       |    |      | 
#        |__    |      |____|    | \  |           |__     |    |      |__
#           \   |      |    |    |  \ |           |       |    |      |               
#        ___|    \___  |    |    |   \|           |      _|_   |___   |___                                            
   
from os import walk ,getcwd


# if qm_mm is mentioned in the command line then the program will use the QM, MM and dihedral values
#     introduced in the cfg file
# if not, the program will use the PES files and the cfg file 
qm_mm , cfg = 0 , 0
#check = sys.argv[-1]
check = 'cfg'
scans = ['nma1.scn']
if check == 'qm_mm' : 
    qm_mm = 1 
    cfg = 'cfg' 
elif check == 'cfg' : 
    cfg = 'cfg' 

# if nosymmetry is mentioned in the cfg file the equvalence scheme will be disabled 
nosymm = 0
if cfg : 
    sf = open('cfg')
    for line in sf : 
        if not '#' in line : 
            if 'nosymmetry' in line : 
                nosymm = 1 
                break    
    sf.close()


# if qm_mm is mentioned in the command line we have to make sure that the entries are valid  
claimed = qm_mm + 0
# lists for the QM and MM energies, energy value per conformer
# QM energies are extracted from the PES files (scn)
# whereas MM energies will be calculated using the parameters stored in the tempar file 
QME , MME , = [] , []

# dictionary for the dihedral values;
# conformer index: the list of the values the ith dihedral takes through out the the conformers   
confdi = {i:{} for i in range(1,999999)}   # 

if qm_mm :
    #__________________________________________________________
    ### extracting the QM and MM energies with the associated
    #   torsional options if existing in the cfg file
    #----------------------------------------------------------
    sf = open('cfg')
    # dictionary for the dihedrals; will be filled with the attributes of each dihedral;
    # the indices of the atoms,   
    dihedrals = {}     # {1:{'num':[1,2,3,4], 2:[],...} , 2:{}, ...}   
    d = 0
    for line in sf :
        if not '#' in line : 
            line = line.split()
            if len(line) > 1 : 
                if 'qm' in line :
                    QME = [float(q) for q in line[1:]]
                    n = len(QME)
                elif 'mm' in line :
                    MME = [float(m) for m in line[1:]]
                    if len(MME) != n :     # invalid entries; unequal lengths of QM and MM energies, then
                        qm_mm = 0          # qm_mm enquiry will be canceled and      
                        break              # no need to go further in the cfg file  
                elif 'dihedral' in line :
                    print line 
                    if len(line)-5 != n :  # invalid entries; unequal lengths of dihedral values and QM energies, then
                        qm_mm = 0          # qm_mm enquiry will be canceled and  
                        break              # no need to go further in the cfg file 
                    d += 1     # counter for the number of dihedrals 
                    dihedrals[d] = {}
                    dihedrals[d]['num'] = [int(float(i)) for i in line[1:5]]
                    values =  [float(i) for i in line[5:]]
                    for i in range(len(values)) :
                        v = values[i] 
                        confdi[i+1][d] = v 
    sf.close()
    confdi = {i:confdi[i] for i in confdi if confdi[i]}


# if one of QME, MME or confdi is empty (means the keyword is absent in the cfg file)
# qm_mm will be ignored and the program will exit
if not (QME and MME and confdi) :
    qm_mm = 0
if qm_mm == 0 and claimed == 1 :
    print '\nYou have requested torsional parameterization but you did not provide'
    print 'valid entries in the configuration file (.cfg)'
    df = open('dihedrals','w')
    df.close()
    cfg = 0
    scans , cc , QME , MME , confdi = [] , {} , [] , [] , {}
    nconf = 0
    sys.exit()


if not qm_mm : 
    #scans = sys.argv[1:-1]        
    n = len(scans)
    print '\nThe folowing {} PES files were detected'.format(n)
    for s in scans : 
        print s + '  ' ,
    print '\n'


#_______________________________________________________________
# building the internal coordinates from the tmp file
#---------------------------------------------------------------
atoms , bonds = atoms_bonds()
uniqs = get_symmetries(atoms)
if nosymm : 
    uniqs = {i:[] for i in atoms}

if not qm_mm :
    angles = get_angles(atoms)
    ureys = get_ureys(angles)
    dihedrals = get_dihedrals(atoms,bonds)
    impropers = get_impropers(atoms)

    pairs4 = []
    for d in dihedrals :
        [i,j,k,l] = dihedrals[d]['num']
        pairs4 += [sorted([i,l])]

    pair234 = [bonds[b]['num'] for b in bonds] + [ureys[u]['num'] for u in ureys] + pairs4
    pairs = []
    for i in atoms :
        for j in atoms :
            if i < j and not [i,j] in pair234:
                pairs += [[i,j]]


if scans and not qm_mm : 
    #_______________________________________________________________
    # extracting the level of theory
    #---------------------------------------------------------------
    sf = open(scans[0])
    theory = 'HF='      # required to extarct the QM energies from the summary of the output file  
    for line in sf :
        if line.startswith(' # ') :
            if 'mp2' in line or 'MP2' in line :
                theory = 'MP2='
            elif 'b3' in line or 'B3' in line :
                theory = 'DFT='
            elif 'CCSD' in line or 'ccsd' in line :
                theory = 'CCSD='
            break
    sf.close()

    #_______________________________________________________________
    ### extracting the stable conformers (coordinates and energies)
    #---------------------------------------------------------------
    cc  = []  # a list to collect the dictionaries of the cartesian coordinates of atoms for each conformer 
    QME = []  # a list to collect the QM energis (in hartrees) from all PES files  
    for f in scans :    # looping over the PES files detected 
        sf = open(f)
        qm = ''    # string to collect the QM energies together 
        for line in sf :
            if 'Standard orientation' in line :
                sf.readline() ; sf.readline() ; sf.readline()    # skipping 3 lines
                line = sf.readline()
                x = {}   # dictionary to collect the cartesian coordinates of atoms 
                for line in sf :
                    if line.startswith(' ----------------') :
                        break
                    else :
                        lin = line.split()
                        x[int(lin[0])] = [float(lin[i]) for i in range(3,6)]
            elif 'Optimized Parameters' in line:      # this indicates that a stable conformer has been reached
                cc += [x] # storing the stable confermer in the cc list  
            elif theory in line:     # this indicats that the summary of the output file has been reached, 
                qm += line
                for line in sf : 
                    if 'RMSD=' in line :
                        qm += line.replace('RMSD=',theory)
                        break
                    else :
                        qm += line
                # now all the QM energies are collectes in one string 
                qm = qm.replace(' ','').replace('\n','').replace('|','').replace('\\','').replace('\r','')
                qm = qm.split(theory)[1].split(',')  
                qm = [float(q) for q in qm]  # this is the list of QM energies ofthe current PES file 
                QME += qm   # appending the QM energies of the current PES file to QME list
                # turning to the next PES file if any
                break
        sf.close()
    nconf = len(QME)
    cc = {i+1:cc[i] for i in range(nconf)}   # converting the cc list to dictionary such that the keys are the indices of the conformers 

    #__________________________________________________________________________________
    ### storing the values of the torsional angles for each conformer in a dictionary
    #   to avoid calculating the dihedral angles so many times during optimization
    #----------------------------------------------------------------------------------
    confdi = {}
    for x in cc : 
        c = cc[x]
        confdi[x] = {}
        for d in dihedrals :  # dihedrals is the dictionary built from the tmp file  
            [i,j,k,l] = dihedrals[d]['num']
            confdi[x][d] = dihedralvalue(c[i], c[j], c[k], c[l])
            # dihedralvalue is the function that calculates the values
            # of the dihedral angle through out the conformers given the cartesian coordinates of its atoms 


# preparing the list of QM energies
if qm_mm or scans : 
    nconf = len(QME)   # the number of conformers 
    minQM = min(QME)   # the minimum energy
    QME = [(q-minQM) for q in QME]  # shifting the QM energies such that the minimum is zero 
    if not qm_mm :    # QM energies introduced in the cfg are supposed to be in kcal/mole
        QME = [q*627.509 for q in QME]    # converting the QM energies from hartrees (gaussian unit) to kcal/mole



#_______________________________________________________________________________________________________
### extracting the torsional constraints and the fitting parameters from the confoiguration file (cnfg)
#-------------------------------------------------------------------------------------------------------
frozper = []    # the list of the frozen periodicities 
froztor = []    # the list of the frozen torsional parameters  
phases , kmax , maxenergy  = 0 , 8 , 12     # default values of torsional constraints, for more datails check below
wmethod , nsteps , tol = 0 , 5000 , 0.0001  # default values for the fitting 
wts = [1 for q in QME]  # initiating the weights dof the conformers to uniformity 


if cfg :    
    sf = open('cfg')    
    for line in sf :
        if not line.startswith('#') :
            
            line = line.split()
            n = len(line)
            if 'fp' in line and n > 5 :
                line = line[1:]
                frozper += [[int(p) for p in line]] # the first 4 entries after fp are the  indexes of the dihedral angle, the rest are the desired periodicities  
            elif 'ft' in line and n == 8:   # each line gives information on a dihedral angle with the desired parameters,
                                            # another (torsional parameters) for the same dihedral angle needs another line  
                line = line[1:]
                tor = [int(i) for i in line[:4]] # the indexes of the dihedral angle  
                tor += [float(line[4])]          # the desired torsional force constant for that dihedral angle 
                tor += [int(line[5])]            # the desired torsional periodicity for that dihedral angle
                tor += [float(line[6])]          # the desired phase shift for that dihedral angle
                froztor += [tor]
            elif 'phases' in line :              # default phases is not mentioned in the cfg file which means pahses = 0 (i.e. phase shaifts takes either 0 or 180) and positive force constants 
                phases = 1                       # if so phase shifts takes any values from -180 to 180 and the force constants are allowed to be positive or negative 
            elif 'kmax' in line and n == 2:      
                kmax = float(line[1])            # the maximum value allowed for the torsional force constants     
            elif 'maxenergy' in line and n == 2: 
                maxenergy = int(float(line[1]))  # the maximum value allowed for QM energies, conformers whose energies above this vale are excluded from fitiing
            elif 'nsteps' in line and n == 2:   
                nsteps = int(float(line[1]))     # the number of least squares iterations, i.e. fitting trials will stop at the assigned number even if satisfactory results are not reached  
            elif 'weight' in line :
                wmethod = line[1]                # there are two weighting methods; b for boltzmann and c for cutoff energy 
                if wmethod == 'b' and n == 3:
                    temp = float(line[2])        # the temperature which will be used for boltzmann distribution (wi = -Ei/kT)
                    wts = [exp(-1.0*q/(0.001987*temp)) for q in QME]
                elif wmethod == 'c' and n == 4: 
                    cutoff = float(line[2])      # the energy below which the conformers will be multiplied by a wighting factor   
                    wt = float(line[3])          # the weighting factor
                    wts = []
                    for q in QME :
                        if q > cutoff :
                            wts += [1]
                        else :
                            wts += [wt]                    
    # sorting froztor and frozper for each dihedral 
    sf.close()       


# assigning zero weights for the conformers whose energies > maxenergy 
for i in range(nconf) :
    if QME[i] > maxenergy :
        wts[i] = 0 
# the fractional weights
wtsum = float(sum(wts))

wts = [w/wtsum for w in wts]
print '\nThe following weights are in effect'
print wts
print


#__________________________________________________________________________________
### using the entries of the configuration file to make the related modifications
#----------------------------------------------------------------------------------
### constraints on the values of the torsional parameters
if phases == 1 :   # the force constants will range from 0 to kmax and the phase shifts will take any value from -180 to 180 during fitting      
    kmin = 0
    deltamax = 180.0
    deltamin = -180.0
    print 'The force constants will range between 0 and ' , kmax
    print 'The phase shifte will range between -180 and 180 degrees '
else:              # the force constants will range from -kmax to kmax and the phase shifts will be limited to either 0 or 180 after fitting  
    kmin = -kmax + 0.0 
    deltamax = 0
    deltamin = 0
    print 'The force constants will range from 0 to ' , kmax
    print 'The phase shifte will be either 0 or 180 degrees'


# assigning the frozen dihedrals and the frozen periodicities
# dihedrals dictionary will be filled with the frozen force constants and periodicities
# assumed we have froztor and frozper their structure will be something like 
#dihedrals =  {1: {'num': [2, 1, 5, 6]}, 2: {'num': [2, 1, 5, 7]}, ...}
#frozper =  [[2, 1, 5, 6, 1, 3], [2, 1, 5, 7, 2]]
#froztor =  [[2, 1, 5, 6, 3.3, 2, 0.0], [2, 1, 5, 6, 4.1, 1, 180.0], [2, 1, 5, 7, 0.1, 3, 0.0]]
for d in dihedrals :
    dihedrals[d]['f'] = []   # initiating force constant for each dihedral
    dihedrals[d]['p'] = []   # initiating periodicities  for each dihedral
    t1 = dihedrals[d]['num']
    pool = []                # to avoid the repeted periodicities
    for f in froztor :
        t2 = f[:4]
        t3 = [t1[3],t1[2],t1[1],t1[0]]        # that's the dihedral angle with its reverse to ensure natching  
        if t1 in [t2,t3]:
            if not f[5] in pool :             # this ensures a new periodicity (even for the same dihedral) (repeted periodicity can be introduced by the user by mistake)
                pool += [f[5]]
                dihedrals[d]['f'] += [f[4:]]  # appending the torsional parameters to the matched dihedral (more than one set of torsional parameters are allowed for the one dihedral)
    dihedrals[d]['f'].sort()                                   
    for f in frozper :                        # it is even allowed to introduce more periodicities than those mentioned in froztor (doublicate periodicities will be excluded)       
        t2 = f[:4]
        t3 = [t1[3],t1[2],t1[1],t1[0]]        # that's the dihedral angle with its reverse to ensure matching
        if t1 in [t2,t3]:
            for p in f[4:] :
                if not p in pool :
                    pool += [p]
                    dihedrals[d]['p'] += [p]         # we only accept the periodicities taht are not assigned in the frozen torsionals, doublicate periodicities are excluded 
    dihedrals[d]['p'].sort()


# the following would be the shape of dihedrals after appending the frozen torsional parameters and the frozen periodicities (without doublicates in the periodicities) 


'''
dihedrals = {1: {'p': [3], 'num': [2, 1, 5, 6], 'f': [[3.1, 1, 0.0], [3.3, 2, 180.0]]},
             2: {'p': [2], 'num': [2, 1, 5, 7], 'f': [[0.1, 3, 0.0]]},
             3 {'p': [], 'num': [3, 1, 5, 6], 'f': [[2.2, 3, 0.0]]},
             ...}
'''




#______________________________________________________________
### equivalencing the dihedrals with correcting the conflicts
#--------------------------------------------------------------
# equivalencing the dihedrals 
def equitorsional (t1,t2) : 
    k = 0
    uniqslist = [[i]+ uniqs[i] for i in uniqs]
    for i in range(4) : 
        for u in uniqslist : 
            if t1[i] in u and t2[i] in u : 
                k += 1 
                break
            elif t1[i] in u and t2[i-4] in u :    # to account for the reverse order
                k += 1 
                break
    if k == 4 : 
        return 1

equidihedrals = {d:[] for d in dihedrals}
tmp = [d for d in dihedrals]
pool = []
for i in tmp : 
    if i in dihedrals : 
        equidihedrals[i] += [i]
        if not i in pool :
            for j in tmp : 
                if j > i and not j in pool: 
                    if equitorsional (dihedrals[i]['num'],dihedrals[j]['num']) :
                        equidihedrals[i] += [j]
                        pool += [j]

print '\nThe program equivalencing scheme has detected the following equivalent torsional angles'
for d in equidihedrals :
    if len(equidihedrals[d]) > 1 : 
        print [dihedrals[q]['num'] for q in equidihedrals[d]]


# eliminating the conflicting frozen torsionals and periodicitiess from equivalent dihedrals
# if the user introduced different forozen torsional parameters for two dihedrals
#     while being equivalent according to onclick both dihedrals will be no longer equivqlent
for i in equidihedrals :                                     # the ith dihedral
    for j in equidihedrals[i] :                              # the jth equivalent dihedral of the ith dihedral
        if j > i :                                           # to exclude the ith dihedral itself 
            if dihedrals[i]['f'] and dihedrals[j]['f'] :     # if both the ith and jth dihedrals have frozen torsional parameters               
                if dihedrals[i]['f'] != dihedrals[j]['f'] :  # and if their frozen torsional parameters are different while being equivalent according the onclick
                    ks = [k for k in equidihedrals[i] if not k in [i,j]]  # i and j will be excluded from equivalence and the remaining dihedrals will be rearranged 
                    k = min(ks)                              # determining the minimum index of ks  
                    equidihedrals[k] = ks
                    equidihedrals[i] = [i]
                    equidihedrals[j] = [j]
            elif dihedrals[i]['p'] and dihedrals[j]['p'] :
                if dihedrals[i]['p'] != dihedrals[j]['p'] :
                    ks = [k for k in equidihedrals[i] if not k in [i,j]]
                    k = min(ks)
                    equidihedrals[k] = ks
                    equidihedrals[i] = [i]
                    equidihedrals[j] = [j]
    

# joining the equivalent dihedrals in order to freeze them all if one of them is frozen 
for i in equidihedrals :
    dihedrals[i]['q'] = equidihedrals[i]
    for j in equidihedrals[i] :
        if j > i : 
            if dihedrals[j]['f'] :       # if one of the dihedrals which are equivalent to the ith dihedral is frozen then its frozen parameters will be copied in the ith dihedral  
                dihedrals[i]['f'] = dihedrals[j]['f']
                break
            elif dihedrals[j]['p'] and not dihedrals[i]['f'] :
                dihedrals[i]['p'] = dihedrals[j]['p']
                break
    dihedrals[i]['w'] = len(dihedrals[i]['p'])

print '\nConsidering both the frozen torsionals and periodicities entries and the equivalencing scheme,'
print 'The equivalent dihedrals are modified to be'
for d in equidihedrals :
    if len(equidihedrals[d]) > 1 : 
        print [dihedrals[q]['num'] for q in equidihedrals[d]]


print 
print '\nThe frozen torsional parameters are modified to be '
for d in dihedrals :    
    if dihedrals[d]['f'] :
        for q in dihedrals[d]['q'] :
            num = dihedrals[q]['num']
            for f in dihedrals[d]['f'] :
                print num[0] , num[1] , num[2], num[3], ' ' , f[0] , f[1] , f[2] 


print '\nThe frozen periodicities are modified to be '
for d in dihedrals :    
    if dihedrals[d]['p'] :
        for q in dihedrals[d]['q'] :
            num = dihedrals[q]['num']
            print num[0] , num[1] , num[2], num[3], ' ' ,
            for p in dihedrals[d]['p'] :
                print p ,
            print 
 


#_____________________________________________________________
### calculating the Molecular mechanics energy surface (MME)
#-------------------------------------------------------------
if not qm_mm :
    #___________________________________________________________________________________________________________________________________________________________

    #        _____________________________________

    #           Molecular Mechanics Energy (MME) 
    #       ______________________________________

    #___________________________________________________________________________________________
    ### extracting the parameters from the tempar file (torsional paramaeters are not included) 
    #-------------------------------------------------------------------------------------------
    sf = open('tempars')
    b , u , a , m , j = 0 , 0 , 0 , 0 , 0 
    bondpars , ureypars , anglepars , improperpars , charge , LJ = {} , {} , {} , {} , {} , {}
    for line in sf :
        if 'bond' in line :
            b += 1
            bondpars[b] = {}
            line = line.split()
            bondpars[b]['num'] = [int(float(line[1])) , int(float(line[2]))]
            bondpars[b]['f'] = float(line[3])
            bondpars[b]['q'] = float(line[4])
        elif 'urey' in line :
            u += 1
            ureypars[u] = {}
            line = line.split()
            ureypars[u]['num'] = [int(float(line[1])) , int(float(line[2]))]
            ureypars[u]['f'] = float(line[3])
            ureypars[u]['q'] = float(line[4])
        elif 'angle' in line :
            a += 1
            anglepars[a] = {}
            line = line.split()
            anglepars[a]['num'] = [int(float(line[1])) , int(float(line[2])) , int(float(line[3]))]
            anglepars[a]['f'] = float(line[4])
            anglepars[a]['q'] = float(line[5])
        elif 'improper' in line :
            m += 1
            improperpars[m] = {}
            line = line.split()
            improperpars[m]['num'] = [int(float(line[1])) , int(float(line[2])) , int(float(line[3])) , int(float(line[4]))]
            improperpars[m]['f'] = float(line[5])
            improperpars[m]['q'] = float(line[6])
        elif 'lj' in line :
            j += 1
            line = line.split()
            LJ[j] = [int(float(i)) for i in line[1:]]
        elif 'charges' in line :
            line = line.split()
            charge =  {i:float(line[i]) for i in range(1,len(line))}


    #_________________________________________________________________________
    ### calculating all types of MME energies except the unfrozen torsionals  
    #-------------------------------------------------------------------------
    counter = 0
    MME = []
    for s in cc : 
        c = cc[s]
        counter += 1
        ene = 0
        
        for b in bondpars : 
            [i,j] = bondpars[b]['num']
            d = dist(c[i],c[j])
            f = bondpars[b]['f']
            q = bondpars[b]['q']
            ene += f*(d-q)**2
            
        for u in  ureypars: 
            [i,j] = ureypars[u]['num']
            d = dist(c[i],c[j])
            f = ureypars[u]['f']
            q = ureypars[u]['q']
            ene += f*(d-q)**2
            
        for a in anglepars : 
            [i,j,k] = anglepars[a]['num']
            p = anglevalue(c[i],c[j],c[k])
            f = anglepars[a]['f']
            q = anglepars[a]['q']
            ene += f*((p-q)*pi/180)**2
            
        for d in dihedrals :
            if dihedrals[d]['f'] :            
                for q in dihedrals[d]['q'] : 
                    [i,j,k,l] = dihedrals[q]['num']                
                    t = dihedralvalue(c[i],c[j],c[k],c[l])
                    for h in dihedrals[d]['f'] :           
                        [f,p,delta] = h
                        ene += f*(1 + cos((p*t-delta)*pi/180))
                        
        for m in improperpars : 
            [i,j,k,l] = improperpars[m]['num']
            d = dihedralvalue(c[i],c[j],c[k],c[l])
            f = improperpars[m]['f']
            q = improperpars[m]['q']
            ene += f*((d-q)*pi/180)**2
            
        for p in pairs : 
            [i,j] = p
            r = dist(c[i],c[j])
            ene += 332.0716 * charge[i] * charge[j] / r
            epsi , ri = LJ[i][0] , LJ[i][1]
            epsj , rj = LJ[j][0] , LJ[j][1]
            eps = sqrt(epsi*epsj)
            rmin = (ri+rj)
            ene += eps * ((rmin/r)**12 -2*(rmin/r)**6)
            
        for p in pairs4 : 
            [i,j] = p
            r = dist(c[i],c[j])
            ene += 332.0716 * charge[i] * charge[j] / r
            epsi , ri = LJ[i][2] , LJ[i][3]
            epsj , rj = LJ[j][2] , LJ[j][3]
            eps = sqrt(epsi*epsj)
            rmin = (ri+rj)
            ene += eps * ((rmin/r)**12 -2*(rmin/r)**6)  
        MME += [ene]

if qm_mm or  scans : 
	minMM = min(MME)
	MME = [m-minMM for m in MME]


#_____________________________________________________________
### the torsional energies which are required to be fitted 
#------------------------------------------------------------
Ts = [QME[i]-MME[i] for i in range(nconf)]
m = min(Ts)
Ts = [t-m for t in Ts]

'''
print 'dihedrals (the dictionary of the dihedrals after extracting the frozen parameters and considering the equivalent duhedrals) ' 
for d in dihedrals : 
    print d , dihedrals[d]
'''


#____________________________________________________________________________________
### preparing the dihedrals which will be introduced to the least squares algorithm 
#------------------------------------------------------------------------------------
### reducing the equivalent dihedrals[d]['q'] by emptying it of its contents  
for d in equidihedrals :
    for q in dihedrals[d]['q'] :
        if q != d :
            dihedrals[q] = {'num':dihedrals[q]['num'], 'q':[], 'p':[], 'f':[], 'w':0}


### reducing the frozen torsionals by emptying it of its contents
for d in equidihedrals :
    if dihedrals[d]['f'] :
        dihedrals[d] = {'num':dihedrals[d]['num'], 'q':[], 'p':[], 'f':1, 'w':0}


### assigning the weight of dihedrals to be fit (that is the length of the equivalent dihedrals, this significantly decreases the time of fitting)
# a list will be added to be filled in with the fitting torsional parameters 
for d in dihedrals :
    if dihedrals[d]['q'] :
        dihedrals[d]['w'] = len(dihedrals[d]['q'])
        dihedrals[d]['t'] = []



### adding the missing periodicities
# in fact we need here to find out a reliable scheme to assign the periodicities for a dihedral angle
for d in dihedrals :
    if dihedrals[d]['q'] and not dihedrals[d]['p'] : 
        [i,j] = dihedrals[d]['num'][1:3]
        ni , nj = len(atoms[i]['cona']) , len(atoms[j]['cona'])
        if sorted([ni,nj]) == [3,4] : 
            dihedrals[d]['p'] = [1,3,6]
        else : 
            dihedrals[d]['p'] = [max(ni,nj)-1]



print '\nThe following periodicities will be in effect'
for d in dihedrals :
    if dihedrals[d]['p'] : 
        print dihedrals[d]['num'] + dihedrals[d]['p']



dfit = {}    # the dihedrals subject to fitting 
for d in dihedrals :
    if dihedrals[d]['w'] :
        dfit[d] = {'q': dihedrals[d]['q'] , 'p':dihedrals[d]['p']}



#______________________________________________________________
### excluding the conformers whose energies exceeds maxenergy
#--------------------------------------------------------------

QME2 , MME2 , Ts2 , wts2  , nconf2 = QME+[] , MME+[] , Ts+[] , wts+[] , len(QME)
QME , MME , Ts , wts = [], [], [], []
for i in confdi : 
    c = i-1
    if QME2[c] <= maxenergy : 
        QME += [QME2[c]+0.0]
        MME += [MME2[c]+0.0]
        Ts += [Ts2[c]+0.0]
        wts += [wts2[c]+0.0] 
nconf = len(QME)


#________________________________________________________
### building the matrix A and the vector V for A*K = V 
#--------------------------------------------------------
V = Ts + []
'''
let the torsional energy of the cth conformer (i.e. Ts[c]) be V_c
the least square fitting will work on

V_c = sum_d sum_n K_dn*(1 + cos(n*A_cd - delta_dn)) 
K_dn is the force constant of the dth dihedral with the nth periodicity
delta_dn is the phase shift of the dth dihedral with the nth periodicity
A_cd is the value of the dth dihedral in the cth conformer

dropping the 1 will not affect the fitting
let delta_dn = 0 because we will only consider the values 0 and 180 for the phase shift (notice: cos(180) = -cos(0) ), then

V_c = sum_d sum_n K_dn*(cos(n*A_cd)

'''

## building A
# spreading the dihedrals with periodicities of all conformers in a matrix
# the equivalent dihedrals will summed to give one entry for each conformer
A = []
i = 0
for c in confdi :
    if QME2[c-1] <= maxenergy : 
        A += [[]]
        for d in dfit :
            for p in dfit[d]['p'] :
                a = 0.0
                for q in dfit[d]['q'] :          # summing the equivalent dihedrals into one force constant         
                    a += cos(p*confdi[c][q]*pi/180.0) 
                A[i] += [a]
        i += 1


print len(A) , len(A[0]) , '=================================================='

#________________________________________________________
### Modifying the matrices for  A'A*K = A'V 
#--------------------------------------------------------
TA = T(A)
AA = MM(TA,A)
AV = Mv(TA,V)


'''
print 'QME = ' , QME 
print 
print 'MME = ' , MME 
print
print 'Ts = ' , Ts 
print
print 'confdi (dihedral values over the conformers) = '
for d in confdi : 
    print d , confdi[d]
print 

print '\nthe dihedrals suject to fitting ' 
for d in dfit :
    print d , dfit[d]

print 
print 'the A matrix built from dfit and confdi' 
for a in A : 
    print a 
print len(A) , len(A[0])
'''


'''
TA = T(A)
AA = MM(TA,A)
AV = Mv(TA,V)
s = 0.05

Kf = LSQR(AA,AV) 
V2 = Mv(A,Kf)



print '============================='
print 'The final force constants are '
for k in Kf : 
    print k
print '============================='
MME = [MME[i]+V2[i] for i in range(nconf)]
m = min(MME)
MME = [M-m for M in MME]

for i in range(nconf) : 
    print i+1 , QME[i] , MME[i]

'''



# function to add the implement the restraint to the equation AA*K = AV
def restraint (AA,AV,Kf,step,s): 
    nK = len(Kf)    
    res = []
    for i in range(nK) : 
        if step == 0 : 
            res += [s]
        else : 
            res += [s/sqrt(Kf[i]**2+0.01)]
    AAr = [[AA[i][j] for j in range(nK)] for i in range(nK)]
    AVr = [AV[i] for i in range(nK)]
    for i in range(nK) :    
        AAr[i][i] += res[i]
        if step == 0 : 
            AVr[i] += res[i]*Kf[i]
    return AAr , AVr 



s = 0.05
## applying the least square fitting over and over to reach the optimized parameters 
Ki = [0 for i in AA]
nK = len(Ki)
Ko =  [i for i in Ki]
Kf =  [i for i in Ki]
for step in range(10) :
    print step

    #AA , AV =  restraint (AA , AV , Kf ,step,s)        
    Kf = LSQR(AA,AV)  
    V2 = Mv(A,Kf)

    residual = sqrt(sum([wts[c]*(V[c]-V2[c])**2 for c in range(nconf)]))
    conv1 = sqrt(sum([(Ko[i]-Kf[i])**2 for i in range(nK)])/nK)
    conv2 = sqrt(sum([wts[c]*(V[c]-V2[c])**2 for c in range(nconf)]))
    if  conv1 < 1e-5 and conv2 < 0.01:
        print 'least squares fitting has converged in  '  , step , ' iterations'       
        break
    elif step == nsteps:
        print 'no convergence after '  , step , ' iterations'
        break  
    else :
        ko = [k for k in Kf]



print 
for k in Kf : 
    print k
print 

m = min(V2)
V2 = [v-m for v in V2]

for i in range(len(V)) : 
    print i+1 , V[i] , V2[i]



print 

MME = [MME[c] + V2[c] for c in range(nconf)]
m = min(MME)
MME = [M-m for M in MME]

for c in range(nconf) : 
    print c+1 , QME[c] , MME[c]





'''
# function to calculate torsional energies after fitting the force constants (V2) 
def tore(A,Kf,V) : 
    nconf = len(V)
    
    # verifying the obtained force constants and phase shifts
    # as we have the Ts which were calculated as the difference between QME and MME (Ts = QME - MME)
    # and as obtained by the least squares fitting we have the Kf which is supposed to be the optimized force constants
    # then we calculates Tf which is supposed to make Ts - Tf minimum
    # Tf is calculated from scratch, i.e. the dihedrals dictionary will be filled by the fit torsional parameters snd then the Tf is calculated normally 
    # here we calculate thr torsional energy using the fit force constants and subtract it from the torsional energy   
    # we will first multiply A*Kf and see whether it works or not 
    # if it doesn't work we will use the function V_c = sum_d sum_n [K_dn (  1 + cos(n phi_cd + delta_dn)  ) ] in dfits 
    # if both work with the same results we may use them both together

    V2 = Mv(A,Kf)
    m = min(V2)
    V2 = [v-m for v in V2]
    m = min(V)
    V = [v-m for v in V]

    
    # this part wil be used in case the the previouis one doesn't work  
    nK = len(Kf)
    delta = [0.0 for i in range(nK)]
    # limiting K[i] to positive values and making delta[i] = 180 if required 
    for i in range(nK):
        if Kf[i] < 0  :
            Kf[i] = - Kf[i]
            delta[i] = 180.0
     
    # distributing Kf and delta over the dihedrals taking the periodicities into account 
    # this process is supposed to be the reverse of the building of A     
    i = 0 
    for d in dfit : 
        for p in dfit[d]['p'] : 
            dfit[d]['t'] += [[Kf[i] , p , delta[i]]]     # we'd better fill in dihedrals like that to make it easy when we collect all torsional parameters  
        i += 1

    # calculating the torsional energy normally using the torsional energy function 
    V2 = []
    for c in confdi : 
        for d in dfit : 
            for q in dfit[d]['q'] : 
                phi = confdi[c][q]*pi/180
                for t in dfit[d]['t'] :                     
                    [k , p ,  theta] = t 
                    V2 += [k* (1 + cos(p*phi-theta) )]
    # m = min(V2)
    # V2 = [V2[c]-m for i in range(nconf)]   # do we really need to set the minimum to zero here or later (i.e. after calculating the final MME) ????????
    residual = sqrt(sum([wts[c]*(V[c]-V2[c])**2 for c in range(nconf)]))
    return V2 , dfit , residual  



## applying the least square fitting over and over to reach the optimized parameters 
Ki = [0 for i in AA]
nK = len(Ki)
Ko =  [i for i in Ki]
Kf =  [i for i in Ki]
for step in range(5001) :
    print step
    AA , AV =  restraint (AA , AV , Kf ,step,s)        
    Kf = LSQR(AA,AV)  
    V2 = Mv(AA,Kf)
    print len(wts) , len(V), len(V2), nconf
    residual = sqrt(sum([wts[c]*(V[c]-V2[c])**2 for c in range(nconf)]))

    #V2 , dfit , residual = tore(AA,Kf,AV)              

    conv1 = sqrt(sum([(Ko[i]-Kf[i])**2 for i in range(nK)]))/nK        
    conv2 = sqrt((sum([(qme[i]-MME[i])**2 for i in range(nconf)]))/nconf)
    if  conv1 < 1e-5 and conv2 < 0.01:
        print 'least squares fitting has converged in  '  , step+1 , ' iterations'       
        break
    elif step == 5000 :
        print 'no convergence after '  , step+1 , ' iterations'
        break  
    else :
        fo = [f for f in ff]

    if residual < 0.01 : 
        print 'least squares fitting has converged in  '  , step+1 , ' iterations'       
        break 
    elif step == 5000 :
        print 'no convergence after '  , step+1 , ' iterations'
    else : 
        Ko = [f for f in Kf]

if 1 :
    MME = [MME[c]+V2[c] for c in range(nconf)]
    m = min(MME)
    MME = [M-m for M in MME]
    for c in range (nconf) : 
        print QME[c] , MME[c]



if 1 :     
    # this part wil be used in case the the previous one doesn't work  
    nK = len(Kf)
    delta = [0.0 for i in range(nK)]
    # limiting K[i] to positive values and making delta[i] = 180 if required 
    for i in range(nK):
        if Kf[i] < 0  :
            Kf[i] = - Kf[i]
            delta[i] = 180.0
     
    # distributing Kf and delta over the dihedrals taking the periodicities into account 
    # this process is supposed to be the reverse of the building of A     
    i = 0 
    for d in dfit : 
        for p in dfit[d]['p'] : 
            dfit[d]['t'] += [[Kf[i] , p , delta[i]]]     # we'd better fill in dihedrals like that to make it easy when we collect all torsional parameters  
        i += 1

    # calculating the torsional energy normally using the torsional energy function 
    V2 = []
    for c in confdi : 
        for d in dfit : 
            for q in dfit[d]['q'] : 
                phi = confdi[c][q]*pi/180
                for t in dfit[d]['t'] :                     
                    [k , p ,  theta] = t 
                    V2 += [k* (1 + cos(p*phi-theta) )]
    # m = min(V2)
    # V2 = [V2[c]-m for i in range(nconf)]   # do we really need to set the minimum to zero here or later (i.e. after calculating the final MME) ????????

MME = [MME[c]+V2[c] for c in range(nconf)]
m = min(MME)
MME = [M-m for M in MME]
for c in range (nconf) : 
    print QME[c] , MME[c]


#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
### keeping the force constant positive
if not phases :
    for d in fitf : 
        for p in fitf[d] : 
            if fitf[d][p][0] < 0 : 
                fitf[d][p][0] = -fitf[d][p][0]
                fitf[d][p][2] = 180


### adding the fitted torsional energies to the MME energies
if qm_mm or scans : 
	MMf = [MME[i]+TEf[i] for i in range(nconf)]
	minMMf = min(MMf)
	MMf = [MMf[i]-minMMf for i in range(nconf)]






## copying the 

print '\nThe torsional parameters (atom1 atom2 atom3 atom4  force constant    multiplicity     phase angle)'
tf = open('dihedrals','w')
for d in fitf : 
    for q in dihedrals[d]['q'] : 
        for p in fitf[d] : 
            fitf[d][p][0] = round(fitf[d][p][0],4)
            fitf[d][p][1] = int(fitf[d][p][1])
            #tf.write('dihedral ')
            for a in dihedrals[q]['num']: 
                print a , '  ',
                tf.write(str(a) + '   ')
            for a in fitf[d][p] : 
                print a , '  ', 
                tf.write(str(a) + '   ')
            print
            tf.write('\n')
tf.close()

print '\nThe Quantum mechanical and the Molecular mechanical energies '
for i in range(nconf) :
    print '%3s %15s %15s'% (i+1 , QME[i] , MMf[i])


print
print 'Torsional Fitting is done  _______________________________________________________'
print '__________________________________________________________________________________'


tf = open('pesdone','w')
tf.close()

'''