### Monte Carlo Simulated Annealing
### 12-1-2018
import sys
from math import *
from random import random , uniform
from functions import *

## entries just for trials until we integrate it mcsa with other parts
# the temporary files temp and tempars gave to be kept to run this script 

print ( '\n\n' )
print ( '__________________________________________________________________________________' )
print ( 'Torsional Fitting ________________________________________________________________' )




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
dvals = {i:{} for i in range(1,999999)}   # 

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
                    nconf = n + 0
                elif 'mm' in line :
                    MME = [float(m) for m in line[1:]]
                    if len(MME) != n :     # invalid entries; unequal lengths of QM and MM energies, then
                        qm_mm = 0          # qm_mm enquiry will be canceled and
                        break              # no need to go further in the cfg file  
                elif 'dihedral' in line :
                    if len(line)-5 != n :  # invalid entries; unequal lengths of dihedral values and QM energies, then
                        qm_mm = 0          # qm_mm enquiry will be canceled and  
                        break              # no need to go further in the cfg file 
                    d += 1     # counter for the number of dihedrals 
                    dihedrals[d] = {}
                    dihedrals[d]['num'] = [int(float(i)) for i in line[1:5]]
                    values =  [float(i) for i in line[5:]]
                    for i in range(len(values)) :
                        v = values[i] 
                        dvals[i+1][d] = v 
    sf.close()
    dvals = {i:dvals[i] for i in dvals if dvals[i]}

# if one of QME, MME or dvals is empty (means the keyword is absent in the cfg file)
# qm_mm will be ignored and the program will exit
if not (QME and MME and dvals) :
    qm_mm = 0


if qm_mm == 0 and claimed == 1 :
    print ( '\nYou have limited the parameterization to torsional fitting ' )
    print ( 'but you did not provide valid entries in the configuration file (.cfg)' )
    df = open('dihedrals','w')
    df.close()
    cfg = 0
    scans , cc , QME , MME , dvals = [] , {} , [] , [] , {}
    nconf = 0
    sys.exit()


scans = []
if not qm_mm : 
    scans = sys.argv[1:-1]
    #scans = ['nma1.scn']
    n = len(scans)
    print ( '\nThe folowing {} PES files were detected'.format(n) )
    for s in scans : 
        print ( s + '  ' , end=' ')
    print ( '\n' )



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
    dvals = {}
    for x in cc : 
        c = cc[x]
        dvals[x] = {}
        for d in dihedrals :  # dihedrals is the dictionary built from the tmp file  
            [i,j,k,l] = dihedrals[d]['num']
            dvals[x][d] = dihedralvalue(c[i], c[j], c[k], c[l])
            # dihedralvalue is the function that calculates the values
            # of the dihedral angle given the cartesian coordinates of its atoms 

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
phases , kmax , maxenergy  = 0 , 5 , 10     # default values of torsional constraints, for more datails check below
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
    sf.close()       

# assigning zero weights for the conformers whose energies > maxenergy 
for i in range(nconf) :
    if QME[i] > maxenergy :
        wts[i] = 0 
# the fractional weights
wtsum = float(sum(wts))

wts = [w/wtsum for w in wts]
print ( '\nThe following weights were in effect' )
print ( wts )
print ()



#__________________________________________________________________________________
### using the entries of the configuration file to make the related modifications
#----------------------------------------------------------------------------------
### constraints on the values of the torsional parameters
if phases == 1 :   # the force constants will range from 0 to kmax and the phase shifts will take any value from -180 to 180 during fitting      
    kmin = 0
    deltamax = 180.0
    deltamin = -180.0
    print ( 'The force constants ranged between 0 and ' , kmax )
    print ( 'The phase shifts ranged between -180 and 180 degrees ' )
else:              # the force constants will range from -kmax to kmax and the phase shifts will be limited to either 0 or 180 after fitting  
    kmin = -kmax+0.0 
    deltamax = 0
    deltamin = 0
    print ( 'The force constants ranged from 0 to ' , kmax )
    print ( 'The phase shifts were either 0 or 180 degrees' )


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


# this would be the shape of dihedrals after appending the frozen torsional parameters and the frozen periodicities (without doublicates in the periodicities) 
'''
dihedrals = {1: {'p': [3], 'num': [2, 1, 5, 6], 'f': [[3.3, 2, 0.0], [4.1, 1, 180.0]]},
             2: {'p': [2], 'num': [2, 1, 5, 7], 'f': [[0.1, 3, 0.0]]},
             3: {'p': [], 'num': [3, 1, 5, 6], 'f': []},
             ...}
'''


#____________________________________________________________________
### equivalencing the dihedrals with correcting the conflicts if any
#--------------------------------------------------------------------
# equivalencing the dihedrals
print ( uniqs )
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

print ( '\nThe program equivalencing scheme has detected the following equivalent torsional angles' )
for d in equidihedrals :
    if len(equidihedrals[d]) > 1 : 
        print ( [dihedrals[q]['num'] for q in equidihedrals[d]] )

# eliminating the conflicting frozen torsionals and periodicitiess from equivalent dihedrals
for i in equidihedrals :
    for j in equidihedrals[i] :
        if j > i :
            if dihedrals[i]['f'] and dihedrals[j]['f'] :
                if dihedrals[i]['f'] != dihedrals[j]['f'] :
                    print ( i , j )
                    ks = [k for k in equidihedrals[i] if not k in [i,j]]
                    k = min(ks)
                    equidihedrals[k] = ks
                    equidihedrals[i] = [i]
            elif dihedrals[i]['p'] and dihedrals[j]['p'] :
                if dihedrals[i]['p'] != dihedrals[j]['p'] :
                    print ( i , j )
                    ks = [k for k in equidihedrals[i] if not k in [i,j]]
                    k = min(ks)
                    equidihedrals[k] = ks
                    equidihedrals[i] = [i]
    

# joining the equivalent dihedrals
for i in equidihedrals :
    dihedrals[i]['q'] = equidihedrals[i]
    for j in equidihedrals[i] :
        if j > i : 
            if dihedrals[j]['f'] :
                dihedrals[i]['f'] = dihedrals[j]['f']
                break
            elif dihedrals[j]['p'] and not dihedrals[i]['f'] :
                dihedrals[i]['p'] = dihedrals[j]['p']
                break
    dihedrals[i]['w'] = len(dihedrals[i]['p'])

print ( '\nTo consider the frozen torsionals and periodicities entries in the configuration file,' )
print ( 'the following equivalent torsionals will be in effect' )
for d in equidihedrals :
    if len(equidihedrals[d]) > 1 : 
        print ( [dihedrals[q]['num'] for q in equidihedrals[d]] )



if not qm_mm :
    #___________________________________________________________________________________________________________________________________________________________

    #        ________________________________

    #           Molecular Mechanics Energy 
    #       __________________________________

    #_______________________________________________________________
    ### The lists of all parameters (non frozen torsional paramaeters are not included) 
    #---------------------------------------------------------------
    # the following MME parameters will be used to calculate the MME

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


    #_____________________________________________________________________
    ### calculating all types of MME energies except un frozen torsionals  
    #---------------------------------------------------------------------
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
        for a in improperpars : 
            [i,j,k,l] = improperpars[a]['num']
            d = dihedralvalue(c[i],c[j],c[k],c[l])
            f = improperpars[a]['f']
            q = improperpars[a]['q']
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
            rmin = (ri+rj)/2
            ene += eps * ((rmin/r)**12 -2*(rmin/r)**6)  
        MME += [ene]

if qm_mm or  scans : 
	minMM = min(MME)
	MME = [m-minMM for m in MME]


### the torsional energies which are required to be fit 
Ts = [QME[i]-MME[i] for i in range(nconf)]
 

### reducing the equivalent dihedrals[d]['q']
for d in equidihedrals :
    for q in dihedrals[d]['q'] :
        if q != d :
            dihedrals[q] = {'num':dihedrals[q]['num'], 'q':[], 'p':[], 'f':[], 'w':0}

### reducing the frozen torsionals
for d in equidihedrals :
    if dihedrals[d]['f'] :
        dihedrals[d] = {'num':dihedrals[d]['num'], 'q':[], 'p':[], 'f':1, 'w':0}


### adding the missing periodicities 
for d in dihedrals :
    if dihedrals[d]['q'] and not dihedrals[d]['p'] : 
        [i,j] = dihedrals[d]['num'][1:3]
        ni , nj = len(atoms[i]['cona']) , len(atoms[j]['cona'])
        if sorted([ni,nj]) == [3,4] : 
            dihedrals[d]['p'] = [1,3,6]
        else : 
            dihedrals[d]['p'] = [max(ni,nj)-1]

print ( '\nThe following periodicities will be in effect' )
for d in dihedrals :
    if dihedrals[d]['p'] : 
        print ( dihedrals[d]['num'] + dihedrals[d]['p'] )

#////////////////////////////////////////////////////////////////////////////////////////////////


#__________________________________________________________________________________________
### The Annealing process
#____________________________________________

# a function to generate a new point in the vicinity of the last accepted point during fitting steps 
def generatepoint(fitc,kmin,kmax,deltamin,deltamax):
    fitn = {}
    for d in fitc : 
        fitn[d] = {}
        for p in dihedrals[d]['p'] : 
            k = fitc[d][p][0] + uniform(-0.5,0.5)
            delta = fitc[d][p][2] +uniform(-10.0,10.0) 
            k = max(min(kmax,k),kmin)
            delta = max(min(deltamax,delta),deltamin)
            fitn[d][p] = [k,p,delta]
    return fitn

### a function to calculate the energies of the torsional angles
def fitenergy(fitn):   
    TF = []
    for v in dvals : 
        te = 0
        for d in fitn :  
            for p in fitn[d] : 
                k , per , delta = fitn[d][p]  
                for t in dihedrals[d]['q'] :
                    f = dvals[v][t]+0.0
                    te += k*(1 + cos((per * f - delta)*pi/180)) 
                    #print ( k , per , f , delta , cos((per * f - delta)*pi/180)  , te )
        TF += [te]
    return TF 

### a function to calculate the square root
def rmse(nconf,QME,MME,TE,wts,maxenergy):  
    TEav = 0
    for i in range (nconf) :
        TEav = TEav + wts[i] * TE[i]
    c = MMav - QMav + TEav 
    mse = 0
    for i in range(nconf) :
        if QME[i] < maxenergy :
            mse = mse + (wts[i]*(QME[i] - (MME[i] + TE[i]) + c)**2)
    rms = sqrt(mse)
    return rms

### a function to update the datapoints durung the MCSA
def update(fitn):
    fitc = {}
    for d in fitn : 
        fitc[d] = {} 
        for p in dihedrals[d]['p'] : 
            k = fitn[d][p][0]+0.0
            delta = fitn[d][p][2]+0.0
            fitc[d][p] = [k,p,delta]
    return fitc

QMav = 0
for i in range (nconf) :
    QMav = QMav + wts[i] * QME[i]

MMav = 0
for i in range (nconf) :    
    MMav = MMav + wts[i] * MME[i]

        
#_________________________
# Annealing parameters
#-------------------------
n = 100                          # Number of cycles
m = 50                           # Number of trials per cycle
na = 0.0                         # Number of accepted solutions
p1 = 0.99                         # Probability of accepting worse solution at the start
p50 = 0.00001                      # Probability of accepting worse solution at the end
t1 = -1.0/log(p1)                # Initial temperature
t50 = -1.0/log(p50)              # Final temperature
frac = (t50/t1)**(1.0/(n-1.0))   # Fractional reduction every cycle
na = na + 1.0
t = t1 + 0   # current temperature
dEav = 0.0
if nsteps : 
    n = int(2*nsteps/m)
print ( '\nThe anneling parameters' )
print ( '    The starting temperature = {} K'.format(t1) )
print ( '    The final temperature = {} K'.format(t50) )
print ( '    The cooling scheme: temp(n+1) =  temp(n)*' , frac )
print ()
#_________________________________
# Initialize the set of fitting parameters and the resulting variables
#---------------------------------
fitc , fitf = {} , {}
for d in dihedrals :
    if dihedrals[d]['q']:
        fitc[d] , fitf[d] = {} , {}         
        for p in dihedrals[d]['p'] :
            fitc[d][p] = [0.0,p,0.0]
            fitf[d][p] = [0.0,p,0.0]

TEc = [0 for e in QME]
TEf = [0 for e in QME]



print ( '\nRoot Mean Squred Error Accepted' )
rmsc = rmse(nconf,QME,MME,TEc,wts,maxenergy)
rmsf = rmsc + 0.0

print ( 'RMSE = ' , rmsf )

#________________________
# the iteration process
#------------------------
for i in range(n):
    for j in range(m):
        # Generate new trial points
        fitn = generatepoint(fitc,kmin,kmax,deltamin,deltamax)
        TEn = fitenergy(fitn) 
        rmsn = rmse(nconf,QME,MME,TEn,wts,maxenergy)
   
        #print ( xn , fn
        dE = abs(rmsn - rmsc)
        if (rmsn > rmsc): # it is a worse solution
            if (i==0 and j==0):
                dEav = dE + 0
            if dEav < 1e-15 :
                dEav = 1e-15
            p = exp(-dE/(dEav * t)) # the probability of acceptance
            if (random() < p):
                accepted = 1
            else:
                accepted = 0
        else:
            accepted = 1
        if accepted:                
            fitc = update(fitn)
            rmsc = rmsn + 0            
            na = na + 1.0  # increment number of accepted solutions            
            dEav = (dEav * (na-1.0) +  dE) / na # update DeltaE_avg 
        if rmsn < rmsf :   
            fitf = update(fitn)
            rmsf = rmsn + 0
            print ( 'RMSE = ' , rmsf )
            for i in range(nconf) :                    
                TEf[i] = TEn[i]+0
              
    t = frac * t  # Lower the temperature for next cycle


### keeping the force constant positive
if not phases :
    for d in fitf : 
        for p in fitf[d] : 
            if fitf[d][p][0] < 0 : 
                fitf[d][p][0] = -fitf[d][p][0]
                fitf[d][p][2] = 180


### adding the fit torsional energies to the MME energies
if qm_mm or scans : 
	MMf = [MME[i]+TEf[i] for i in range(nconf)]
	minMMf = min(MMf)
	MMf = [MMf[i]-minMMf for i in range(nconf)]


print ( '\nThe torsional parameters (atom1 atom2 atom3 atom4  force constant    multiplicity     phase angle)' )
tf = open('dihedrals','w')
for d in fitf : 
    for q in dihedrals[d]['q'] : 
        for p in fitf[d] : 
            fitf[d][p][0] = round(fitf[d][p][0],4)
            fitf[d][p][1] = int(fitf[d][p][1])
            #tf.write('dihedral ')
            for a in dihedrals[q]['num']: 
                print ( a , '  ', end=' ')
                tf.write(str(a) + '   ')
            for a in fitf[d][p] : 
                print ( a , '      ', end=' ')
                tf.write(str(a) + '   ')
            print ()
            tf.write('\n')
tf.close()

print ( '\nThe Quantum mechanical and the Molecular mechanical energies ' )
j = 0
for i in range(nconf) :
    if QME[i] < maxenergy : 
        j += 1
        print ( '%3s %15s %15s'% (j , QME[i] , MMf[i]) )


print ()
print ( 'Torsional Fitting is done  _______________________________________________________' )
print ( '__________________________________________________________________________________' )


tf = open('pesdone','w')
tf.close()
