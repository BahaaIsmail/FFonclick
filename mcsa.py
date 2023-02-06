### Monte Carlo Simulated Annealing
### 12-5-2018
import sys
from math import *
from random import random , uniform
from functions import *

print ('Torsional Parameterization .. >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

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
'''
alternative = sys.argv[-1]
if alternative == '1' :
    alternative = 1
else :
    alternative = 0
'''
alternative = 1
claimed = alternative + 0

if alternative : 
    #_______________________________________________________________
    ### extracting the QM and MM energies with the associated
    #   torsional options if existing in the configuration file
    #---------------------------------------------------------------
    sf = open('cfg')
    QME , MME , = [] , []
    ignored = []
    dihedrals = {}
    dvals = {i:{} for i in range(1,999999)}
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
                    if len(MME) != n :
                        MME = []
                        alternative = 0
                        break
                elif 'torsional' in line :
                    if len(line)-5 != n :
                        alternative = 0
                        break
                    d += 1
                    dihedrals[d] = {}
                    dihedrals[d]['a'] = [int(float(i)) for i in line[1:5]]
                    values =  [float(i) for i in line[5:]]
                    for i in range(len(values)) :
                        v = values[i] 
                        dvals[i+1][d] = v 
    sf.close()
    dvals = {i:dvals[i] for i in dvals if dvals[i]}
print (alternative)

if not (QME and MME and dvals) :
    alternative = 0


scans = 0
if not alternative : 
    #scans = sys.argv[1:-1]
    scans = ['nma1.scn' ,  'nma2.scn' , 'nma3.scn']
    print (scans)

if alternative == 0 and claimed == 1 :
    print ('\nYou have requested the alternative way for torsional parameterization ')
    print ('but you did not provide valid entries in the configuration file (.cfg)')
    print ('the PES files will be used if provided, otherwise the process will be halted')
    if not scans :
        sys.exit()



#_______________________________________________________________
# building the internal coordinates from the tmp file
#---------------------------------------------------------------
atoms , bonds = atoms_bonds()
uniqs = get_symmetries(atoms)

if not alternative :
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


if not alternative : 
    #_______________________________________________________________
    # extracting the level of theory
    #---------------------------------------------------------------
    sf = open(scans[0])
    theory = 'HF'
    for line in sf :
        if line.startswith(' # ') :
            if 'mp2' in line or 'MP2' in line :
                theory = 'MP2'
            elif 'b3' in line or 'B3' in line :
                theory = 'DFT'
            elif 'CCSD' in line or 'ccsd' in line :
                theory = 'CCSD'
            break
    sf.close()

    #_______________________________________________________________
    ### extracting the stable conformers (coordinates and energies)
    #---------------------------------------------------------------
    cc , QME = [] , []
    for f in scans :
        sf = open(f)
        for line in sf :
            if 'Standard orientation' in line :
                sf.readline() ; sf.readline() ; sf.readline() 
                line = sf.readline()
                x = {}
                for line in sf :
                    if line.startswith(' ----------------') :
                        break
                    else :
                        lin = line.split()
                        x[int(lin[0])] = [float(lin[i]) for i in range(3,6)]
            elif 'Optimized Parameters' in line:
                cc += [x]
            elif '1|1|' in line or '\1\1' in line: 
                theo = '|HF=' 
                if theory == 'MP2' :
                    theo = '|MP2='
                elif theory == 'CCSD' :
                    theo ='|CCSD=' 
                qm = '' 
                for line in sf : 
                    if '@' in line :
                        qm += line
                        break 
                    else :
                        qm += line
                qm = qm.replace('\r\n','').replace('\n','').replace(theo,'!').replace('|RMSD=','!')
                
                qm = qm.split('!')[1].replace(' ','').split(',')
                qm = [float(q) for q in qm]
                QME += qm
                break
        sf.close()
    nconf = len(QME)
    cc = {i+1:cc[i] for i in range(nconf)}

    #__________________________________________________________________________________
    ### storing the values of the torsional angles for each conformer in a dictionary
    #   to avoid calculating the dihedral angles so many times during optimization
    #----------------------------------------------------------------------------------
    dvals = {}
    for x in cc : 
        c = cc[x]
        dvals[x] = {}
        for d in dihedrals :
            [i,j,k,l] = dihedrals[d]['num']
            dvals[x][d] = dihedralvalue(c[i], c[j], c[k], c[l])

nconf = len(QME)
minQM = min(QME)
QME = [(q-minQM) for q in QME]
if not alternative :
    QME = [q*627.509 for q in QME]


#__________________________________________________________________________
### extracting the fitting paramaters from the confoiguration file (cnfg)
#--------------------------------------------------------------------------

### extracting the fitting parameters
sf = open('cfg')
frozper = []
froztor = []
repeat , phases , kmax , maxenergy , wmethod = 1 , 0 , 5 , 10 , 0
wts = [1 for q in QME]
nsteps = 5000

for line in sf :
    if not line.startswith('#') :
        line = line.split()
        n = len(line)
        if 'fp' in line and n > 5 : 
            line = line[1:]
            frozper += [[int(p) for p in line]]
        elif 'ft' in line and n == 8:
            line = line[1:]
            tor = [int(i) for i in line[:4]]
            tor += [float(line[4])]
            tor += [int(line[5])]
            tor += [float(line[6])]
            froztor += [tor]
        elif 'repeat' in line and n == 2:
            repeat = int(float(line[1]))
        elif 'phases' in line and n == 2 :
            phases = int(float(line[1]))
            if phases != 0 :
                phases = 1 
        elif 'kmax' in line and n == 2:
            kmax = float(line[1])
        elif 'maxenergy' in line and n == 2:
            maxenergy = int(float(line[1]))
        elif 'nsteps' in line and n == 2:
            nsteps = int(float(line[1]))
        elif 'weight' in line :
            wmethod = line[1]
            if wmethod == 'b' and n == 3:
                temp = float(line[2])
                wts = [exp(-1.0*q/(0.001987*temp)) for q in QME]
            elif wmethod == 'c' and n == 4: 
                cutoff = float(line[2])
                wt = float(line[3])
                wts = []
                for q in QME :
                    if q > maxenergy :
                        wts += [0.1]
                    elif q > cutoff :
                        wts += [1]
                    else :
                        wts += [wt]
sf.close()       


### using the entries of the configuration file to make the related modifications

# the fractional weights
wtsum = sum(wts)
wts = [w/wtsum for w in wts]
print ('\nThe following weights are in effect')
print (wts)
print()

### constraints on the values of the torsional parameters
if phases == 1 :    
    kmin = 0
    deltamax = 180.0
    deltamin = -180.0
    print ('The force constants will range between 0 and ' , kmax)
    print ('The phase shifte will range between -180 and 180 degrees ')
else:
    kmin = -kmax+0.0
    deltamax = 0
    deltamin = 0
    print ('The force constants will range between 0 and ' , kmax)
    print ('The phase shifte will either 0 or 180 degrees')


# assigning the frozen dihedrals and the frozen periodicities
for d in dihedrals :
    dihedrals[d]['f'] = []
    dihedrals[d]['p'] = []
    t1 = dihedrals[d]['num']
    
    pool = []    # to avoid the repeted periodicities
    for f in froztor :
        t2 = f[:4]
        t3 = [t1[3],t1[2],t1[1],t1[0]]
        if t1 in [t2,t3]:
            if not f[5] in pool :
                pool += [f[5]]
                dihedrals[d]['f'] += [f[4:]]
    dihedrals[d]['f'].sort()
    if not dihedrals[d]['f'] : 
        for f in frozper :        
            t2 = f[:4]
            t3 = [t1[3],t1[2],t1[1],t1[0]]
            if t1 in [t2,t3]:
                dihedrals[d]['p'] = f[4:]
        dihedrals[d]['p'].sort()


# equivqlencing the dihedrals 
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

print ('\nThe equivalencing sheme has detected the following torsional angles to be equivalent')
for d in equidihedrals :
    if len(equidihedrals[d]) > 1 : 
        print ([dihedrals[q]['num'] for q in equidihedrals[d]])

# eliminating the conflicting frozen torsionals and periodicitiess from equivalent dihedrals
for i in equidihedrals :
    for j in equidihedrals[i] :
        if j > i :
            if dihedrals[i]['f'] and dihedrals[j]['f'] :
                if dihedrals[i]['f'] != dihedrals[j]['f'] :
                    print (i , j)
                    ks = [k for k in equidihedrals[i] if not k in [i,j]]
                    k = min(ks)
                    equidihedrals[k] = ks
                    equidihedrals[i] = [i]
            elif dihedrals[i]['p'] and dihedrals[j]['p'] :
                if dihedrals[i]['p'] != dihedrals[j]['p'] :
                    print (i , j)
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

print ('\nTo consider the frozen torsionals and periodicities entries in the configuration file,')
print ('the following equivalent torsionals will be in effect')
for d in equidihedrals :
    if len(equidihedrals[d]) > 1 : 
        print ([dihedrals[q]['num'] for q in equidihedrals[d]])



'''        

bondpars = {
1:{'a':[1, 2],   'f': 333.2696,  'q':1.0924}, 
2:{'a':[1, 3],   'f':333.2696,   'q':1.0924}, 
3:{'a':[1, 4],   'f':333.2696,   'q':1.0924}, 
4:{'a':[1, 5],   'f':217.9119,   'q':1.5157}, 
5:{'a':[5, 6],   'f':688.4041,   'q':1.2327}, 
6:{'a':[5, 7],   'f':365.0012,   'q':1.3674}, 
7:{'a':[7, 8],   'f':464.3642,   'q':1.0101}, 
8:{'a':[7, 9],   'f':257.2852,   'q':1.4499}, 
9:{'a':[9, 10], 'f':333.2696,   'q':1.0924},
10:{'a':[9, 11], 'f':333.2696,   'q':1.0924},
11:{'a':[9, 12], 'f':333.2696,   'q':1.0924}}


ureypars = {
1:{'a':[2, 3],      'f':4.9969,      'q':1.7781},
2:{'a':[2, 4],      'f':4.9969,      'q':1.7781},
3:{'a':[2, 5],      'f':8.1751,      'q':2.134},
4:{'a':[3, 4],      'f':4.9969,      'q':1.7781},
5:{'a':[3, 5],      'f':8.1751,      'q':2.134},
6:{'a':[4, 5],      'f':8.1751,      'q':2.134},
7:{'a':[10, 11],   'f':4.9969,      'q':1.7781},
8:{'a':[10, 12],   'f':4.9969,      'q':1.7781},
9:{'a':[11, 12],  'f':4.9969,      'q':1.7781}}

anglepars = {
1:{'a':[2, 1, 3],       'f':11.7833,    'q':108.9635},
2:{'a':[2, 1, 4],       'f':11.7833,    'q':108.9635},
3:{'a':[2, 1, 5],       'f':18.6423,    'q':110.1937},
4:{'a':[3, 1, 4],       'f':11.7833 ,     'q':108.9635},
5:{'a':[3, 1, 5],       'f':18.6423,    'q':110.1937},
6:{'a':[4, 1, 5],       'f':18.6423,    'q':110.1937},
7:{'a':[1, 5, 6],       'f':35.6106,    'q':121.9856},
8:{'a':[1, 5, 7],       'f':44.9334,    'q':114.9637},
9:{'a':[6, 5, 7],       'f':35.2025,    'q':123.0315},
10:{'a':[5, 7, 8],       'f':17.2267,    'q':118.8443},
11:{'a':[5, 7, 9],       'f':41.5148,    'q':122.1051},
12:{'a':[8, 7, 9],       'f':16.7079,    'q':118.9779},
13:{'a':[7, 9, 10],     'f':17.4576,    'q':109.7152},
14:{'a':[7, 9, 11],     'f':17.4576,    'q':109.7152},
15:{'a':[7, 9, 12],     'f':17.4576,    'q':109.7152},
16:{'a':[10, 9, 11],   'f':11.7833,    'q':108.9635},
17:{'a':[10, 9, 12],   'f':11.7833,    'q':108.9635},
18:{'a':[11, 9, 12],   'f':11.7833,    'q':108.9635}}

improperpars = {
1:{'a':[5, 1, 7, 6],      'f':43.5061,   'q': 0.0000}}


charge = {
1:   -0.265, 
2:    0.090,
3:    0.090,
4:    0.090,
5:    0.511,
6:   -0.517,
7:   -0.471,
8:    0.307,
9:   -0.105,
10:  0.090, 
11:  0.090,
12:  0.090}
'''

       
LJ = {
1:   [ -0.078, 	2.05, 	-0.01, 	1.9      ],
2:   [ -0.024, 	1.34, 	-0.024, 	1.34    ],
3:   [ -0.024, 	1.34, 	-0.024, 	1.34    ],
4:   [ -0.024, 	1.34, 	-0.024, 	1.34    ],
5:   [ -0.11, 	2.0,  	-0.11, 	2.0      ],
6:   [ -0.12, 	1.7,  	-0.12, 	1.4      ],
7:   [ -0.2 ,	1.85, 	-0.2, 	1.55    ],
8:   [ -0.046, 	0.2245, 	-0.046, 	0.2245],
9:   [ -0.078, 	2.05, 	-0.01, 	1.9      ],
10: [ -0.024, 	1.34, 	-0.024, 	1.34    ],
11: [ -0.024, 	1.34, 	-0.024, 	1.34    ],
12: [ -0.024, 	1.34, 	-0.024, 	1.34    ]} 

if not alternative :
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
            rmin = (ri+rj)
            ene += eps * ((rmin/r)**12 -2*(rmin/r)**6)  
        MME += [ene]

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

print ('\nThe following periodicities will be in effect')
for d in dihedrals :
    if dihedrals[d]['p'] : 
        print (dihedrals[d]['num'] + dihedrals[d]['p'])

#////////////////////////////////////////////////////////////////////////////////////////////////


#__________________________________________________________________________________________
### The Annealing process

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
                    #print k , per , f , delta , cos((per * f - delta)*pi/180)  , te
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
p1 = 0.7                         # Probability of accepting worse solution at the start
p50 = 0.001                      # Probability of accepting worse solution at the end
t1 = -1.0/log(p1)                # Initial temperature
t50 = -1.0/log(p50)              # Final temperature
frac = (t50/t1)**(1.0/(n-1.0))   # Fractional reduction every cycle
na = na + 1.0
t = t1 + 0   # current temperature
dEav = 0.0
if nsteps : 
    n = int(nsteps/m)
print ('\nThe anneling parameters')
print ('    The starting temperature = {} K'.format(t1))
print ('    The final temperature = {} K'.format(t50))
print ('    The cooling scheme: temp(n+1) =  temp(n)*' , frac)
print()
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

for k in range(repeat) :
    print ('\n============================================')
    print ('Trial torsional fitting number {}'.format(k+1))


    rmsc = rmse(nconf,QME,MME,TEc,wts,maxenergy)
    rmsf = rmsc + 0.0

    print ('RMSE = ' , rmsf)

    #________________________
    # the iteration process
    #------------------------
    for i in range(n):
        for j in range(m):
            # Generate new trial points
            fitn = generatepoint(fitc,kmin,kmax,deltamin,deltamax)
            TEn = fitenergy(fitn) 
            rmsn = rmse(nconf,QME,MME,TEn,wts,maxenergy)
       
            #print xn , fn
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
                print ('RMSE = ' , rmsf)
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


    ### adding the fitted torsional energies to the MME energies
    MMf = [MME[i]+TEf[i] for i in range(nconf)]
    minMMf = min(MMf)
    MMf = [MMf[i]-minMMf for i in range(nconf)]

    print ('RMSE = ' , rmsf)

    print ('\nThe Quantum mechanical and the Molecular mechanical energies ')
    for i in range(nconf) :
        print ('%3s %15s %15s'% (i , QME[i] , MMf[i]))



print ('\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> .. Torsional parameterization ended successfully \n\n\n\n')
    
    


