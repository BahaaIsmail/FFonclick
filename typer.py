#!/usr/bin/python

#!/usr/bin/python


### this typer containes the tested and verified fixar.py
### find it in the directory "/home/bahaa/work_fix/fixar_tested_and_verified_2017-11-11/fixar.py"

### typer 29-8-2016

# extracting atom, numbers and cartesian coordinates from mol2 file
import sys
from math import sqrt , cos , acos , sin , atan2 , pi ,  ceil 

print ( '\n\n' )
print ( '__________________________________________________________________________________' )
print ( 'Lennard-Jones Parameterization ___________________________________________________' )


cfg = sys.argv[-1]
if cfg == '0' :
    cfg = 0


#system('python geom.py ' + mol2)

################################################################################
################################################################################
################################################################################
################################################################################
### ifuncs
################################################################################
################################################################################
################################################################################
################################################################################

def conic (i,ic) :
    # turns a list of atoms in numbers into a string of their labels
    # sometimes i need to coine a string to an array of atoms so tha
    # input: 2,[1,4,5,3]
    # output: 'CCNCO'
    conas = [atom[i]['a']] + sorted([atom[j]['a'] for j in ic ])
    return ''.join(conas)


def strcona (L): 
    # turns a list of atoms in numbers into a string of their labels
    return ''.join(sorted([atom[j]['a'] for j in L]))


def shape(i) : 
    return atom[i]['a']+str(atom[i]['ncon'])+str(int(atom[i]['nbs']))

def attr(i) :
    # returns some of the atributes of the atom i
    # inoput: the index of the atom
    # output: its symbol, inexes of thr connected atoms, a string of the symbols 
    #         of the connected atoms and the number of the connected atoms
    return atom[i]['a'] , atom[i]['cona'] , \
           atom[i]['con'] , atom[i]['ncon']

def intercept (L1, L2) :
    # returns the commone indices between two lists
    # input: [[1,5,4] , [2,4,5]]
    # ouput: [5,4]
    return [i for i in L1 if i in L2]


def diff (L1,L2): 
    # see: def intercept
    return [i for i in L1 if not i in L2]


def got (i,L,x) :
    # this functions searches for certain connectivities, a connectivity here involves  
    # bond order followed by atom symbol   
    # example; is the atom i has a double bond with a nitrogen atom
    # input: atom index, a list to search in and a list of expected connectivities 
    # output: the list of connected symbols preceded by the bondorder  
    gots = []   
    for j in range(atom[i]['ncon']): 
        s = atom[i]['cona'][j]
        if s in L and str(atom[i]['ord'][j]) + atom[s]['a'] in x : 
            gots += [s]
    return gots


def add(i,L) : 
    if not i in L : 
        L += [i]
    

def adds (L1,L2) : 
    L2 += [i for i in L1 if not i in L2]

   
def hetero(L) : 
    hets = []
    for i in L : 
        if not atom[i]['a'] in ['C','H'] :
            hets += [i]
    return hets
        



def onecyc(L,N): 
    # L is the list of atoms tested wether they are in one cycle with size n in the list n
    Ls = sorted(L)
    for c in cycle: 
        if cycle[c]['n'] in N : 
            if intercept(Ls,cycle[c]['num']) == Ls :
                return 1 



def cg(i,x):
    # assigns the atom type provided it is not already assigned
    if not atom[i]['cg']:
        atom[i]['cg'] = x
        if i in uniqs : 
            for j in uniqs[i] : 
                atom[j]['cg'] = x
                
def cgd(L,d): 
    # d is a dictionary with the atom symbols are the keys and the types are the definition 
    # calles cg() to assign the atom type for key index in a dictionary 
    for j in L :
        x = atom[j]['a']
        if x in d:
            cg(j,d[x]) 
            


def cgs (L , x) :
    for j in L : 
        cg(j,x)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

'''
### getting the mol2 and the srt file
srt = mol2.replace('.mol2','')
sf = open(srt)

typs = {}
penalty = {}
i = 1
for line in sf :
    if line.startswith('ATOM') :
        line = line.split()
        typs[i] = line[2]
        if not 'LP' in typs[i] : 
            #print ( line[5].strip('!') )
            penalty[i] = float(line[5].strip('!'))
        i += 1
sf.close()
'''

sf = open('tmp')
line = sf.readline().split()
atom = {i+1:{'a':line[i] , 'ord':[] , 'cona':[] , 'ncon':0, 'c':0 , 'cg':''} for i in range(len(line))}
nats = len(atom)

# modifying Cl, Br, Al, He , Ne and DUM to not be conflicted
# these names will be corected at the end of the typing steps
d = {'Cl':'L','Br':'B','Al':'A','He':'M','Ne':'E','DUM':'D','DUMM':'U','Se':'S', 
     'CL':'L','BR':'B','AL':'A','HE':'M','NE':'E'}
for i in atom :
    a = atom[i]['a']
    if a in d :
        atom[i]['a'] = d[a]
        
bondlist = []
bond = {}       
k = 0   
for line in sf :    
    k += 1
    bond[k] = {}
    [i,j,bo] = line.split() 
    bond[k]['num'] = sorted([int(i) , int(j)])
    if bo in ['ar','Ar','AR','co2','CO2','Co2','1.5'] :
        bo = 1.1
    elif bo == '2.5' :
        bo = 2.1
    else :
        bo = int(float(bo))
    bond[k]['ord'] = bo
    [i,j] = bond[k]['num']
    atom[i]['ord'] += [bo]
    atom[j]['ord'] += [bo]
    atom[i]['cona'] += [j]
    atom[j]['cona'] += [i]
    bon = sorted([i , j])
    bondlist += [bon]
    bond[k]['ord'] = bo
    

### ===================================
### more attributes for atom dictionary
for i in atom :
    if atom[i]['ord'] :
        atom[i]['o'] = max(atom[i]['ord'])   # the max order
    ncon = len(atom[i]['cona'])
    atom[i]['ncon'] = ncon
    nbs = sum(atom[i]['ord'])     # make sure
    
    if int(nbs) == nbs : 
        nbs = int(nbs)
    atom[i]['nbs'] = round(nbs,1) 
    atom[i]['con'] = conic(i , atom[i]['cona'])  
 


################################################################################
################################################################################
################################################################################
################################################################################

###### symmetries
### this is the symmetrizer part of structure
## cis-trans, axial-equatorial, resonance configurations are not in effect
vascular = {}
for i in atom :
    gener = sorted([atom[j]['con'],j] for j in atom[i]['cona']) 
    vascular[i] = [atom[i]['con']]+ [j[0] for j in gener]
    pool = [i]+[j[1] for j in gener] 
    while len(pool) < nats : 
        gener = [sorted([[atom[j]['con'],j] for j in atom[k[1]]['cona'] if not j in pool]) for k in gener]
        gener = [j for k in gener for j in k ]
        vascular[i] += [j[0] for j in gener]
        pool += [j[1] for j in gener]

uniqs = {}   # unique atoms as keys and the symmetric atoms as definitions 
pool = []
for i in atom:
    uniqs[i] = [j for j in atom if j != i and vascular[i] == vascular[j]]


uniq = [i for i in uniqs] # 









###############################################################################
################################################################################
################################################################################
################################################################################





### ===================================
### extracting frozen and equivalent 
# Lennard-Lones parameters from the 
# configuration file if exist
nosymm = 0
if cfg : 
    modif = 0
    LJ = {}
    eqlj = {}
    sf = open('cfg')
    for line in sf : 
        if not '#' in line : 
            if 'nosymmetry' in line : 
                nosymm = 1
            if 'flj' in line : 
                line = line.split()
                if len(line) == 6 : 
                    LJ[int(line[1])] = [float(j) for j in line[2:]]
            elif 'elj' in line : 
                line = line.split()
                if len(line) > 2 : 
                    line = sorted([int(j) for j in line[1:]])
                    eqlj[line[0]] = line[1:]
        if LJ or eqlj : 
            modif = 1
    sf.close()

if nosymm : 
    uniqs = {i:[] for i in atom}

### ===================================
### merging the equivalent atoms in the 
# configuration file with those detected 
# by the equivalencing scheme of Onclick
if cfg : 
    for i in uniqs : 
        for j in eqlj : 
            if intercept([j]+eqlj[j] , [i]+uniqs[i]) : 
                ks = [k for k in [j]+eqlj[j] if not k in [i]+uniqs[i]]
                uniqs[i] += ks 
    # joining the interceptions in uniqs after adding the eqlj
    s = [i for i in uniqs]
    s.reverse()
    for i in s : 
        for j in s : 
            if j < i :
                if intercept([j]+uniqs[j] , [i]+uniqs[i]) : 
                    uniqs[j] += [k for k in [i]+uniqs[i] if not k in uniqs[j] ]
                    uniqs[i] = [0]
    uniqs = {i:[j for j in uniqs[i] if j != 0 and j != i] for i in uniqs if uniqs[i] != [0]}

 
### =============================================
### freezing the atoms equivalet to frozen atoms
if cfg : 
    s = [i for i in LJ]
    for i in s :

        for j in uniqs : 
            if i in [j]+uniqs[j] : 
                for k in [j]+uniqs[j] : 
                    LJ[k] = LJ[i]



if cfg and modif: 
    print ()
    print ( 'According to the  the configuration file the program has modified the equivalences to be ' )
    for i in uniqs :
        if uniqs[i] : 
            eq = [i] + uniqs[i]
            print ( eq )
    print ()




###############################################################################
################################################################################
################################################################################
################################################################################


####============================================================================
#          _____________         ________________
#         /\            \       /\               \ 
#         \ \    ________\      \ \______    _____\
#          \ \   \_______/       \/_____/\   \____/
#           \ \   \________             \ \   \
#            \ \           \             \ \   \
#             \ \________   \             \ \   \
#              \/_______/\   \             \ \   \
#                 ______\_\   \             \ \   \
#                /\            \             \ \   \
#                \ \____________\             \ \___\
#                 \/____________/              \/___/
####============================================================================
# the angles, dihedrals, impropers, ureys, cycles

def addbonds (ics) :
    if ics :        
        for i in ics :
            nic = len(ics[i]['num'])
            ic = ics[i]['num']
            ics[i]['bon'] = []
            ics[i]['ord'] = []
            for j in range(nic-1) :
                bon = sorted([ic[j],ic[j+1]])
                for k in bond :
                    if bon == bond[k]['num'] :
                        ics[i]['bon'] += [k]
                        ics[i]['ord'] += [bond[k]['ord']] 
                    
                    
## creating anglelist 
anglelist = []
for i in atom :
    for j in atom[i]['cona'] :
        for k in atom[i]['cona'] :
            if k > j :
                anglelist += [[j , i , k]]

anglelist.sort()
nangs = len(anglelist)
angle = {i+1:{'num':anglelist[i]} for i  in range (nangs)}
addbonds(angle)


### creating dihedrallist
dihedrallist = []
for i in bondlist :
    for j in atom[i[0]]['cona'] :
        for k in atom[i[1]]['cona'] :
            if j != i[1] and k != i[0]:
                if k > j :
                    dihedrallist += [[j] + i + [k]]
                elif k < j :
                    dihedrallist += [[k , i[1] , i[0] , j]]
dihedrallist.sort()
ndihs = len(dihedrallist)
dihedral = {i+1:{'num':dihedrallist[i]} for i  in range (ndihs)}
addbonds(dihedral)

### getting ureys from angle
ureylist = []                # ureys = [[1 , 3] , [] , [] , ..]  in terms of atom numbers
for a in anglelist :
    us = [a[0] , a[2]]
    if not us in ureylist :
        ureylist += [us]
nurs = len(ureylist)



### creating improper from bond
improperlist = []
for i in atom :
    if atom[i]['a'] == 'C' and  atom[i]['ncon'] == 3 :
        k = atom[i]['cona']
        for j in range (3) :
            if atom[k[j]]['a'] == 'O' and \
               atom[i]['ord'][j] == 2 :
                improperlist += [[i] + k] 
nimps = len(improperlist)
improper = {i+1:{'num':improperlist[i]} for i  in range (nimps)}
addbonds(improper)





###### cycles
cyclelist = []
def intercept (ic1, ic2) :
    return [i for i in ic1 if i in ic2]

def diff (ic1, ic2) :
    return [i for i in ic1 if not i in ic2]

def identical (ic1 , ic2 ):
    return sorted(ic1) == sorted(ic2)     


def momcyc (ic):
    for cyc in cyclelist :        
        shared = intercept(ic,cyc)
        if identical (cyc , shared):
            return  1
        elif len(ic) > len(cyc) and len(shared) > 3 : 
            return 1

def addcyc (ic) :
    if not momcyc(ic) :
        return 1
    
               
# 3-mem
for a in anglelist : 
    if a[-1] in atom[a[0]]['cona'] and addcyc(a): 
        cyclelist += [a]

# 4-mem
for d in dihedrallist : 
    if d[-1] in atom[d[0]]['cona'] and addcyc(d): 
        cyclelist += [d]
     

# 5-mem cycles :
for d in dihedrallist :
    ic1 = atom[d[0]]['cona']
    ic2 = atom[d[-1]]['cona'] 
    fifth = intercept(ic1,ic2)
    fifth = diff(fifth , d)
    for j in fifth :
        if addcyc (d+[j]) :
            cyclelist += [d+[j]]


# 6-mem and 7-mem cycles
def cyc67 (ics1 , ics2):
    cycs = []
    for d in ics1:
        d1 = atom[d[0]]['cona']
        d2 = atom[d[-1]]['cona']
        for i in ics2 :
            if not intercept(d,i):
                b = [i[0] , i[-1]]
                j = intercept(d1,b)
                k =  intercept(d2,b)
                if j and k and j[0] != k[0] :
                    cyc = d + [k[0],j[0]]
                    if len(i) == 3 :
                        cyc = d + [k[0],i[1],j[0]]
                    cycs += [cyc]
    return cycs 
               
cycs67 = cyc67(dihedrallist,bondlist) + cyc67(dihedrallist,anglelist)
for cyc in cycs67 :
    if addcyc(cyc) :
        cyclelist += [cyc] 



ncycs = len(cyclelist)
## cycle list rearrangements
for c in range(ncycs) :
    ic = cyclelist[c]   
    m = ic.index(min(ic))
    icm = ic[m]
    ic = ic[m+1:]+ic[:m]
    if ic[0] > ic[-1] :
        ic.reverse()
    cyclelist[c] = [icm] + ic

      
cycle = {i+1:{'num':cyclelist[i]} for i  in range (ncycs)}
addbonds(cycle)
for i in cycle :   # to take the last atom into account
    a1 = cycle[i]['num'][0]
    a2 = cycle[i]['num'][-1]
    b = sorted ([a1,a2])
    for j in bond :
        if b == bond[j]['num'] :
            cycle[i]['bon'] += [j]
            cycle[i]['ord'] += [bond[j]['ord']]
            break

for i in cycle :
    cycle[i]['n'] = len(cycle[i]['num'])

cycle4 = [c for c in cycle if cycle[c]['n'] == 4]
cycle5 = [c for c in cycle if cycle[c]['n'] == 5]

## assigning the atom[i]['c']
for c in cycle :     
    n = cycle[c]['n']   
    for i in cycle[c]['num']: 
        if not atom[i]['c'] : 
            atom[i]['c'] = n



################################################################################
################################################################################
################################################################################
################################################################################
# in case at least one aromatic bond exists the program will extract the aromatic 
# rings from the input file  

fixaromatics = 0
for b in bond : 
    if bond[b]['ord'] == 1.1 : 
        fixaromatics = 1 
        break 
    
aromatics = {c:0 for c in cycle}
if fixaromatics: 
    for c in cycle : 
        ords = round(sum(cycle[c]['ord']),1)
        if ords in [5.5 , 6.6 , 7.7]: 
            aromatics[c] = 1

#### fixing the aromatic bonds if existing
######## fixing the ar with only one option

if fixaromatics :
    Cb = []
    Nb = []
    for c in aromatics :         
        for s in aromatics : 
            if s > c : 
                cyc = cycle[c]['num']
                cys = cycle[s]['num']                
                for i in intercept(cyc,cys):
                    br = 1
                    for j in atom[i]['cona'] : 
                        if not (j in cyc or j in cys):
                            br = 0
                            break
                    if br : 
                        if atom[i]['a'] == 'C' : 
                            add(i,Cb)
                        elif atom[i]['a'] == 'N' : 
                            add(i,Nb) 

    def casea(i) : 
        if i in Cb : 
            return 'Cb'
        if i in Nb : 
            return 'Nb'
        return atom[i]['a'] + str(int(atom[i]['nbs']))            
        
    def caseb(b):
        return sorted([casea(i) for i in bond[b]['num']]) 

    def npi(c): 
        d = []
        for b in cycle[c]['bon']: 
            if bond[b]['ord'] == 2 :
                adds(bond[b]['num'],d)
        for i in cycle[c]['num']: 
            if atom[i]['a'] == 'C' and (atom[i]['nbs'] < 4 or i in Cb): 
                add(i,d)            
            if atom[i]['a'] == 'N' and atom[i]['nbs'] < 3 : 
                add(i,d)
        h1 = [i for i in cycle[c]['num'] if i in Nb and not i in d]
        h2 = [i for i in cycle[c]['num'] if atom[i]['a'] in ['O','S'] or (casea(i) == 'N3' and not i in d)]
        return len(d) + len(h1) + 2*len(h2) 

    npis = {0:[], 1:[], 2:[], 3:[], 4:[], 5:[], 6:[], 7:[]}    
    for c in aromatics :        
        n = npi(c) - 6 
        if n in [1,3,5,7] : 
            if cycle[c]['n'] == 5 : 
                n = n + 1
            elif cycle[c]['n'] == 7: 
                n = n - 1
        if n <= 0 : 
            npis[0] += [c]        
        else : 
            npis[n] += [c]
   
    
    def update(b,o): 
        bond[b]['ord'] = o 
        #print ( bond[b]['num'] , o )
        [i,j] = bond[b]['num']
        
        k = atom[i]['cona'].index(j)
        atom[i]['ord'][k] = o
        atom[i]['nbs'] = sum(atom[i]['ord'])
        atom[i]['o'] = max(atom[i]['ord'])
        
        l = atom[j]['cona'].index(i)
        atom[j]['ord'][l] = o
        atom[j]['nbs'] = sum(atom[j]['ord'])
        atom[j]['o'] = max(atom[j]['ord'])
    
        for c in aromatics : 
            if b in cycle[c]['bon']: 
                k = cycle[c]['bon'].index(b)
                cycle[c]['ord'][k] = o
    def N1(i): 
        for b in bond : 
            if i in bond[b]['num'] and bond[b]['ord'] == 1.1 :
                update(b,1)
                
    def Bij(i,j,o): 
        for b in bond : 
            if bond[b]['ord'] == 1.1 : 
                if sorted([i,j]) == bond[b]['num']: 
                    update(b,o)
                    break 

    def A(): 
        fixedout = 0
        for step in range(99):    
            fixedin = 0
            for i in atom : 
                if 1.1 in atom[i]['ord'] : 
                    L = ['O2.1','O2.2','S2.1','S2.2','P5.1','P5.2','S6.1','S6.2']
                    if atom[i]['nbs'] - atom[i]['ncon'] >= 1 or  atom[i]['a'] + str(round(atom[i]['nbs'],1)) in L : 
                        N1(i)                          
                        fixedout = 1
                        fixedin = 1   
                        break
                    elif atom[i]['a'] + str(round(atom[i]['nbs'],1)) in ['C3.1','N2.1'] : 
                        for b in bond : 
                            if bond[b]['ord'] == 1.1 and i in bond[b]['num']:                                 
                                update(b,2)
                                fixedout = 1
                                fixedin = 1
                                break
            if not fixedin : 
                return fixedout
            
    def SO():
        fixed = 0
        for a in ['O','S']: 
            for i in atom : 
                if atom[i]['a'] == a and atom[i]['nbs'] == 1.1:                     
                    j = atom[i]['cona'][0]
                    ates = len([k for k in atom[j]['cona'] if atom[k]['a'] in ['O','S'] and atom[k]['nbs'] == 1.1])
                    if ates > 1 :
                        Bij(i,j,2)
                        N1(j)
                        fixed = 1
        
        for c in aromatics : 
            if c in npis[6] : 
                cyc = cycle[c]['num']
                for i in cycle : 
                    if atom[i]['ncon'] == 3 :
                        j = diff(atom[i]['cona'],cyc)[0] 
                        if atom[j]['a'] in ['O','S'] and atom[j]['ncon'] == 1 and atom[j]['nbs'] == 1.1 : 
                            Bij(i,j,1)
                            fixed = 1
            return fixed
            
    for itera in range(99):        
        fixed = 0        
        if not SO() :
            if not A():
                break


fixaromatics = 0
for b in bond : 
    if bond[b]['ord'] == 1.1 : 
        fixaromatics = 1 
        break 
    

if fixaromatics: 
    for c in npis[0] :        
        if 1.1 in cycle[c]['ord'] : 
            for i in cycle[c]['num'] :
                if 1.1 in atom[i]['ord'] and casea(i) == 'N3': 
                    N1(i)
                    A()
            for b in cycle[c]['bon'] : 
                if caseb(b) == ['Nb','Nb']: 
                    update(b,1)
                    A()
 
                            
    stat = {c:{'npi':0 , 'nC':0} for c in aromatics if 1.1 in cycle[c]['ord']}
    for n in npis : 
        for c in npis[n] : 
            if 1.1 in cycle[c]['ord']:
                stat[c]['npi'] = n 
    for c in aromatics :         
        if 1.1 in cycle[c]['ord'] : 
            n = 0
            for i in cycle[c]['num'] : 
                if 1.1 in atom[i]['ord'] : 
                    if casea(i) in ['C3','N2','Cb']: 
                        n += 1
            stat[c]['nC'] = n

    for np in npis : 
        for c in npis[np]:         
            if 1.1 in cycle[c]['ord']: 
                npi = stat[c]['npi']
                nC  = stat[c]['nC']
                nd = min(3,int((npi + nC)/2))
                if cycle[c]['n'] == 6 and nd == 3 :  
                    for b in cycle[c]['bon'] :
                        if bond[b]['ord'] == 1.1 : 
                            update(b,2)
                            A()                            
                else : 
                    def doub(b) : 
                        cas = caseb(b)
                        nn = cas.count('N3') 
                        nc = cas.count('C3') + cas.count('N2') + cas.count('Cb')
                        nNb = cas.count('Nb')
                        return nn , nc , nNb
                        
                    if nd == 1 : 
                        fixed = 0
                        for l in [['C3','C3'],['C3','N2'],['C3','N3'],['N2','N2'],['N2','N3'],
                                  ['C3','Cb'], ['N2','Cb'], ['N3','Cb'],
                                  ['C3','Nb'], ['N2','Nb'], ['N3','Nb']]:
                            for b in cycle[c]['bon']: 
                                if bond[b]['ord'] == 1.1: 
                                    if caseb(b) == l  :                        
                                        update(b,2)
                                        A()
                                        fixed = 1
                                        break
                            if fixed : 
                                break
                    elif nd == 2 :    
                        #print ( cycle[c]['num'] , npi , nC , nd , 'hhhhhhhhhhhhhh'
                        nNbo = 9999999
                        b1f = 0                    
                        for b1 in cycle[c]['bon']: 
                            if bond[b1]['ord'] == 1.1 : 
                                nn1 , nc1 , nNb1 = doub(b1)
                                for b2 in cycle[c]['bon']: 
                                    if bond[b2]['ord'] == 1.1 : 
                                        if b2 > b1 and not intercept(bond[b1]['num'],bond[b2]['num']):
                                            nn2 , nc2 , nNb2 = doub(b2)
                                            nn = nn1+nn2
                                            nc = nc1+nc2
                                            nNb = nNb1+nNb2
                                            #print ( bond[b1]['num'] , bond[b2]['num'] , nn , nc , nNb 
                                            if nc == nC and nn == npi and nNb < nNbo :                                                
                                                b1f = b1
                                                b2f = b2
                                                nNbo = nNb + 0
                        if b1f : 
                            update(b1f,2)
                            update(b2f,2)
                            A()
                    elif cycle[c]['n'] == 7:                     
                        nNbo = 9999999
                        b1f = 0
                        for b1 in cycle[c]['bon']: 
                            if bond[b1]['ord'] == 1.1 : 
                                nn1 , nc1 , nNb1 = doub(b1)
                                for b2 in cycle[c]['bon']: 
                                    if bond[b2]['ord'] == 1.1 : 
                                        if not intercept(bond[b1]['num'],bond[b2]['num']):
                                            nn2 , nc2 , nNb2 = doub(b2)
                                            for b3 in cycle[c]['bon']: 
                                                if bond[b3]['ord'] == 1.1 : 
                                                    if not intercept(bond[b1]['num'],bond[b2]['num']) and not intercept(bond[b2]['num'],bond[b3]['num']):
                                                        nn3 , nc3 , nNb3 = doub(b3)                                
                                                        nn = nn1+nn2+nn3
                                                        nc = nc1+nc2+nc3
                                                        nNb = nNb1+nNb2+nNb3
                                                        if nc == nC and nn == npi and nNb < nNbo :
                                                            b1f = b1
                                                            b2f = b2
                                                            b3f = b3
                                                            nNbo = nNb + 0
                        if b1f : 
                            update(b1f,2)
                            update(b2f,2)
                            update(b3f,2)
                            A()
    for b in bond : 
        if bond[b]['ord'] == 1.1 : 
            if caseb(b) in [['Cb','Cb'], ['Cb','CN']] : 
                update(b,2)
                A()
            elif caseb(b) == ['Nb','Nb']: 
                update(b,1)
                A()
    for b in bond : 
        if bond[b]['ord'] == 1.1 : 
            update(b,1)
            A()
        

    
################################################################################
################################################################################
################################################################################
################################################################################
### aromatic and planar rings 

# rings (not containing ncon = 4 or atoms other than C,N,O,S,P  )

ring = {c:{'r':0 , 'h':[], 'd':[], 'x':[]} for c in cycle if cycle[c]['n'] in [5,6,7]}
for c in ring :
    itisring = 1
    cyc = cycle[c]['num']
    for i in cyc :        
        a , cona , con , ncon = attr (i)
        if not a in ['C','O','S','N','N+','P'] : 
            itisring = 0
        elif ncon == 4 : 
            itisring = 0
        elif atom[i]['ord'].count(2) + atom[i]['ord'].count(3) > 1 : 
            itisring = 0
            break
    if itisring :
        ring[c]['r'] = 1
        for i in cyc :
            a , cona , con , ncon = attr(i)
            ancon = shape(i)
            if ancon in ['O22','S22','N33','N+33'] :
                ring[c]['h'] += [i]
            elif ncon < atom[i]['nbs']: 
                if got(i , intercept(cyc,cona) , ['3C','2C','2N','2N+','2P']) :
                    if not i in ring[c]['d'] :
                        ring[c]['d'] += [i]                        
                else :
                    if not i in ring[c]['x'] :
                        ring[c]['x'] += [i]

## removing nonrings from ring 
ring  = {c:ring[c] for c in ring if ring[c]['r'] and len(ring[c]['h']) < cycle[c]['n']}
## extracting planars from ring before excluding nonaromatics
## planar rings will be further reduced after knowing furans 
### planars
planars = []
for c in ring :
    if cycle[c]['n'] == 5 and len(ring[c]['h']) < 2 :
        add (c,planars)


for c in planars : 
    print ( cycle[c]['num'] )

################################################################################
################################################################################
################################################################################
################################################################################
### in case no aromatic bonds exixting the program will predict the aromatic rings with the following scheme 
if not fixaromatics : 
    # removing the exocyclic double bonds in nonrings
    # removing nonaromatic rings 
    ringlist = [c for c in ring]
    for step in range(99) : 
        updated = 0
        for c in ringlist :
            for i in ring[c]['x'] : 
                j = diff(atom[i]['cona'],cycle[c]['num'])[0]
                if not [s for s in ring if j in ring[s]['d']] : # cycles in wich the exodouble cond is a double bond
                    ring[c]['x'].remove(i)
                    updated = 1
            d = len(ring[c]['d'])
            h = len(ring[c]['h'])
            x = len(ring[c]['x'])
            if d + 2*h + x < 6 or (cycle[c]['n'] < 6 and d + h > 6): 
                ring[c]['r'] = 0 
                updated = 1
        ring  = {c:ring[c] for c in ring if ring[c]['r']}
        ringlist = [c for c in ring]
        if not updated : 
            break
    
    
    # preparing the rings for detemining the aromatics
    # by deleting the strained rings (i.e. which share anothe rring wit more than 2 atoms)
    def strained(c) : 
        for s in cycle: 
            if s != c and len(intercept(cycle[c]['bon'],cycle[s]['bon'])) > 1 : 
                return 1 
    ringlist = [c for c in ring]
    for c in ringlist : 
        if strained(c) :         
            ring[c]['r'] = 0
    ringlist = [c for c in ring]
                           
    ################################################################################
    # collecting the aromatic rings for the current step 
    
    # permanently aromatic rings  
    aromatics = {c:0 for c in ring if ring[c]['r']}
    for c in ring : 
        if not ring[c]['x']: 
            d = len(ring[c]['d'])
            h = len(ring[c]['h'])
            if d + 2*h == 6 : 
                aromatics[c] = 1
    
    #print ( 'permanent aromatics = ' ,  aromatics )
    
    def dhx(c): 
        cyc = cycle[c]['num']
        x = ring[c]['x']
        d = ring[c]['d']
        h = ring[c]['h']        
        a = []  # hetero atom in the bridge
        r = []  # exocyclic double bond in an aromatic ring
        y = []  # exocyclic double bond and inring double bond in the same bridge
        for s in ring : 
            if s != c and aromatics[s] :            
                cys = cycle[s]['num']
                b = intercept(cyc,cys)
                for i in b : 
                    if i in h :
                        add(i,a)                            
                    elif i in x:
                        d6 = cycle[s]['n'] == 6 and len(ring[s]['d']) == 6  
                        if d6 : 
                            add(i,d)
                        else : 
                            add(i,r)
        return d,h,a,x,r 
    
    
    def lendhx(d,h,a,x,r) : 
        d = len(d)
        a = len(a)
        h = len(h)-a
        r = len(r)
        x = len(x)-r
        return d,h,a,x,r 
       
       
    def npi(d,h,a,r) :
        np = d + 2*h + a  
        
        #print ( c , d , 2*h , a , r , np
        #if np + a > 6 : 
        #    return 0
        if np <= 6 and np + a + r >= 6 :
            return 1
        elif cycle[c]['n'] == 7 and d == 6 and (h == 1 or a == 1) : 
            return 1
        return 0
    
        
    cand = [c for c in aromatics if not aromatics[c]]
    #print ( 'cand = ' , cand 
    for itr in range(99) : 
        #print ( itr , '_______________________'
        updated = 0    
        for c in cand :
            d,h,a,x,r = dhx(c)
            d,h,a,x,r = lendhx(d,h,a,x,r)
            getaromatic = npi(d,h,a,r) 
            if getaromatic != aromatics[c] : 
                aromatics[c] = getaromatic
                updated = 1
        if not updated : 
            break


################################################################################
################################################################################
################################################################################
################################################################################


'''
            _________    _____    _____
           /        /   /    /   /    /
          /      __    /____/   /____/
         /        /   /  \     /
        /________/   /    \   / 
   
'''
# functional groups


# adding "+" to positive nitrogenes and '-' to the negative oxygens and sulfurs
for i in atom : 
    if atom[i]['a'] =='N' and atom[i]['nbs'] == 4 : 
        atom[i]['a'] = 'N+'
    elif atom[i]['a'] in  ['O','S'] and atom[i]['nbs'] == 1 :
        atom[i]['a'] += '-'


for x in ['N+','O-','S-'] :
    for i in uniqs : 
        if atom[i]['a'] == x or [j for j in uniqs[i] if atom[j]['a'] == x] : 
            atom[i]['a'] = x
            for j in uniqs[i] : 
                atom[j]['a'] == x


# correcting 'con' definitions
for i in atom :
    atom[i]['con'] = conic(i , atom[i]['cona']) 


# the resonating nitrogens in gaunidine, gaunidinium and amidinium
guans = []
for i in atom : 
    if shape(i) == 'C34' : 
        cona = atom[i]['cona']   
        k = got(i,cona,['2N+']) 
        if k  :
            k = k[0]
            for j in cona : 
                if shape(j) == 'N33' : 
                    atom[j]['a'] = 'N+'
                    adds([i,j,k],guans)


# correcting 'con' definitions after guans
for i in atom :
    atom[i]['con'] = conic(i , atom[i]['cona']) 

# carbonyls and amides
carbonyl = []
amide = []
ester = []
for i in atom : 
    if shape(i) == 'C34' : 
        cona = atom[i]['cona']     
        if got(i,cona,['2O','2S','2O-','2S-']) : 
            add(i,carbonyl)
            for j in cona : 
                s = shape(j)
                if s == 'N33' : 
                    add(j,amide)
                elif s == 'O22' and not 'H' in atom[j]['con']: 
                    add(j,ester)


# conjugated systems of double bonds
conj = []
conjs = []
for d in dihedral :
    if dihedral[d]['ord'] == [2,1,2] :
        conjs += [dihedral[d]['num']]
        for i in dihedral[d]['num'] :
           if atom[i]['a'] == 'C' : 
               add(i,conj)
                  

# hetero atoms
heteros = []
for i in atom : 
    if shape(i) in ['O22','O-22','O-11','S22','S-22','S-11','N33'] : 
        heteros += [i]
    


# N+33 which can't share the +ve in amidinium in order to keep aromaticity
forbidden = [] 
for c in aromatics :
    d,h,a,x,r = dhx(c)
    d,h,a,x,r = lendhx(d,h,a,x,r)
    if d + 2*h + 2*a + r == 6 :         
        adds([i for i in ring[c]['h']],forbidden) 
    for i in ring[c]['d'] : 
        if shape(i) == 'N+34' and i in guans :            
            j = got(i,atom[i]['cona'],['2C'])[0]
            k = diff(atom[j]['cona'],cycle[c]['num'])[0]
            if shape(k) == 'N+33':                     
                if [s for s in aromatics if s != c and i in cycle[s]['num']] : 
                    a += 1
                else : 
                    h += 1
                if [s for s in aromatics if j in ring[s]['x']] :
                    r += 1                
                if not npi(d-1,h,a,r) : 
                    add(k,forbidden)


aromatics = [c for c in aromatics if aromatics[c]]    
#print ( 'aromatics : ' , aromatics

# classification of cyclic atoms
def collectar(n,L) : 
    R = [c for c in L if cycle[c]['n'] == n]
    A = []
    for c in R :
        for i in cycle[c]['num'] :
            add(i,A)
    return R , A
    
benzenes, benza = collectar(6,aromatics)
aromatic7 , aroma7 = collectar(7,aromatics)
furans, furas = collectar(5,aromatics)
fura = [i for i in furas if atom[i]['a'] in ['O','S']]
cycle5 , cyca5 = collectar(5,cycle)
planars = [c for c in planars if not c in aromatics]
planars , planas = collectar(5,planars)
plana = [i for i in planas if shape(i) == 'N33']
allaroma = furas + benza + aroma7

                   
def addmap(m) : 
    tf = open('typlist','a')
    tf.write(m+'\n')
    tf.close()
        
    
    
################################################################################
################################################################################
################################################################################
################################################################################

###### symmetries
### this is the symmetrizer part of structure
## cis-trans, axial-equatorial, resonance configurations are not in effect
vascular = {}
for i in atom :
    gener = sorted([atom[j]['con']+str(atom[j]['c']),j] for j in atom[i]['cona']) 
    vascular[i] = [atom[i]['con']+str(atom[i]['c'])]+ [j[0] for j in gener]
    pool = [i]+[j[1] for j in gener] 
    while len(pool) < nats : 
        gener = [sorted([[atom[j]['con']+str(atom[j]['c']),j] for j in atom[k[1]]['cona'] if not j in pool]) for k in gener]
        gener = [j for k in gener for j in k ]
        vascular[i] += [j[0] for j in gener]
        pool += [j[1] for j in gener]

uniqs = {}   # unique atoms as keys and the symmetric atoms as definitions 
pool = []
for i in atom:
    uniqs[i] = [j for j in atom if j != i and vascular[i] == vascular[j]]
uniq = [i for i in uniqs] #  



uniqs_short = {}     # to remove the repetitions
pool = []
for i in atom:
    if not i in pool :
        temp = [j for j in atom if j > i and vascular[i] == vascular[j]]
        uniqs_short[i] = temp
        pool += temp

 
if nosymm : 
    uniqs_short = {i:[] for i in atom}

print ()
print ( 'The program equivalencing scheme has detected the following equivalent atoms ' )
for i in uniqs_short :
    if uniqs_short[i] : 
        eq = [i] + uniqs_short[i]
        print ( eq )
print ()


























'''_______________                   _________
   \______________\   |\     |\     |\ _______\
        | |           | |    | |    | |     | |
        | |           | |____| |    | |_____| |
        | |            \|_____\|    | |______\|
        | |                  | |    | |      
        | |           |\     | |    | |      
        | |           | |____| |    | |      
         \|            \|_____\|     \|      
                    
         
'''
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

#### triples groups
for i in atom :
    a , cona , con , ncon = attr (i)
    if atom[i]['a'] == 'N' and atom[i]['ncon'] == 1:
        cg(i,'NG1T1')
    elif atom[i]['a'] == 'C' and atom[i]['ncon'] == 2:
        if got(i,atom[i]['cona'],['3N']) :
            cg(i,'CG1N1')
        elif got(i,atom[i]['cona'],['3C']): 
            if 'H' in atom[i]['con']: 
                cg(i,'CG1T2')  
            else : 
                cg(i,'CG1T1')
        elif len (got(i,atom[i]['cona'],['2O','2N','2N+','2S','2P'])) == 2 :
            cg(i,'CG2O7')
            cgd(atom[i]['cona'],{'O':'OG2D5','N':'NG2D1','S':'SG2D1'})
        cgd(atom[i]['cona'],{'H':'HGPAM1'})

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
### benzenes 

# biphenyl 
def bip(R):
    def excluded(i,j) : 
        for c in benzenes + aromatic7+furans: 
            cyc = cycle[c]['num']
            if i in cyc and j in cyc : 
                return 1        
    L = []
    for c in R : 
        cyc = cycle[c]['num']        
        for i in cyc :
            if atom[i]['ncon'] == 3  and atom[i]['c'] > 4:  
                j = diff(atom[i]['cona'],cyc)[0]
                if got(i,[j],['1C','1N','1N+']) :                     
                    for s in R : 
                        if s != c : 
                            cys = cycle[s]['num']
                            if j in cys and (not j in cyc) and (not i in cys) : 
                                if not excluded(i,j): 
                                    add([i,j],L)                       
                                    add([j,i],L)
    for c in benzenes :
        cyc = cycle[c]['num']
    return L
biphenyl = bip(benzenes)
biphenyla = []
for b in biphenyl : 
    adds(b,biphenyla)

for b in biphenyla : 
    cgd(biphenyla,{'C':'CG2R67' , 'N':'NG2R67' , 'N+':'NG2R67'})

for i in benza : 
    if atom[i]['c'] == 5: 
        cgd([i],{'C':'CG2RC0' , 'N+':'NG2RC0'})
        if atom[i]['a'] == 'N' and  i in plana or i in furas : 
            cg(i,'NG2RC0')
            

# carbonyl and thiocarbonyl in benzene 
CG2R62 = []
H2 = []
for c in benzenes : 
    cyc = cycle[c]['num']
    H62 = 0
    for i in cyc:
        j = got(i,atom[i]['cona'],['2O','2S','2O-','2S-'])
        if j :
            cgd([i],{'C':'CG2R63'})
            H62 = 1
            if 'O' in atom[j[0]]['a'] : 
                adds(cyc,CG2R62)
                cg(j[0],'OG2D4')
    if H62 :
        for i in cyc :
            if atom[i]['a'] =='C' and atom[i]['ncon'] == 3 : 
                H2 += [k for k in atom[i]['cona'] if atom[k]['a'] == 'H']

# fluorinated Carbons
fluorinated = []
for i in benza : 
    if 'F' in atom[i]['con']: 
        fluorinated += [i]
        if atom[i]['a'] == 'C' and atom[i]['con'].count('N') < 2 : 
            cg(i,'CG2R66')


# protonated and carbonyl benzenes
protonated = []
for c in benzenes : 
    cyc = cycle[c]['num']    
    H3 = 0
    for i in cyc : 
        if atom[i]['a'] == 'N+' and not i in forbidden: 
            adds(cyc,protonated)
            H3 = 1
            break
        elif atom[i]['a'] == 'C' and [j for j in atom[i]['cona'] if shape(j) == 'N+34' and j in benza and onecyc([i,j],[6,7])]: 
            add(i,protonated)
            # should it be in one cycle?
    if H3 : 
        for i in cyc :
            if atom[i]['a'] == 'C'  :                
                cgd(atom[i]['cona'],{'H':'HGR63'})




# amidine Carbons
for i in benza :     
    if atom[i]['a'] == 'C' and atom[i]['con'].count('N') > 1 : 
        if not i in protonated :  
            cg(i,'CG2R64')

       
cgd(CG2R62+protonated,{'C':'CG2R62'})

# Nitrogens
NG2R62 = []
for i in benza : 
    if not atom[i]['cg']: 
        if 'N' in atom[i]['a']: 
            if atom[i]['ncon'] == 3 : 
                cg(i,'NG2R61')
            else : 
                for j in atom[i]['cona'] : 
                    if j in benza and 'N' in atom[j]['a'] : 
                        add(i,NG2R62)
                        break
                    else : 
                        for k in atom[j]['cona']: 
                            if k != i and k in benza and 'N' in atom[k]['a'] : 
                                add(i,NG2R62)
                                break
                        if i in NG2R62 : 
                            break

cgs(NG2R62,'NG2R62')
for i in benza: 
    if 'N' in atom[i]['a'] : 
        cg(i,'NG2R60')

# the remaining atoms in benza (Carbons and Oxygens)
cgd(benza,{'C':'CG2R61'})


# H on benzenes
for i in H2 :
    if not atom[i]['cg'] : 
        cg(i,'HGR62')
for i in benza : 
    if atom[i]['a'] == 'C' :
        cona = atom[i]['cona']
        hs = [j for j in cona if atom[j]['a'] in 'N+OSPFBLI' or 'F' in atom[j]['con']] 
        if hs :
            cgd(cona,{'H':'HGR62'})
        else : 
            cgd(cona,{'H':'HGR61'})
			
'''
for i in benza : 
    a = atom[i]['a'] 
    if a in 'NOSP' : 
        for j in atom[i]['cona'] :
            if atom[j]['a'] == 'C' : 
                cgd(atom[j]['cona'],{'H':'HGR62'})
    elif a == 'C' : 
        con = atom[i]['con']
        if 'F' in con or 'B' in con or 'L' in con or 'I' in con : 
            for j in atom[i]['cona'] :
                if atom[j]['a'] == 'C' : 
                    cgd(atom[j]['cona'],{'H':'HGR62'})
      
for i in benza : 
    if atom[i]['a'] == 'C' : 
        cgd(atom[i]['cona'],{'H':'HGR61'})
'''

# Oxygens in 6-mem 
for i in atom : 
    if not atom[i]['cg'] and atom[i]['a'] == 'O' and atom[i]['c'] == 6 : 
        for j in atom[i]['cona'] :
            if atom[j]['nbs'] > atom[j]['ncon']:
                cg(i,'OG3R60')
                break
        if not atom[i]['cg']:
            cg(i,'OG3C61')




################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
### cycle5 and aromatic7


### bipyrrol
bipyrrols = bip(cycle5)
bipyrrol = []
for b in bipyrrols: 
    [i,j] = b
    if not atom[i]['cg'] and not atom[j]['cg'] and not onecyc(b,[5]): 
        bipyrrol += [b + [c for c in cycle5 if j in cycle[c]['num']]]

for b in bipyrrol : 
    [i,j,c] = b
    if atom[i]['ncon'] == 3 :
        a = atom[i]['a']
        if a == 'C': 
            cg(i,'CG2R57')
        elif 'N' in a :  
            a = atom[j]['a']
            if a == 'C' and atom[j]['ncon'] == 3 and not j in guans : 
                cg(i,'NG2R57')
            elif a == 'N' : 
                if j in amide or j in furas: 
                    cg(i,'NG2R57')

         
### aromatic 7-mem  
for i in aroma7 : 
    if atom[i]['a'] == 'C' : 
        if atom[i]['c'] == 5 : 
            cg(i,'CG2RC7')
        else : 
            cg(i,'CG2R71')
        cgd(atom[i]['cona'],{'H':'HGR71'})			
    

# O , N and S of fura 
cgd(furas,{'O':'OG2R50','S':'SG2R50'})


## N+34 and N23 with three connectivities in 5-mem
for c in cycle5 : 
    cyc = cycle[c]['num']
    for i in cyc : 
        if not atom[i]['cg']:
            if 'N' in atom[i]['a'] : 
                
                s = shape(i)
                if s in ['N+34','N+33'] : 
                    cg(i,'NG2R52')
                elif s == 'N23': 
                    cg(i,'NG2R50')

### N33 in plana
cgd(plana+furas,{'N':'NG2R51'})



# oxygen in cycle5
for i in cyca5 : 
	if atom[i]['a'] == 'O' : 
		cg(i,'OG3C51')  
  

### C=N (involving those adjacent to another heteroatom)
CG2R52 = []   # to make sure there is no heteroatom in another cycle
for c in cycle5: 
    cyc = cycle[c]['num']
    for i in cyc :
        if atom[i]['a'] == 'C' and atom[i]['ncon'] == 3 and not atom[i]['cg'] :                        
            cona = atom[i]['cona']
            con = atom[i]['con']
            cycona = intercept(cona,cyc)             
            cyconah = got(i,cycona,['1O','1S','1N','1N+'])
            if got(i , cona , ['2N','2N+']):                    
                if con.count('N') > 1:                     
                    cg(i,'CG2R53')                    
                elif cyconah :                    
                    cg(i,'CG2R53')  
                elif got(i , cycona , ['2N','2N+']): 
                    CG2R52 += [i]                     
            elif got(i,cona,['2O','2S','2O-','2S-']) and cyconah :                
                cg(i,'CG2R53')             

cgs(CG2R52,'CG2R52')

## the default C=C in 5mem
for c in cycle5: 
    cyc = cycle[c]['num'] 
    for i in cyc : 
        if atom[i]['a'] == 'C' and atom[i]['ncon'] == 3 and not atom[i]['cg'] :              
            cycona = intercept(atom[i]['cona'],cyc)
            if got(i,cycona,['2C','2N','2N+','2P']) :
                cg(i,'CG2R51')
                

# exocyclic double bonds
CG251O = []
CG25C1 = []
CG2D1 = []
for c in cycle5 : 
    cyc = cycle[c]['num']
    for i in cyc : 
        cona = atom[i]['cona']
        if atom[i]['a'] == 'C' and not atom[i]['cg'] and atom[i]['ncon'] == 3: 
            j = diff(cona,cyc)                      
            j = got(i,j,['2C','2N','2N+','2O','2S','2O-','2S-','2P'])
            if j :                 
                if not (onecyc([i]+j,[3,4,5,6,7]) and atom[j[0]]['c'] != 5):
                    cycona = diff(cona,j)
                    if [k for k in cycona if atom[k]['a'] in 'O-N+S-'] :  
                        add(i,CG251O)                
                    elif i in conj : 
                        add(i,CG25C1)
                    elif atom[j[0]]['a'] == 'C': 
                        add(i,CG2D1)
cgs(CG251O,'CG251O')
cgs(CG25C1,'CG25C1')
cgs(CG2D1,'CG2D1')

                    
for i in cyca5 : 
    if atom[i]['a'] == 'C' and atom[i]['ncon'] == 3 : 
        con = atom[i]['con']
        n = len(got(i,atom[i]['cona'],['2N','2N+']))
        if n and con.count('N') + con.count('O') + con.count('S') > n :
            cg(i,'CG2R53')
        else : 
            cg(i,'CG2R51')


### H on cycle5
for i in cyca5 :   
    con = atom[i]['con']
    if 'H' in con and 'CG2R5' in atom[i]['cg'] :      
        cona = atom[i]['cona']   # check H-C~C-F
        het = con.count('O') or con.count('S') or con.count('N') or con.count('P')  
        if i in guans:                
            cgd(cona,{'H':'HGR53'})
        elif het:
            cgd(cona,{'H':'HGR52'})
        else:
            cgd(cona,{'H':'HGR51'})






################################################################################
################################################################################
################################################################################
################################################################################
# guanidine/guanidinium and amidinium , amides

for i in atom : 
    if not atom[i]['cg'] : 
        if shape(i) == 'C34' and got(i,atom[i]['cona'],['2N','2N+']) : 
            n = atom[i]['con'].count('N') 
            if n == 3 : 
                cg(i,'CG2N1')
            elif n == 2 : 
                cg(i,'CG2N2')
cgd(guans,{'N+':'NG2P1'})



################################################################################
################################################################################
################################################################################
################################################################################
# sp3 atoms in cycles 3,4,5

### the bridgehead C (C in bicyclic systems containing at least one 5-mem or smaller ring)
# and the aliphatic C in 3 and 4-mem cycles
bridgehead = []
for i in atom : 
    if atom[i]['a'] == 'C' and atom[i]['ncon'] == 4 : 
        cs = [cycle[c]['n'] for c in cycle if i in cycle[c]['num']]
        if len(cs) > 1 and intercept(cs,[3,4,5]) : 
            add(i,bridgehead)


for i in atom :        
    if not atom[i]['cg']: 
        if atom[i]['a'] == 'C' and atom[i]['ncon'] == 4 : 
            if i in bridgehead: 
                cg(i,'CG3RC1')                
            elif atom[i]['c'] == 3 : 
                cg(i,'CG3C31')
            elif atom[i]['c'] == 4 : 
                cg(i,'CG3C41')
            elif atom[i]['c'] == 5 : 
                con = atom[i]['con'] 
                nH = con.count('H')
                if nH and 'N+' in con : 
                    nH += 2
                cg(i,'CG3C5'+str(nH))                
        elif atom[i]['a'] == 'O': 
            if atom[i]['c'] == 3 : 
                cg(i,'OG3C31')
            elif atom[i]['c'] in [4,5] : 
                cg(i,'OG3C51')


### the amide N-C=O in 4,5-mem cycles
for c in cycle :   
    n = cycle[c]['n'] 
    if n in [4,5] :
        cyc = cycle[c]['num']    
        for i in cyc :            
            if i in carbonyl : 
                cg(i,'CG2R53')                    
                N33 = [j for j in atom[i]['cona'] if shape(j) == 'N33']
                cgs(N33,'NG2R'+str(n)+'3')


################################################################################
################################################################################
################################################################################
################################################################################
### nitro, phosphate, sulfate

# nitro
for i in atom : 
    if shape(i) in ['N+34','N35','N33','N23']: 
        cona = atom[i]['cona'] 
        O = [j for j in cona if 'O' in atom[j]['a'] and atom[j]['ncon'] == 1] 
        if len(O) > 1 : 
            cg(i,'NG2O1')
            cgs(O,'OG2N1')

# phosphate
for i in atom : 
    a , cona , con , ncon = attr (i)
    if a == 'P':
        nO = len([j for j in cona if (atom[j]['a'] in ['O-','S-'] and atom[j]['nbs'] == 1)])
        #nO += len([j for j in cona if (atom[j]['a'] in ['S'] and atom[j]['ncon'] == 2)])
        nP =  len(got(i,cona,['2N+','2C','2S','2P']))
        q = 5-atom[i]['nbs']-nO+nP
        if q > -1 :             
            cg(i,'PG0')
        elif q == -1 : 
            cg(i,'PG1')
        else : 
            cg(i,'PG2')
        O = [j for j in cona if 'O' in atom[j]['a'] and atom[j]['ncon'] == 1] 
        if 'S' in con : 
            cgs(O,'OG2S1')
        else : 
            cgs(O,'OG2P1')
    

## sulfer on phosphate  (SG2P1 and SG2P2)
for i in atom :
    if not atom[i]['cg'] and 'S' in atom[i]['a'] and atom[i]['ncon'] < 3 and 'P' in atom[i]['con']:         
        k = [j for j in atom[i]['cona'] if 'PG' in atom[j]['cg']]
        if k : 
            nf = 0
            for j in k : 
                n = atom[j]['con'].count('S')
                if n > nf : 
                    nf = n + 0                                
            atom[i]['cg'] = 'SG2P'+str(min(2,nf)) 


# sulfate
for i in atom : 
    a , cona , con , ncon = attr (i)
    if 'S' in a and (ncon > 2 or ncon == 1): 
        if ncon == 3 : 
            cg(i,'SG3O3')
        elif ncon == 4 : 
            nO = len([j for j in cona if atom[j]['a'] == 'O-'])
            nP =  len(got(i,cona,['2N+','2C','2S','2P']))
            q = 6-atom[i]['nbs']-nO+nP
            if q < 0 : 
                cg(i,'SG3O1')
            else : 
                cg(i,'SG3O2')
        O = [j for j in cona if 'O' in atom[j]['a'] and atom[j]['ncon'] == 1] 
        cgs(O,'OG2P1')
                

# -O- in phosphate and sulfate 
for i in atom :
    if atom[i]['cg'] in ['PG0','PG1','PG2','SG3O1','SG3O2','SG3O3']: 
        for j in atom[i]['cona']: 
            if atom[j]['a'] == 'O' and atom[j]['ncon'] == 2 : 
                con = atom[j]['con']
                if con in ['OPP','OPS','OPS-','OSS','OSS-']: 
                    cg(j,'OG304')
                elif 'H' in con: 
                      cg(j,'OG311')
                else : 
                    cg(j,'OG303')




         

################################################################################
################################################################################
################################################################################
################################################################################
### carbonyl and thiocarbonyl

for i in carbonyl :  
    if not atom[i]['cg']:               
        a , cona , con , ncon = attr(i) 
        k = got(i,cona,['2O','2S','2O-','2S-'])[0]
        nhet = len([j for j in cona if atom[j]['a'] in 'N+O-S-'])
        if con.count('O') + con.count('S') == 3 :
            cg(i,'CG2O6')   
        elif nhet == 3 :
            cg(i,'CG2O6')            
        elif con.count('O') + con.count('S') == 2  and '-' in con :
            cg(i,'CG2O3')
            if 'O-' in con : 
                O = [j for j in cona if shape(j) in ['O12','O-12','O-11']]
                cgs(O,'OG2D2')
        elif nhet < 3 and [j for j in cona if j !=k and 'SG3O' in atom[j]['cg']] : 
            cg(i,'CG2O3')
        elif con.count('O') + con.count('S') == 2  : 
            cg(i,'CG2O2') 
        elif 'N' in con :
            cg(i,'CG2O1') 
        elif 'H' in con : 
            cg(i,'CG2O4') 
        else : 
            cg(i,'CG2O5') 
        if 'O-' in con :
            O = [j for j in cona if shape(j) in ['O12','O-12','O-11']]
            cgs(O,'OG2D2')

  

# amides 
for i in amide : 
    nH = atom[i]['con'].count('H')
    cg(i,'NG2S'+str(nH))

# phosphoramidate
for i in atom : 
    c = atom[i]['cg']
    if 'PG' in c or (c == 'SG3O1' and 'O-' in atom[i]['con']):         
        N = [j for j in atom[i]['cona'] if shape(j) == 'N33'] 
        cgs(N,'NG2S3')




################################################################################
################################################################################
################################################################################
################################################################################
# N+

for i in atom : 
    if atom[i]['a'] == 'N+': 
        if atom[i]['ncon'] < 4 : 
            cg(i,'NG2P1')
        elif atom[i]['ncon'] == 4 : 
            nH = min(3,atom[i]['con'].count('H'))
            cg(i,'NG3P'+str(nH))
            # the attached methyls in ammonium salt 
            if not 'H' in atom[i]['con'] :
                for j in atom[i]['cona'] :                
                    if [atom[j]['a'] , atom[j]['ncon']] == ['C',4] and atom[j]['con'] == 'CHHHN+':
                        if not 'F' in atom[j]['con']:
                            cgd(atom[j]['cona'] , {'H':'HGP5'}) # should it be CH3?


################################################################################
################################################################################
################################################################################
################################################################################
# sulfur   
for i in atom:
    if not atom[i]['cg'] and 'S' in atom[i]['a']: 
        a , cona , con , ncon = attr(i)        
        if ncon == 2 :
            if 'O' in con : 
                cg(i,'SG301')                
            elif 'H' in con : 
                cg(i,'SG311')
            elif con.count('S') > 1:
                cg(i,'SG301')
            else : 
                cg(i,'SG311')
        elif ncon == 1 :            
            if got(i,cona,['2S']) :
                cg(i,'SG301')
            elif [j for j in cona if atom[j]['cg'] == 'SG3O1']: 
                cg(i,'SG302')
            elif got(i,cona,['2N+']) : 
                j = cona[0]
                if not j in benza + furas + aroma7 : 
                    cg(i,'SG311')
            if not atom[i]['cg']  :
                cg(i,'SG2D1')



for i in atom:
    if not atom[i]['cg'] and 'S' in atom[i]['a']: 
        a , cona , con , ncon = attr(i)        
        if ncon == 1 :            
            if a == 'S-' : 
                cg(i,'SG302')




################################################################################
################################################################################
################################################################################
################################################################################
### oxygens

for i in atom:
    if not atom[i]['cg'] and 'O' in atom[i]['a'] :
        a , cona , con , ncon = attr(i)
        if ncon == 1 :
            j = cona[0]
            if atom[i]['nbs'] == 1:                
                if atom[j]['a'] == 'S' and atom[j]['nbs'] == 1 : 
                    cg(i,'OG2P1')
            else : 
                if atom[j]['con'] in ['CCCO','CCCO-']: 
                    cg(i,'OG2D3')
                elif atom[j]['a'] == 'S' : 
                    cg(i,'OG2P1')                    
                else : 
                    cg(i,'OG2D1')                
        elif ncon == 2:
            nH = con.count('H')
            if nH == 2 :
                cg(i,'OGTIP3')
                cgs(cona,'HGTIP3')
            elif nH == 1 :
                cg(i,'OG311')
            elif i in ester :
                cg(i,'OG302')                
            else :
                cg(i,'OG301')


for i in atom:
    if not atom[i]['cg'] and 'O' in atom[i]['a'] :
        a , cona , con , ncon = attr(i)
        if ncon == 1 :
            j = cona[0]
            if atom[i]['nbs'] == 1:                
                cg(i,'OG312')

        
            
################################################################################
################################################################################
################################################################################
################################################################################
### C=N and C=C



for i in atom : 
    if not atom[i]['cg'] : 
        if atom[i]['a'] == 'N' and atom[i]['ncon'] == 2 : 
            cg(i,'NG2D1')

CG2D1O = []  
CG2DC1 = []
CG2D1 = []  
for i in atom :
    if not atom[i]['cg'] :
        a , cona , con , ncon = attr(i)
        if [a,ncon] == ['C',3]:
            nH = con.count('H')
            if nH == 2 :
                if i in conj :
                    cg(i,'CG2DC3')
                else :
                    cg(i,'CG2D2')
            else:
                print ( i , atom[i]['ord'] )
                k = got(i,cona,['2C','2N','2N+','2P'])[0]                                
                if [j for j in cona if atom[j]['a'] in 'O-N+S-' and j != k] :                     
                    add(i,CG2D1O)
                elif i in conj : 
                    add(i,CG2DC1)
                else :
                    add(i,CG2D1)

cgs(CG2DC1,'CG2DC1')                      
cgs(CG2D1O,'CG2D1O')                  
cgs(CG2D1,'CG2D1') 

# modefying 1-1 in conjugated systems(e.g. CG2DC1-CG2DC1)
rep = ['CG2DC1','CG2D1O','CG25C1','CG251O']
for d in conjs :  
    if atom[d[1]]['cg'] in rep and atom[d[2]]['cg'] in rep:
        atom[d[2]]['cg'] = atom[d[2]]['cg'].replace('1','2')
        if atom[d[3]]['cg'] in rep :
            atom[d[3]]['cg'] = atom[d[3]]['cg'].replace('1','2')


            
################################################################################
################################################################################
################################################################################
################################################################################
### amine


# N in cyca5
for i in cyca5 : 
    if atom[i]['a'] == 'N' : # needs mol2
        cg(i,'NG3C51') 


amine = [i for i in atom if atom[i]['a'] == 'N' and atom[i]['ncon'] == 3 and not atom[i]['cg']]
# external amine
for i in amine :
    if atom[i]['con'].count('H') == 2 and [j for j in atom[i]['cona'] if j in benza] : 
        cg(i,'NG2S3')
        cgd(atom[i]['cona'],{'H':'HGP4'})

# N in Hydrazine 
cgs([j for j in amine if atom[j]['con'].count('N') > 1],'NG3N1') 


## amine groups
amine = [i for i in amine if not atom[i]['cg']]
C44 = [i for i in atom if atom[i]['a'] == 'C' and atom[i]['ncon'] == 4]
for i in amine :
    a , cona , con , ncon = attr (i)
    if [a,ncon] == ['N', 3] and not atom[i]['cg']:
        nH = str(con.count('H'))        
        if con in ['NCCC','NCCH','NCHH','NHHH'] :
            nH = str(con.count('H'))
            cg(i,'NG3'+nH+'1')
            if not i in  allaroma :  
                cgd (cona , {'H':'HGPAM' + nH})            
        else :
            atom[i]['cg'] = 'NG3' + nH + '1'     
        # the attached methyls in special cases of methylamine
        exc = [j for j in cona if not atom[j]['con'] in ['CHHHN' , 'HN']]
        if not exc :
            cgc = 'CG3AM' + nH
            cgh = 'HGAAM' + nH
            cgd (cona , {'C':cgc})
            for j in cona :
                cgd (atom[j]['cona'] , {'H':cgh})
                


################################################################################
################################################################################
################################################################################
################################################################################
### aliphatic carbon

for i in atom :
    if not atom[i]['cg'] :
        if atom[i]['a'] == 'C' and atom[i]['ncon'] == 4 :            
            a , cona , con , ncon = attr(i)
            nH = con.count('H')            
            # the following priority is taken from paramchem
            if 'N+' in con and nH :
                cg(i,'CG3'+str(nH)+'4') 
            if 'S-' in con and nH :
                cg(i,'CG323')
            if not atom[i]['cg'] and 'F' in con :                
                d = {1:'CG322', 2:'CG312', 3:'CG302'}
                nF = con.count('F')
                cg(i,d[nF])
            if not atom[i]['cg'] :
                nH = min(3,nH)
                cg(i,'CG3'+str(nH)+'1')


################################################################################
################################################################################
################################################################################
################################################################################
### halos
for i in atom : 
    if atom[i]['a'] == 'I': 
        cg(i,'IGR1')

for i in atom :
    a = atom[i]['a']
    if a in ['F','L','B','I'] :        
        d = {'L' : {1:'CLGA1' , 2:'CLGA1' , 3:'CLGA3' , 4:'CLGA3'} ,
             'B' : {1:'BRGA1' , 2:'BRGA2' , 3:'BRGA3' , 4:'BRGA3'} ,
             'F' : {1:'FGA1'  , 2:'FGA2'  , 3:'FGA3'  , 4:'FGA3' } }   
        if not atom[i]['cona'] :
            atom[i]['cg'] = d[a][1]
        else :
            j = atom[i]['cona'][0]
            aj , cona , con , ncon = attr(j)
            if a == 'F' and aj != 'C':                
                cg(i,'FGP1')                
            elif j in benza or j in furas or j in aroma7:                
                cgd([i],{'F':'FGR1' , 'L':'CLGR1' , 'B':'BRGR1' , 'I':'IGR1'})
            else :
                nX = len([k for k in cona if atom[k]['a'] in 'FLBI'])
                cgd([i] , {'F':d['F'][nX] , 'L':d['L'][nX], 'B':d['B'][nX]})


################################################################################
################################################################################
################################################################################
################################################################################
### misc

for i in atom :
    a = atom[i]['a']
    d = {'A':'ALG1' , 'D':'DUM' , 'U':'DUMM' , 'M':'HE' , 'E':'NE'}
    if a in d:
        cg(i,d[a])
  
################################################################################
################################################################################
################################################################################
################################################################################
### Hydrogens
for i in atom :
    if atom[i]['a'] == 'H' and not atom[i]['cg']:
        j = atom[i]['cona'][0]
        a , cona , con , ncon = attr(j)
        if [a,ncon] == ['C',3] :
            if got(j,cona,['2N+']) or  j in carbonyl: 
                cg(i,'HGR52')
                
for i in atom :
    if atom[i]['a'] == 'H' and not atom[i]['cg']:
        j = atom[i]['cona'][0]
        a , cona , con , ncon = attr(j)
        if [a,ncon] == ['C',3] :                
            if got(j,cona,['2N','2C']):
                nH = con.count('H')
                if nH == 1 :
                    cg(i,'HGA4') 
                elif nH == 2:
                    cgd(cona,{'H':'HGA5'})                                    

        #elif atom[j]['cg'] == 'NG2S3': 
        #    for k in cona : 
        #        if '2R6' in atom[k]['cg']: 
        #            cg(i,'HGP4')
        #    if not atom[i]['cg']:
        #        cg(i,'HGP1')    
                    
        elif a == 'S' :
            cg(i,'HGP3')
        elif a == 'N+': 
            cg(i,'HGP2')      
        elif a != 'C' :
            cg(i,'HGP1')
        elif a == 'C' and ncon == 4 :             
            nF = con.count('F') 
            if nF : 
                n = str(nF+5)
                atom[i]['cg'] = 'HGA'+n
            else : 
                nH = con.count('H') 
                if nH == 1 : 
                    atom[i]['cg'] = 'HGA1' 
                elif nH == 2 : 
                    atom[i]['cg'] = 'HGA2'
                else : 
                    atom[i]['cg'] = 'HGA3' 


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
'''
errors = []
if 1 in typs : 
    for i in atom :
        onc = atom[i]['cg']
        pc = typs[i]    
        #print ( i , onc , pc , onc == pc
        match = ''
        if onc != pc :
            match = pc    
        n = len(onc)
        if n < 6 :
            cn = 6-n
            onc += cn*' '   
        oncpc = onc+pc
        if pc == onc.replace('1O','2O') or pc == onc.replace('2O','1O') or \
            pc == onc.replace('C1','C2') or pc == onc.replace('C2','C1') : 
                match = ''        
        if match and penalty[i] < 50000 :
            shapes = shape(i) +'\t'+ ''.join(sorted([shape(j) + '\t' for j in atom[i]['cona']]))
            errors += [str(i) +'\t'+ onc +'\t'+ pc + '\t' + str(int(penalty[i])) + '    ' + shapes +  '\n']
            #errors += [str(i) +'\t'+ onc +'\t'+ pc + '\n']

print ( '==============='
if errors : 
    for r in errors : 
        print ( r.strip('\n')

def err() : 
    if errors :   
        tf = open('errors','a')          
        tf.write('==================== '+srt+' ==================\n')
        for r in errors : 
            tf.write(r)
        tf.write('\n\n\n')
    
        tf.close()

err()

'''


cgenff = {
'HGA1':	[-0.045,	1.34	,	-0.045,	1.34  ],
'HGA2':	[-0.035,	1.34	,	-0.035,	1.34  ],
'HGA3':	[-0.024,	1.34	,	-0.024,	1.34  ],
'HGA4':	[-0.031,	1.25	,	-0.031,	1.25  ],
'HGA5':	[-0.026,	1.26	,	-0.026,	1.26  ],
'HGA6':	[-0.028,	1.32	,	-0.028,	1.32  ],
'HGA7':	[-0.03,	1.3	,	-0.03	,	1.3   ],
'HGAAM0':	[-0.028,	1.28	,	-0.028,	1.28  ],
'HGAAM1':	[-0.028,	1.28	,	-0.028,	1.28  ],
'HGAAM2':	[-0.04,	1.26	,	-0.04	,	1.26  ],
'HGP1':	[-0.046,	0.2245,	-0.046,	0.2245],
'HGP2':	[-0.046,	0.2245,	-0.046,	0.2245],
'HGP3':	[-0.1	,	0.45	,	-0.1	,	0.45  ],
'HGP4':	[-0.046,	0.2245,	-0.046,	0.2245],
'HGP5':	[-0.046,	0.7	,	-0.046,	0.7   ],
'HGPAM1':	[-0.009,	0.875	,	-0.009,	0.875	],
'HGPAM2':	[-0.01,	0.875	,	-0.01	,	0.875	],
'HGPAM3':	[-0.012,	0.87	,	-0.012,	0.87	],
'HGR51':	[-0.03,	1.3582,	-0.03	,	1.3582],
'HGR52':	[-0.046,	0.9	,	-0.046,	0.9	],
'HGR53':	[-0.046,	0.7	,	-0.046,	0.7	],
'HGR61':	[-0.03,	1.3582,	-0.03	,	1.3582],
'HGR62':	[-0.046,	1.1	,	-0.046,	1.1	],
'HGR63':	[-0.046,	0.9	,	-0.046,	0.9	],
'HGR71':	[-0.03,	1.3582,	-0.03	,	1.3582],
'!HGTIP3':	[-0.046,	0.2245,	-0.046,	0.2245],
'CG1T1':	[-0.167,	1.84	,	-0.167,	1.84	],
'CG1T2':	[-0.1032,	1.9925,	-0.1032,	1.9925],
'CG1N1':	[-0.18,	1.87	,	-0.18	,	1.87	],
'CG2D1':	[-0.068,	2.09	,	-0.068,	2.09	],
'CG2D2':	[-0.064,	2.08	,	-0.064,	2.08	],
'CG2D1O':	[-0.068,	2.09	,	-0.068,	2.09	],
'CG2D2O':	[-0.068,	2.09	,	-0.068,	2.09	],
'CG2DC1':	[-0.068,	2.09	,	-0.068,	2.09	],
'CG2DC2':	[-0.068,	2.09	,	-0.068,	2.09	],
'CG2DC3':	[-0.064,	2.08	,	-0.064,	2.08	],
'CG2N1':	[-0.11,	2	,	-0.11	,	2	],
'CG2N2':	[-0.11,	2	,	-0.11	,	2	],
'CG2O1':	[-0.11,	2	,	-0.11	,	2	],
'CG2O2':	[-0.098,	1.7	,	-0.098,	1.7	],
'CG2O3':	[-0.07,	2	,	-0.07	,	2	],
'CG2O4':	[-0.06,	1.8	,	-0.06	,	1.8	],
'CG2O5':	[-0.09,	2	,	-0.09	,	2	],
'CG2O6':	[-0.07,	2	,	-0.07	,	2	],
'CG2O7':	[-0.058,	1.563	,	-0.058,	1.563	],
'CG2R51':	[-0.05,	2.1	,	-0.05	,	2.1	],
'CG2R52':	[-0.02,	2.2	,	-0.02	,	2.2	],
'CG2R53':	[-0.02,	2.2	,	-0.02	,	2.2	],
'CG2R57':	[-0.05,	2.1	,	-0.05	,	2.1	],
'CG25C1':	[-0.068,	2.09	,	-0.068,	2.09	],
'CG25C2':	[-0.068,	2.09	,	-0.068,	2.09	],
'CG251O':	[-0.068,	2.09	,	-0.068,	2.09	],
'CG252O':	[-0.068,	2.09	,	-0.068,	2.09	],
'CG2R61':	[-0.07,	1.9924,	-0.07	,	1.9924],
'CG2R62':	[-0.09,	1.9	,	-0.09	,	1.9	],
'CG2R63':	[-0.1	,	1.9	,	-0.1	,	1.9	],
'CG2R64':	[-0.04,	2.1	,	-0.04	,	2.1	],
'CG2R66':	[-0.07,	1.9	,	-0.07	,	1.9	],
'CG2R67':	[-0.07,	1.9924,	-0.07	,	1.9924],
'CG2RC0':	[-0.099,	1.86	,	-0.099,	1.86	],
'CG2R71':	[-0.067,	1.9948,	-0.067,	1.9948],
'CG2RC7':	[-0.099,	1.86	,	-0.099,	1.86	],
'CG301':	[-0.032,	2	,	-0.01	,	1.9	],
'CG302':	[-0.02,	2.3	,	-0.02	,	2.3	],
'CG311':	[-0.032,	2	,	-0.01	,	1.9	],
'CG312':	[-0.042,	2.05	,	-0.042,	2.05	],
'CG314':	[-0.031,	2.165	,	-0.01	,	1.9	],
'CG321':	[-0.056,	2.01	,	-0.01	,	1.9	],
'CG322':	[-0.06,	1.9	,	-0.06	,	1.9	],
'CG323':	[-0.11,	2.2	,	-0.11	,	2.2	],
'CG324':	[-0.055,	2.175	,	-0.01	,	1.9	],
'CG331':	[-0.078,	2.05	,	-0.01	,	1.9	],
'CG334':	[-0.077,	2.215	,	-0.01	,	1.9	],
'CG3C50':	[-0.036,	2.01	,	-0.01	,	1.9	],
'CG3C51':	[-0.036,	2.01	,	-0.01	,	1.9	],
'CG3C52':	[-0.06,	2.02	,	-0.01	,	1.9	],
'CG3C53':	[-0.035,	2.175	,	-0.01	,	1.9	],
'CG3C54':	[-0.059,	2.185	,	-0.01	,	1.9	],
'CG3C31':	[-0.056,	2.01	,	-0.01	,	1.9	],
'CG3C41':	[-0.065,	2.02	,	-0.01	,	1.9	],
'CG3RC1':	[-0.032,	2	,	-0.01	,	1.9	],
'CG3AM0':	[-0.07,	1.97	,	-0.07	,	1.97	],
'CG3AM1':	[-0.078,	1.98	,	-0.078,	1.98	],
'CG3AM2':	[-0.08,	1.99	,	-0.08	,	1.99	],
'NG1T1':	[-0.18,	1.79	,	-0.18	,	1.79	],
'!NG1D1':	[-0.035,	2.03	,	-0.035,	2.03	],
'NG2D1':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG2S0':	[-0.2	,	1.85	,	-0.0001,	1.85	],
'NG2S1':	[-0.2	,	1.85	,	-0.2	,	1.55	],
'NG2S2':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG2S3':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG2O1':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG2P1':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG2R43':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG2R50':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG2R51':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG2R52':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG2R53':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG2R57':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG2R60':	[-0.06,	1.89	,	-0.06	,	1.89	],
'NG2R61':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG2R62':	[-0.05,	2.06	,	-0.05	,	2.06	],
'NG2R67':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG2RC0':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG301':	[-0.035,	2	,	-0.035,	2	],
'NG311':	[-0.045,	2	,	-0.045,	2	],
'NG321':	[-0.06,	1.99	,	-0.06	,	1.99	],
'NG331':	[-0.07,	1.98	,	-0.07	,	1.98	],
'NG3C51':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG3N1':	[-0.06,	2.05	,	-0.06	,	2.05	],
'NG3P0':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG3P1':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG3P2':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'NG3P3':	[-0.2	,	1.85	,	-0.2	,	1.85	],
'OG2D1':	[-0.12,	1.7	,	-0.12	,	1.4	],
'OG2D2':	[-0.12,	1.7	,	-0.12	,	1.7	],
'OG2D3':	[-0.05,	1.7	,	-0.12	,	1.4	],
'OG2D4':	[-0.12,	1.7	,	-0.12	,	1.7	],
'OG2D5':	[-0.165,	1.692	,	-0.165,	1.692	],
'OG2N1':	[-0.12,	1.7	,	-0.12	,	1.7	],
'OG2P1':	[-0.12,	1.7	,	-0.12	,	1.7	],
'OG2R50':	[-0.12,	1.7	,	-0.12	,	1.7	],
'OG3R60':	[-0.1	,	1.65	,	-0.1	,	1.65	],
'OG301':	[-0.1	,	1.65	,	-0.1	,	1.65	],
'OG302':	[-0.1	,	1.65	,	-0.1	,	1.65	],
'OG303':	[-0.1	,	1.65	,	-0.1	,	1.65	],
'OG304':	[-0.1	,	1.65	,	-0.1	,	1.65	],
'OG311':	[-0.1921,	1.765	,	-0.1921,	1.765	],
'OG312':	[-0.12,	1.75	,	-0.12	,	1.75	],
'OG3C31':	[-0.1	,	1.65	,	-0.1	,	1.65	],
'OG3C51':	[-0.1	,	1.65	,	-0.1	,	1.65	],
'OG3C61':	[-0.1	,	1.65	,	-0.1	,	1.65	],
'!OGTIP3':	[-0.1521,	1.7682,	-0.1521,	1.7682],
'SG2D1':	[-0.565,	2.05	,	-0.565,	2.05  ],
'SG2R50':	[-0.45,	2	,	-0.45	,	2	],
'SG311':	[-0.45,	2	,	-0.45	,	2	],
'SG301':	[-0.38,	1.975	,	-0.38	,	1.975	],
'SG302':	[-0.47,	2.2	,	-0.47	,	2.2	],
'SG3O1':	[-0.47,	2.1	,	-0.47	,	2.1	],
'SG3O2':	[-0.35,	2	,	-0.35	,	2	],
'SG3O3':	[-0.35,	2	,	-0.35	,	2	],
'FGA1':	[-0.135,	1.63	,	-0.135,	1.63	],
'FGA2':	[-0.105,	1.63	,	-0.105,	1.63	],
'FGA3':	[-0.097,	1.6	,	-0.097,	1.6	],
'FGP1':	[-0.097,	1.6	,	-0.097,	1.6	],
'FGR1':	[-0.12,     1.7	,	-0.12	,	1.7	],
'CLGA1':	[-0.343,	1.91	,	-0.343,	1.91	],
'CLGA3':	[-0.31,	1.91	,	-0.31	,	1.91	],
'BRGA1':	[-0.48,	1.97	,	-0.48	,	1.97	],
'BRGA2':	[-0.53,	2.05	,	-0.53	,	2.05	],
'BRGA3':	[-0.54,	2	,	-0.54	,	2	],
'!DUM':	[0	,	0	,	0	,	0	],
'!HE'	:	[-0.02127,	1.48	,	-0.02127,	1.48	],
'!NE'	:	[-0.08545,	1.53	,	-0.08545,	1.53	],
'PG0'	:	[-0.585,	2.15	,	-0.585,	2.15	],
'PG1'	:	[-0.585,	2.15	,	-0.585,	2.15	],
'PG2'	:	[-0.585,	2.15	,	-0.585,	2.15	],
'ALG1':	[-0.65,	2	,	-0.65	,	2	],
'LPH'	:	[0	,	0	,	0	,	0	],
'CLGR1':	[-0.23,	1.86	,	-0.23	,	1.86	],
'BRGR1':	[-0.32,	1.98	,	-0.32	,	1.98	],
'IGR1':	[-0.52,	2.24	,	-0.52	,	2.24	],
'SG2P1':	[-0.6308,	2.0937,	-0.6308,	2.0937],
'OG2S1':	[-0.1423,	1.6796,	-0.1423,	1.6796],
'SG2P2':	[-0.6199,	2.0546,	-0.6199,	2.0546]
}

print ( 'Lennard Jones Parameters (atom cgenff_type epsilon Rmin epsilon4 Rmin4)' )
### writing the parameters in a temporary output file 
tf = open('tempars','a')
for i in atom :    
    [a,b,c,d] = cgenff[atom[i]['cg']]
    print ( i , end=' ' )
    print (('%8s%8s%8s%8s%8s'%(atom[i]['cg'],a,b,c,d)) )
    tf.write('lj%5s%8s%8s%8s%8s\n'%(i,a,b,c,d)) 
tf.close()



print ()
print ( 'Lennard-Jones Parameterization is done  __________________________________________' )
print ( '__________________________________________________________________________________' )


tf = open('ljdone','w')
tf.close()