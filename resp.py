#!/usr/bin/env python

####=================================================================================================
#      ________      ________      ________       ________
#      \  ____ \     \  _____\     \  _____\      \  ____ \
#       \ \   \ \     \ \           \ \            \ \   \ \
#        \ \___\ \     \ \____       \ \_____       \ \___\ \
#         \  __  _\     \  ___\       \_____ \       \  _____\
#          \ \ \ \       \ \                \ \       \ \
#           \ \ \ \       \ \______      ____\ \       \ \
#            \_\ \_\       \_______\     \______\       \_\
#
####=================================================================================================
print ( '\n\n' )
print ( '__________________________________________________________________________________' )
print ( 'Charge Fitting   _________________________________________________________________' )
print ()

from os import system
import sys
from functions import *

#==========================================================================
#==========================================================================

espf = sys.argv[1]
cfg = sys.argv[-1]
if cfg == '0' : 
    cfg = 0
'''
espf = 'nma.esp'
cfg = 'cfg'
'''

print ( '\nThe Electrostatic Potential file: {}  was detected\n'.format(espf) )

print ( 'STEP: extracting the input data for RESP' )


atoms , bonds = atoms_bonds()
nat = len(atoms)
uniqs = {i:[] for i in atoms}
if atoms: 
    uniqs = get_symmetries(atoms)
nat = len(atoms)


## extracting the cartesian coordinate of the charge centers and the ESP points and the electrostatic potential  
sf = open(espf)
for line in sf : 
    if 'Charge =' in line : 
        Q = float(line.split()[2])
        break

for line in sf :
    if 'Electrostatic Properties Using The SCF Density' in line :
        sf.readline() ; sf.readline() ; sf.readline()
        break

acc = {}     # the cartesian coordinates of the atomic centers
for i in atoms :
    line = sf.readline().split()        
    center = [float(line[-3]) , float(line[-2]), float(line[-1])]
    acc[i] = center

ecc = []     # the cartesian coordinates of the esp points
for line in sf :
    if 'ESP Fit Center' in line :
        line = line.split()
        ecc += [[float(line[-3]) , float(line[-2]) , float(line[-1])]]
    else :
        break

esp = []      # the values of the esp point
for line in sf :
    if 'Electrostatic Properties (Atomic Units)' in line :
        for line in sf :
            if 'Fit' in line :
                line = line.split()
                esp += [float(line[2])]
            elif 'Axes restored to original set' in line:
                break
print ( '\n{} electrostatic potential points has been deteted\n'.format(len(esp)) )
sf.close()


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


equic = []
frozen = {}
frozeng = []
keep = []
nosymm = 0
if cfg : 
    sf = open('cfg')
    for line in sf :
        if not line.startswith('#') :
            if 'nosymmetry' in line : 
                nosymm = 1
            line = line.split()
            n = len(line)
            if n > 1 : 
                if 'fc' in line : 
                    frozen[int(line[1])] = float(line[2])
                elif 'gc' in line : 
                    line = line[1:]
                    frozeng += [[int(i) for i in line[:-1]] + [float(line[-1])]]
                elif 'ec' in line and n > 2: 
                    line = line[1:]
                    equic += [sorted([int(i) for i in line])]
                elif 'kc' in line :
                    line = line[1:]
                    keep = [int(i) for i in line]
    
    sf.close()
frozeng += [[i for i in atoms] + [Q]]

sf.close()


if nosymm : 
    uniqs = {i:[] for i in atoms}
 
print ()
print ( 'The program equivalencing scheme has detected the following equivalent charge centers ' )
for i in uniqs :
    if uniqs[i] : 
        eq = [i] + uniqs[i]
        print ( eq )
print ( '\n' )
   

#!!!!! here we need to define the methyl hydrogen because it should be allowed for fitting in the second stage 


'''
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




#==========================================================================
#==========================================================================     
### the constraints
# ther is four types of constraints
# 1 symmetric atoms
# 2 two atoms constrained to have the same value 
# 3 an atom constrained to particular charge
# 4 a group of atoms constrained to have a total charge


equic = [] # atoms constrained to have the same value during fittting
                # provided by the user in successive lines like that :
                # 1 2
                # 4 7 9
                                                               
symmetric = {} # determines the atoms which are constrained to have
               # the same charge (detected by the program)

frozen = {}   # the atomic centers which are constrained to have a fixed value
              # for stage1 they are provided by the user in the configuaration file
              # they are provided as successive lines of 'atom charge' like
              # fc 1 -0.25563
              # fc 2  0.58742
              # for stage2 they are the hetero atoms and the bound hydrogens in addition to those provided by the user
              # in either satge frozen atoms will be excluded  

frozeng = [[i for i in atoms] + [Q]]
#frozeng += [[1 , 2 , 0.9]]
#frozeng = []
### the frozen groups
# groups of atoms are constrained to have a fixed summ for their charges
# these groups provided in the configuration file if required
# they are provided as successive lines of 'atom1  atom2  atom3 ... charge' like
# 1 2 0.2596
# 1 3 4 -0.154875
# the total charge as the sum of all charges is the default element iprint ( 'equivalent = ' , equic
print ( 
print ( 'frozen charges = ' , frozen
print ( 
print ( 'group charges = ' , frozengn the frozengroups  
'''
def modcons (samc , symm , froz , frozg) :
    # modefying the constraints if ther is any of cnflicts or additional equivalences
    # the you have to stick to the following sequence 

    nfrozg = len(frozg)

    ### disabeling the equivalence of the unequal frozen atoms 
    for i in froz :
        for j in froz :
            if i < j and froz[i] != froz[j]: 
                if i in symm and j in symm[i] : 
                    k  = min([a for a in symm[i] if not a == j])
                    ks = [a for a in symm[i] if not a in [j,k]]
                    symm[k] = ks
                    symm[i] = []
                    symm[j] = []
                    print ( 'the equivalence between the charge centers {} and {} is disabled to keep them frozen'.format(i,j) )
                else : 
                    for k in symm : 
                        if i in symm[k] and j in symm[k] : 
                            symm[k] = [a for a in symm[k] if not a in [i,j]]
                            print ( 'the equivalence between the charge centers {} and {} is disabled to keep them frozen'.format(i,j) )
                            break 


    ### disabeling the frozen groups whose charge centers are individually frozen  
    frozgi = [f for f in frozg]
    for f in frozgi : 
        g = f[:-1]
        z = [i for i in g if i in froz]
        if len(z) == len(g) : 
            frozg.remove(f)
            print ( 'The frozen group charge {} is disabled to keep its individual frozen charge centers in effect'.format(f) )


    ### disabeling the equivalence of the unequal frozen groups 
    frozgi = [f for f in frozg]
    def equigroups(g1,g2) : 
        for s in symm : 
            if len(intercept(g1,[s]+symm[s])) != len(intercept(g2,[s]+symm[s])) : 
                return 0
        return 1            
    
    for i in frozgi :
        for j in frozgi :
            if j > i and len(i) == len(j) and i[-1] != j[-1] : 
                if equigroups(i,j) :
                    print ( 'The equivalences within the group charges {} and {} are disabled to keep the frozen group charges in effect'.format(i,j) )
                    for s in symm : 
                        sij = intercept(i , [s]+symm[s]) + intercept(j , [s]+symm[s]) 
                        m = min([k for k in [s]+symm[s] if not k in sij])
                        ms = [k for k in [s]+symm[s] if not k in [m]+sij]
                        symm[m] = ms
                        if s != m : 
                            symm[s] = []
    nfrozg = len(frozg)



    ### appending the 'equic' atoms to the symmetric atoms
    symmj = [i for i in symm]
    for i in samc :    
        for j in symmj :
            if intercept(i,[j]+symm[j]): 
                symm[j] += [k for k in i if not k in [j]+symm[j]]
                

    # joining the intercepting symmetries
    s = [i for i in symm]
    s.reverse()
    for i in s :
        for j in s :
            if j < i :
                if intercept([i]+symm[i] , [j]+symm[j]):
                    symm[j] += [k for k in [i]+symm[i] if not k in symm[j] ]
                    symm[i] = [0]
    symm = {i:[j for j in symm[i] if j != 0 and j != i] for i in symm if symm[i] != [0]}


    ### freezing the atoms equivalent to frozen atoms
    fro = [f for f in froz]
    for f in fro :
        for s in symm :        
            if f in [s]+symm[s] : 
                for k in [s]+symm[s] : 
                    froz[k] = froz[f]    
                break 

    return symm , froz , frozg , nfrozg

#______________________________________________________________________________
#______________________________________________________________________________
### preparing A and B to solve Aq = B for q

print ( 'STEP: input data processing ' )
def distesp(i,j) :  return 1/dist(acc[i],ecc[j])
def dod (v,w): return sum([v[i]*w[i] for i in v])
nesp = len(esp)
nat1 = nat+1



#atoms1 = [i for i in range(1,nat1)]
R = {i:{j+1:distesp(i,j) for j in range(nesp) } for i in atoms} # Rij the reciprocal of the distance between the jth

# building the A0 and B0 
B0 = {i:dot(esp,lstd(R[i])) for i in atoms}   #  
A0 = {i:{j:dod(R[i],R[j]) for j in atoms if j >= i} for i in atoms}
for i in atoms :
    for j in atoms :
        if j > i :
            A0[j][i] = A0[i][j] 


def A0B0(A0 , B0 , frozg , sym) :
    nfrozg = len(frozg)
    ### adding the constraints to the A0 and B0    
    for i in range(nfrozg):
        B0[nat1+i] = frozg[i][-1]
        A0[nat1+i] = {}
        sym [nat1+i] = []
        for j in atoms :
            if j in frozg[i] :
                A0[nat1+i][j] = 1
                A0[j][nat1+i] = 1
            else :
                A0[nat1+i][j] = 0
                A0[j][nat1+i] = 0

    for i in range(nfrozg) :
        for j in range(nfrozg) :
            if i == j :
                A0[nat1+i][nat1+j] = 1e-10
            else :
                A0[nat1+i][nat1+j] = 0
    return A0 , B0 , sym

### now the shape of the A0 and B0 is
'''
    _                                  _    _    _         _    _
   |  A11    A12 ...... A1n   1   1   1 |  |  q1  |       |  B1  |
   |  A21    A22 ...... A2n   1   1   0 |  |  q2  |       |  B2  |
   |  .      .   .      .     1   0   1 |  |  .   |       |  .   |
   |  .      .     .    .     1   0   1 |  |  .   |       |  .   |
   |  .      .       .  .     1   0   0 |  |  .   |   =   |  .   |
   | An1    An2 ......  Ann   1   0   0 |  |  qn  |       |  Bn  |
   |  1      1  1 1 1 1 1     E   0   0 |  |  L   |       |  Q   |
   |  1      1  0 0 0 0 0     0   E   0 |  |  L   |       |  Q1  |
   |  1      0  1 1 0 0 0     0   0   E |  |  L   |       |  Q2  |
   |_                                  _|  |_    _|       |_    _|
   
  
'''
### calculating A and B with eliminating the frozen charges and adding the symmetries  
def AB (Ar,Br,froz,sym,s): 
    cod = {}
    k = 0
    for i in Br : 
        if i in froz : 
            cod[i] = -1
        elif i in sym : 
            cod[i] = k
            for j in sym[i] : 
                cod[j] = k
            k += 1
    codes = [i for i in range(k)]
    B = [0 for i in codes]
    A = [[0 for j in codes] for i in codes]
    for i in Br : 
        ci = cod[i]
        if ci >= 0 : 
            B[ci] += Br[i]
            for j in Br : 
                cj = cod[j]
                if cj >= 0 : 
                    A[ci][cj] += Ar[i][j]
                else : 
                    B[ci] -= froz[j]*Ar[i][j]
    return A , B 


### adding the restraints to the diagonal of A (no b for hyperbolic)
def restraint (A0,B0,qf,s,step): 
    res = {}
    for i in atoms : 
        if atoms[i]['a'] == 'H' :
            res[i] = 0.0
        else : 
            if step == 0 : 
                res[i] = s
            else : 
                res[i] = s/sqrt(qf[i-1]**2+0.01)
    Ar = {i:{j:A0[i][j] for j in A0} for i in A0}
    Br = {i:B0[i] for i in B0}  
    for i in res :    
        Ar[i][i] += res[i]
        if step == 0 : 
            Br[i] += res[i]*qf[i-1]
    return Ar , Br 



def qall(qf,symm , froz):
    qn = {}
    k = 0    
    for i in atoms :
        if i in froz :
            qn[i] = froz[i]
        elif i in symm :
            qn[i] = qf[k]
            for j in symm[i] :
                qn[j] = qf[k]
            k += 1
    qn = [qn[i] for i in qn]
    return qn
                

   
def fit(A0,B0,qi,sym,froz,s):
    nq = max(atoms)
    qo =  [i for i in qi]
    qf =  [i for i in qi]
    for step in range(101) :
        Ar , Br =  restraint (A0 , B0 , qf ,s,step)
        A , B = AB (Ar,Br,froz,sym,s)
        qf = LSQR(A,B)        
        qf = qall(qf,sym,froz)
        conv = sqrt(sum([(qo[i-1]-qf[i-1])**2 for i in atoms]))/nq        
        if  conv < 1e-5 :
            print ( 'resp fitting has converged in  '  , step+1 , ' iterations' )       
            return qf
        elif step == 100 :
            print ( 'no convergence after '  , step+1 , ' iterations' )
            return qf 
        else :
            qo = [q for q in qf]
        

# restraint options
# 1 - make restraints = 0 for hydrogen atoms and non zero for heavy atoms (selected)
# 2 - make restraints non zero for all atoms 


#_______________
# first stage
# frozen atoms and frozen groups are in effect 
# equics and symmetry are in effect
print ( '\nSTEP: first stage fitting' )
symm = {i:[] for i in atoms}
samc = []
froz = {i:frozen[i] for i in frozen}

if froz : 
    print ()
    print ( 'According to the configuration file, the following atoms are frozen ' )
    print ( froz )

symm , froz , frozg , nfrozg = modcons (samc , symm , froz , frozeng)

A0 , B0 , symm  = A0B0 (A0 , B0 , frozeng , symm)

q0 = [] 
for i in atoms :
    if i in frozen : 
        q0 += [frozen[i]]
    else : 
        q0 += [0]
           
qf = fit (A0,B0,q0,symm,frozen,0.0005)


#_______________
# second stage
# all constraints are in effect
# all atoms except H and C are frozen in this stage
print ( '\nSTEP: second stage fitting' )
symm = uniqs
samc = [s for s in equic]
froz = {i:qf[i-1] for i in atoms if i in frozen or i in keep or not (atoms[i]['a'] == 'C' or atoms[i]['con'] == 'HC')}  # make sure whether methyl 
symm , froz , frozg , nfrozg = modcons (samc , symm , froz , frozeng)

A0 , B0 , symm  = A0B0 (A0 , B0 , frozeng , symm)

qf = fit (A0,B0,qf,symm,froz,0.001)


print ()
print ( 'According to the  the configuration file the program has modified the equivalences to be ' )
for i in symm :
    if symm[i] : 
        eq = [i] + symm[i]
        print ( eq )

print ()
print ( 'With program modification scheme, the frozen atoms were corrected to be ' )
print ( froz )

print ( 'and the frozen group charge were corrected to be ' )
for f in frozg :
    print ( f )


print ()
print ( 'The final charges ' )
k = 1
for q in qf :
    print ( k , q )    
    k += 1
print () 
tf = open('tempars','a')
tf.write('charges ')
for q in qf :
    tf.write(str(round(q,4)) + '   ')
tf.write('\n')
tf.close()

'''
tf = open('charges.txt','w')
for q in qf :
    tf.write(str(round(q,4)) + '   ')
tf.write('\n')
tf.close()
'''

print ()
print ( 'Charge Fitting is done  __________________________________________________________' )
print ( '__________________________________________________________________________________' )


tf = open('espdone','w')
tf.close()