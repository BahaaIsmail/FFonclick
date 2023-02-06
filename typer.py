#!/usr/bin/python


#   _______________                   _________
#   \______________\   |\     |\     |\ _______\
#        | |           | |    | |    | |     | |
#        | |           | |____| |    | |_____| |
#        | |            \|_____\|    | |______\|
#        | |                  | |    | |      
#        | |           |\     | |    | |      
#        | |           | |____| |    | |      
#         \|            \|_____\|     \|      
                   

### this typer containes the tested and verified fixar.py
### find it in the directory "/home/bahaa/work_fix/fixar_tested_and_verified_2017-11-11/fixar.py"

### typer 29-8-2016

def cg(i,x):
    # assigns the atom type provided it is not already assigned
    if not atom[i]['cg']:
        atom[i]['cg'] = x


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

from functions import *

atom , bond = atoms_bonds()
angle = get_angles(atom)
dihedral = get_dihedrals(atom,bond)
cycle = get_cycles(atom,bond,angle,dihedral)
uniqs = get_symmetries(atom)
uniq = [i for i in uniqs]


# modefying Cl, Br, Al, He , Ne and DUM to not be conflicted
# these names will be corected at the end of the typing steps
d = {'Cl':'L','Br':'B','Al':'A','He':'M','Ne':'E','DUM':'D','DUMM':'U','Se':'S', 
     'CL':'L','BR':'B','AL':'A','HE':'M','NE':'E'}
for i in atom :
    a = atom[i]['a']
    if a in d :
        atom[i]['a'] = d[a]

### ===================================
### more attributes for atom dictionary
for i in atom :
    if atom[i]['ord'] :
        atom[i]['o'] = max(atom[i]['ord'])   # the max order
    nbs = sum(atom[i]['ord'])     # make sure    
    if int(nbs) == nbs : 
        nbs = int(nbs)
    atom[i]['nbs'] = round(nbs,1) 
    atom[i]['con'] = conic(i , atom)
    atom[i]['c'] = 0
    atom[i]['cg'] = ''


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
### aromatic and planar rings


# rings (not containing ncon = 4 or atoms other than C,N,O,S,P  )

ring = {c:{'r':0 , 'h':[], 'd':[], 'x':[]} for c in cycle if cycle[c]['n'] in [5,6,7]}
for c in ring :
    itisring = 1
    cyc = cycle[c]['num']
    for i in cyc :        
        a , cona , con , ncon = attr (i,atom)
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
            a , cona , con , ncon = attr (i,atom)
            ancon = shape(i,atom)
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

################################################################################
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
        if d + 2*h + x < 6 or d + h > 6: 
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
#print 'permanent aromatics = ' ,  aromatics

        
cand = [c for c in aromatics if not aromatics[c]]
#print 'cand = ' , cand 
for itr in range(99) : 
    #print itr , '_______________________'
    updated = 0    
    for c in cand :
        getaromatic = 0
        cyc = cycle[c]['num']
        x = ring[c]['x']
        d = len(ring[c]['d'])
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
                    elif i in x :
                        add(i,r)             
        
        a = len(a)
        h = len(h)-a
        r = len(r)
        x = len(x)-r
        base = d + 2*h            
        npi = base + a + r
        #print c , d , 2*h , a , r , npi
        if npi + a >= 6 : 
            if base <= 6 and h+a < 3: 
                getaromatic = 1
            elif npi <= 6 : 
                getaromatic = 1
            elif npi-r < 6 or npi+a-r < 6 : 
                getaromatic = 1
        if getaromatic != aromatics[c] : 
            aromatics[c] = getaromatic
            updated = 1
    if not updated : 
        break

aromatics = [c for c in aromatics if aromatics[c]]    
#print 'aromatics : ' , aromatics
  

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
plana = [i for i in planas if shape(i,atom) == 'N33']


# adding "+" to positive nitrogenes and '-' to the negative oxygens and sulfurs
for i in atom : 
    if atom[i]['a'] =='N' and atom[i]['nbs'] == 4 : 
        atom[i]['a'] = 'N+'
    elif atom[i]['a'] in  ['O','S'] and atom[i]['nbs'] == 1 :
        atom[i]['a'] += '-'

# correcting 'con' definitions
for i in atom :
    atom[i]['con'] = conic(i , atom) 


# the resonating nitrogens in gaunidine, gaunidinium and amidinium
guans = []
for i in atom : 
    if shape(i,atom) == 'C34' : 
        cona = atom[i]['cona']   
        k = got(i,cona,['2N+']) 
        if k  :
            k = k[0]
            for j in cona : 
                if shape(j,atom) == 'N33' : 
                    atom[j]['a'] == 'N+'
                    adds([i,j,k],guans)

# correcting 'con' definitions after guans
for i in atom :
    atom[i]['con'] = conic(i , atom) 

# carbonyls and amides
carbonyl = []
amide = []
ester = []
for i in atom : 
    if shape(i,atom) == 'C34' : 
        cona = atom[i]['cona']     
        if got(i,cona,['2O','2S']) : 
            add(i,carbonyl)
            for j in cona : 
                s = shape(j,atom)
                if s == 'N33' : 
                    add(j,amide)
                elif s == 'O22' and not 'H' in atom[j]['con']: 
                    add(j,ester)
                    

# conjugated systems of double bonds
conj = []
conjs = []
for d in dihedral :
    [i,j,k,l] = dihedral[d]['num']
    bons = [sorted([i,j]) , sorted([j,k]) , sorted([j,l])]
    ords = []
    for b1 in bons :
        for b2 in bond :            
            if b1 == bond[b2]['num'] :
                ords += [bond[b2]['ord']]
                break             
    if ords == [2,1,2] :
        conjs += [dihedral[d]['num']]
        for i in dihedral[d]['num'] :
           if atom[i]['a'] == 'C' : 
               add(i,conj)


defn = {i:'' for i in atom}
def repo(i) : 
    defn[i] += ' ! reparameterization is required'


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

#### triples groups
for i in atom :
    a , cona , con , ncon = attr (i,atom)
    if atom[i]['a'] == 'N' and atom[i]['ncon'] == 1:
        cg(i,'NG1T1')
        if con != 'NC': 
            repo(i) 
    elif atom[i]['a'] == 'C' and atom[i]['ncon'] == 2:
        if got(i,atom[i]['cona'],['3N']) :
            cg(i,'CG1N1')
            if con != 'CCN' or atom[i]['c']: repo(i)                
        elif got(i,atom[i]['cona'],['3C']): 
            if 'H' in atom[i]['con']: 
                cg(i,'CG1T2')  
                if con != 'CCH' or atom[i]['c'] : repo(i) 
            else : 
                cg(i,'CG1T1')
                if con != 'CCC' or atom[i]['c'] : repo(i) 
        elif len (got(i,atom[i]['cona'],['2O','2N','2N+','2S','2P'])) == 2 :
            cg(i,'CG2O7')
            cgd(atom[i]['cona'],{'O':'OG2D5','O-':'OG2D5','N':'NG2D1','S':'SG2D1','S-':'SG2D1'})
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
        for c in benzenes : 
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
                        if s > c : 
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
for c in benzenes : 
    cyc = cycle[c]['num']
    H2 = 0
    for i in cyc:
        j = got(i,atom[i]['cona'],['2O','2S'])
        if j :
            cg(i,'CG2R63')
            H2 = 1
            if atom[j[0]]['a'] == 'O' : 
                adds(cyc,CG2R62)
                cg(j[0],'OG2D4')
    if H2 : 
        H2 = []
        for i in cyc : 
            cgd(atom[i]['cona'],{'H':'HGR62'})

        


# fluorinated Carbons
fluorinated = []
for i in benza : 
    if 'F' in atom[i]['con']: 
        fluorinated += [i]
cgs(fluorinated,'CG2R66')


# amidine Carbons
for i in benza :     
    if atom[i]['a'] == 'C' and atom[i]['con'].count('N') > 1 : 
        cg(i,'CG2R64')
        

# protonated and carbonyl benzenes

for c in benzenes : 
    cyc = cycle[c]['num']
    H3 = 0
    for i in cyc : 
        if atom[i]['a'] == 'N+' : 
            adds(cyc,CG2R62)
            H3 = 1
            break
    if H3 : 
        for i in cyc : 
            cgd(atom[i]['cona'],{'H':'HGR63'})
        
cgd(CG2R62,{'C':'CG2R62'})

# Nitrogens
for i in benza : 
    if not atom[i]['cg']: 
        if 'N' in atom[i]['a']: 
            if atom[i]['ncon'] == 3 : 
                cg(i,'NG2R61')
            else : 
                cons = ''.join([atom[j]['a']+atom[j]['con'] for j in atom[i]['cona']])
                if 'N' in cons or 'O' in cons or 'S' in cons or 'P' in cons : 
                    cg(i,'NG2R62')
                else : 
                    cg(i,'NG2R60')

# the remaining atoms in benza (Carbons and Oxygens)
cgd(benza,{'C':'CG2R61','O':'OG3R60'})


# H on benzenes
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
            cg(i,'NG2R57')
            
### aromatic 7-mem  
for i in aroma7 : 
    if atom[i]['a'] == 'C' : 
        if atom[i]['c'] == 5 : 
            cg(i,'CG2RC7')
        else : 
            cg(i,'CG2R71')        
    

# O , N and S of fura 
cgd(furas,{'O':'OG2R50','S':'SG2R50','N':'NG2R51'})

## N+34 and N23 with three connectivities in 5-mem
for c in cycle5 : 
    cyc = cycle[c]['num']
    for i in cyc : 
        if not atom[i]['cg']:
            if 'N' in atom[i]['a'] : 
                s = shape(i,atom)
                if s in ['N+34','N+33'] : 
                    cg(i,'NG2R52')
                elif s == 'N23': 
                    cg(i,'NG2R50')

### N33 in plana
cgs(plana,'NG2R51')

## amide 5-mem cycles
for i in carbonyl : 
    if atom[i]['c'] == 5 : 
        am = [k for k in atom[i]['cona'] if k in amide]
        cgs(am,'NG2R53')

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
            elif got(i,cona,['2O','2S','2P']) and cyconah :                
                cg(i,'CG2R53')             
cgs(CG2R52,'CG2R52')



# exocyclic double bonds
for c in cycle5 : 
    cyc = cycle[c]['num']
    for i in cyc : 
        cona = atom[i]['cona']
        if atom[i]['a'] == 'C' and not atom[i]['cg'] and atom[i]['ncon'] == 3: 
            j = diff(cona,cyc)                      
            j = got(i,j,['2C','2N','2N+','2O','2S','2P'])
            if not onecyc([i]+j,[3,4,5,6,7]):
                cycona = diff(cona,j)
                if [k for k in cycona if not atom[k]['a'] in ['C','H']] :  
                    cg(i,'CG251O')                
                elif i in conj : 
                    cg(i,'CG25C1')

## the default C=C in 5mem
for c in cycle5: 
    cyc = cycle[c]['num'] 
    for i in cyc : 
        if atom[i]['a'] == 'C' and atom[i]['ncon'] == 3 and not atom[i]['cg'] :  
            cona = atom[i]['cona']
            if got(i,cona,['2C','2N','2N+','2P']) :
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
# guanidine/guanidinium and amidinium

for i in atom : 
    if not atom[i]['cg'] : 
        if shape(i,atom) == 'C34' and atom[i]['con'].count('N') == 3 : 
            cg(i,'CG2N1')
            if not 'N+' in con : 
                repo(i)
        elif atom[i]['a'] == 'C' and i in guans : 
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
    
[cycle[c]['n'] for c in cycle if i in cycle[c]['num'] and cycle[c]['n'] in [3,4,5]]
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
for i in atom : 
    if not atom[i]['cg'] : 
        if atom[i]['a'] == 'C' and atom[i]['c'] in [3,4] and i in carbonyl : 
            cona = atom[i]['cona']
            N33 = [j for j in cona if shape(j,atom) == 'N33']


################################################################################
################################################################################
################################################################################
################################################################################
### nitro, phosphate, sulfate

# nitro
for i in atom : 
    if shape(i,atom) in ['N+34','N35']: 
        cona = atom[i]['cona'] 
        O = [j for j in cona if 'O' in atom[j]['a'] and atom[j]['ncon'] == 1] 
        if len(O) > 1 : 
            cg(i,'NG2O1')
            cgs(O,'OG2N1')

# phosphate
for i in atom : 
    a , cona , con , ncon = attr (i,atom)
    if a == 'P': 
        nO = len([j for j in cona if atom[j]['a'] == 'O-']) 
        q = 5-atom[i]['nbs']-nO
        if q > -1 :             
            cg(i,'PG0')
        elif q == 1 : 
            cg(i,'PG1')
        else : 
            cg(i,'PG2')
        O = [j for j in cona if 'O' in atom[j]['a'] and atom[j]['ncon'] == 1] 
        if 'S' in con : 
            cgs(O,'OG2S1')
        else : 
            cgs(O,'OG2P1')
        N = [j for j in cona if shape(j,atom) == 'N33'] 
        cgs(N,'NG2S3')
    

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
    a , cona , con , ncon = attr (i,atom)
    if a == 'S' : 
        if ncon == 3 : 
            cg(i,'SG3O3')
        elif ncon == 4 : 
            if [j for j in cona if atom[j]['a'] == 'O-'] or atom[i]['nbs'] > 6 : 
                cg(i,'SG3O1')
            else : 
                cg(i,'SG3O2')
                

# -O- in phosphate and sulfate 
for i in atom :
    if atom[i]['cg'] in ['PG0','PG1','PG2','SG3O1','SG3O2','SG3O3']: 
        for j in atom[i]['cona']: 
            if atom[j]['a'] == 'O' and atom[j]['ncon'] == 2 : 
                con = atom[j]['con']
                if con in ['OPP','OPS','OSS']: 
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
        a , cona , con , ncon = attr (i,atom) 
        k = got(i,cona,['2O','2S'])[0]
        nhet = len([j for j in cona if not atom[j]['a'] in 'CH'])
        if nhet == 3 : 
            cg(i,'CG2O6')
            if '-' in con : 
                cg(k,'OG2D2')
        elif '-' in con : 
            cg(i,'CG2O3')
            cg(k,'OG2D2')
        elif con.count('O') + con.count('S') == 2  : 
            cg(i,'CG2O2') 
        elif 'N' in con :
            cg(i,'CG2O1') 
        elif 'H' in con : 
            cg(i,'CG2O4') 
        else : 
            cg(i,'CG2O5') 
            cg(k,'OG2D3')
        if not atom[k]['cg']: 
            atom[k]['cg'] = 'OG2D1'

# amides 
for i in amide : 
    nH = atom[i]['con'].count('H')
    cg(i,'NG2S'+str(nH))


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
                    if [atom[j]['a'] , atom[j]['ncon']] == ['C',4]:
                        if not 'F' in atom[j]['con']:
                            cgd(atom[j]['cona'] , {'H':'HGP5'}) # should it be CH3?


################################################################################
################################################################################
################################################################################
################################################################################
# sulfur   
for i in atom:
    if not atom[i]['cg'] and 'S'in atom[i]['a']: 
        a , cona , con , ncon = attr (i,atom)        
        if ncon == 2 :
            if con.count('S') > 1 or 'O' in con:
                cg(i,'SG301')
            else : 
                cg(i,'SG311')
        elif ncon == 1 :
            if atom[i]['nbs'] == 1 or got(i,cona,'2N+'):                    
                cg(i,'SG302')
            else :
                cg(i,'SG2D1')


################################################################################
################################################################################
################################################################################
################################################################################
### oxygens

for i in atom:
    if not atom[i]['cg'] and 'O' in atom[i]['a'] :
        a , cona , con , ncon = attr (i,atom)
        if ncon == 1 :
            j = cona[0]
            if atom[i]['nbs'] == 1:                
                if atom[j]['a'] == 'S' and atom[j]['nbs'] == 1 : 
                    cg(i,'OG2P1')
                else : 
                    cg(i,'OG312')
            else : 
                if atom[j]['con'] == 'CCCO': 
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
                if atom[i]['c'] == 6 :
                    for j in cona :
                        if atom[j]['nbs'] > atom[j]['ncon']:
                            cg(i,'OG3R60')
                            break
                    if not atom[i]['cg']:
                        cg(i,'OG3C61')
                else :
                    cg(i,'OG301')
        
            
################################################################################
################################################################################
################################################################################
################################################################################
### C=N and C=C

for i in atom : 
    if not atom[i]['cg'] : 
        if atom[i]['a'] == 'N' and atom[i]['ncon'] == 2 : 
            cg(i,'NG2D1')
    
       
for i in atom :
    if not atom[i]['cg'] :
        a , cona , con , ncon = attr (i,atom)
        if [a,ncon] == ['C',3]:
            nH = con.count('H')
            if nH == 2 :
                if i in conj :
                    cg(i,'CG2DC3')
                else :
                    cg(i,'CG2D2')
            else:
                if [j for j in cona if atom[j]['a'] in 'ON+S'] : 
                    cg(i,'CG2D1O')
                elif i in conj : 
                    cg(i,'CG2DC1')
                else :
                    cg(i,'CG2D1')

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

# N in Hydrazine 
cgs([j for j in amine if atom[j]['con'].count('N') > 1],'NG3N1') 


## amine groups
amine = [i for i in amine if not atom[i]['cg']]
for i in amine :
    a , cona , con , ncon = attr (i,atom)
    if [a,ncon] == ['N', 3] and not atom[i]['cg']:
        nH = str(con.count('H'))        
        if con in ['NCCC','NCCH','NCHH','NHHH'] :
            nH = str(con.count('H'))
            cg(i,'NG3'+nH+'1')
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
            a , cona , con , ncon = attr (i,atom)
            nH = con.count('H')            
            # the following priority is taken from paramchem
            if 'N+' in con and nH :
                cg(i,'CG3'+str(nH)+'4') 
            if not atom[i]['cg'] and \
               [j for j in cona if atom[j]['a'] == 'S' and atom[j]['nbs'] == 1] and 'H' in con:
                nH = con.count('H')
                if nH == 2 :
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
            aj , cona , con , ncon = attr (j,atom)
            if a == 'F' and aj != 'C':                
                cg(i,'FGP1')                
            elif j in benza or j in furas or j in aroma7:                
                cgd([i],{'F':'FGR1' , 'L':'CLGR1' , 'B':'BRGR1' , 'I':'IGR1'})
            else :
                nX = con.count('F') + con.count('L') + con.count('B') + con.count('I') 
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
        a , cona , con , ncon = attr (j,atom)
        if [a,ncon] == ['C',3] :
            if got(j,cona,['2N+']): 
                cg(i,'HGR52')
            elif got(j,cona,['2N','2C']):
                nH = con.count('H')
                if nH == 1 :
                    cg(i,'HGA4') 
                elif nH == 2:
                    cgd(cona,{'H':'HGA5'})                                    
            elif j in carbonyl:
                cg(i,'HGR52')
        elif atom[j]['cg'] == 'NG2S3': 
            for k in cona : 
                if '2R6' in atom[k]['cg']: 
                    cg(i,'HGP4')
            if not atom[i]['cg']:
                cg(i,'HGP1')    
                    
        elif a == 'S' :
            cg(i,'HGP3')
        elif a == 'N+': 
            cg(i,'HGP2')      
        elif a != 'C' :
            cg(i,'HGP1')  


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


for i in atom :
    print (i , atom[i]['cg'])


'''

errors = []
if 1 in typs : 
    for i in atom :
        onc = atom[i]['cg']
        pc = typs[i]    
        #print i , onc , pc , onc == pc
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
        if penalty[i] < 15 and match and not atom[i]['a'] == 'H' :
            errors += [str(i) +'\t'+ onc +'\t'+ pc + '\n']
        elif penalty[i] > 15 : 
            errors += [str(i) +'\t'+ onc +'\t'+ pc +'\t'+ str(penalty[i]) + '\n']
            


print '==============='
if errors : 
    for r in errors : 
        print r.strip('\n')

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
