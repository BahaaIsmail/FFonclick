
from os import walk ,getcwd , system
from functions import *
import sys

hessf = sys.argv[1]
#hessf = 'nma.fch'


print ( '\n\n' )
print ( '__________________________________________________________________________________' )
print ( 'Harmonic Parameterization ________________________________________________________' )



print ( '\nThe Hessian file: {} was detected\n'.format(hessf) )



cfg = 0 
path = getcwd()
for root , dirs ,files in walk(path) : 
    for f in files : 
        if f == 'cfg' : 
            cfg = 1 
            break 
   

atoms , bonds = atoms_bonds()
angles = get_angles(atoms)
ureys = get_ureys(angles)
impropers = get_impropers(atoms)
symmetries = get_symmetries(atoms)
nat = len(atoms)

nosymm = 0
if cfg : 
    sf= open('cfg')
    for line in sf : 
        if not '#' in line : 
            if 'nosymmetry' in line : 
                nosymm = 1
                break    
    sf.close()

if nosymm : 
    symmetries = {i:[] for i in atoms}
    



sf = open(hessf)
#____________________________________________________________________________
# freq level of theory to know the scaling factor
sf.readline()
line = sf.readline()
if 'HF' in line :
    scale = 0.80156209 
elif 'DFT' in line :
    scale = 0.927369
elif 'MP2' in line :
    scale = 0.889249
else : 
    scale = 1
scale = scale**2

print ( 'The scaling factor {} will be used'.format(scale) )

#____________________________________________________________________________
# total charge
for line in sf :
    if line.startswith('Charge'):
        totalcharge = int(line.split()[2])
        break

#____________________________________________________________________________
# cartesian coordinates (of the stable conformer which was used to calculate freq)
for line in sf :
    if ('Current cartesian coordinates') in line:
        break
c = []     # the cartesian coordinates used for the harmonic parameters
for line in sf :
    if 'Force Field' in line :
        break
    else :
        line = line.split()
        c += [float(line[i])*0.52917721092 for i in range(len(line))]
cc = {}
for i in atoms:
    j = 3*(i-1)
    cc[i] = c[j:j+3]

for a in angles :
    [i,j,k] = angles[a]['num']
    ci,cj,ck = cc[i], cc[j], cc[k]
    angles[a]['v'] = anglevalue(ci,cj,ck)

for m in impropers :
    [i,j,k,l] = impropers[m]['num']
    ci,cj,ck,cl = cc[i], cc[j], cc[k], cc[l]
    impropers[m]['v'] = dihedralvalue(ci,cj,ck,cl) 
#____________________________________________________________________________
# atomic masses (in case frequencies are needed to be calculated)
for line in sf :
    if ('Real atomic weights') in line:
        break
mass = []
for line in sf :
    if 'Atom fragment info' in line :
        break
    else :
        line = line.split()
        mass += [float(line[i]) for i in range(len(line))]




####=================================================================================================
#        ________________       ______________   
#       /\               \     /\             \  
#       \ \    ___________\    \ \    _________\ 
#        \ \   \__________/     \ \   \________/ 
#         \ \   \_______         \ \   \
#          \ \          \         \ \   \
#           \ \    ______\         \ \   \
#            \ \   \_____/          \ \   \
#             \ \   \                \ \   \__________
#              \ \   \                \ \             \
#               \ \___\                \ \_____________\
#                \/___/                 \/_____________/
####=================================================================================================
#____________________________________________________________________________
## Hessian matrix
sf = open(hessf)
for line in sf : 
    if 'Cartesian Force Constants' in line : 
        break
fcs = []
for line in sf : 
    if 'Dipole Moment' in line : 
        break
    else :
        l = line.split()
        fcs += [float(l[i])*627.503/0.52917721092/0.52917721092 for i in range(len(l))]
sf.close()
nat3 = nat*3
G = [[0 for j in range(nat3)] for i in range(nat3)]
k = 0
for i in range(nat3): 
    for j in range(i+1): 
        G[i][j] = fcs[k]
        G[j][i] = fcs[k]
        k += 1

H = {i:{j:[[0,0,0],[0,0,0],[0,0,0]] for j in atoms} for i in atoms}
for i in H : 
    for j in H :
        r = (i-1)*3
        c = (j-1)*3
        for k in range(3): 
            for l in range(3) : 
                H[i][j][k][l] = -G[r+k][c+l] 


####################################################################################################
## calculating and storing the equilibrium values {b} and the {B} matrices of the bonds
####################################################################################################
# the equilibrium values in Angestroms
b = {i:{j:0 for j in atoms} for i in atoms}
# the B matrices of bonds
B = {i:{j:0 for j in atoms} for i in atoms}
    # (note B matrices are the effective direction of each atom that increases the value of the internal coordinate)
    # B = dq/dx   (q is the internal coordinate and x is a catesian coordinate)

bond_list = [bonds[i]['num'] for i in bonds]
urey_list = [ureys[i]['num'] for i in ureys] 
stick_list = bond_list + urey_list
for stick in stick_list :
    [i,j] = stick  
    xi , xj = cc[i] , cc[j]    
    v = disp(xi,xj)
    n = norm(v)
    b[i][j] = b[j][i] = dist(xi,xj)      # the distance between the two atoms (in ANgstroms)
    B[i][j] = [v[k]/n for k in range(3)] # the effective direction = the B matrix
    B[j][i] = [-k for k in B[i][j]]      # the oposit direction


####################################################################################################
## calculating and storing the eigenvalues and eigenvectors of each submatrix h = H[i][j]  
####################################################################################################
# {L} : the eigenvalues 
L = {i:{j:0 for j in b[i]} for i in b} 
# {V} : the eigenvectors
V = {i:{j:0 for j in b[i]} for i in b}   

print ( '\nThe following atom-atom interactions has been found to be unstable,' )
print ( 'so the associated force constants will not be provided' )
unstable = []
for stick in stick_list : 
    [i,j] = stick
    L[i][j] = eigenval(H[i][j])
    L[j][i] = eigenval(H[j][i])
    if L[i][j] and L[j][i] : 
        V[i][j] = eigenvec(H[i][j],L[i][j]) 
        V[j][i] = eigenvec(H[j][i],L[j][i])
    else:        
        unstable += [stick]
        b[i][j] = 0
        b[j][i] = 0
        if sorted([i,j]) in bond_list : 
            print ( sorted([i,j]) )

bond_list = [stick for stick in bond_list if not stick in unstable]
urey_list = [stick for stick in urey_list if not stick in unstable]
stick_list = bond_list + urey_list
b = {i:{j:b[i][j] for j in b[i] if b[i][j]} for i in b}
B = {i:{j:B[i][j] for j in b[i] if b[i][j]} for i in b}


####################################################################################################
## calculating the force constants of the bonds by both Wilson and projection methods
####################################################################################################
# the force constants of the bonds in kcal mol^(-1) A^(-2)
# conv : converts from kJ/mol/Bohr^2 to kcal/mol/A^2, (division by 2 to be compatible with charmm    

def projection(i,j):  # calculates the force constant of a bond by projection method 
    l = L[i][j]
    u = B[i][j]
    v = V[i][j]
    return sum([l[k]*abs(dot(u,v[k])) for k in [0,1,2]]) 


#@@@@@@@@@@@@@@@@@@@@
print ( '\n\nThe force constant of each individual internal coordinate' )
bondpars = {}
print ( '\n\nBond (atom1, atom2, force constant)' )
for a in bonds:  
    [i,j] = bonds[a]['num']
    f1 , f2 = projection(i,j) , projection(j,i)   
    fp = scale*(f1+f2)/4
    print ( i , '\t' , j , '\t' , round(fp,2) )    
    bondpars[a] = [i , j , fp , b[i][j]]
    #system ('perl hesspl '  + hessf + ' ' + str(i) + ' ' + str(j))
    #system ('perl hesso.pl ' + hessf + ' ' + str(i) + ' ' + str(j))

ureypars = {}
print ()
print ( '\n\nUrey-Bradley (atom1, atom2, force constant)' )
for u in ureys:  
    [i,j] = ureys[u]['num']    
    if L[i][j] and L[j][i] : 
        f1 , f2 = projection(i,j) , projection(j,i)   
        fp = scale*(f1+f2)/4
        print ( i , '\t' , j , '\t' , round(fp,2) )
        ureypars[u] = [i , j , fp , b[i][j]]

####################################################################################################
## calculating the force constants of the angles by both Wilson and projection methods
####################################################################################################
angle_list = [angles[i]['num'] for i in angles] 
def projection(i,j,k):  # calculates the force constant of a bond by projection method 
    U = uv(cross(B[j][k],B[j][i]))    
    li , lk = L[i][j] , L[k][j]    
    bi , bk = cross(B[i][j],U) , cross(U,B[k][j])
    f = factor(i,j,k)
    
    vi , vk = V[i][j] , V[k][j]
    fi = f/(b[i][j]**2)/sum([li[n]*abs(dot(bi,vi[n])) for n in [0,1,2]]) 
    fk = f/(b[k][j]**2)/sum([lk[n]*abs(dot(bk,vk[n])) for n in [0,1,2]])
    f1 = 1/(fi+fk)

    
    vi , vk = V[j][i] , V[j][k]
    fi = f/(b[i][j]**2)/sum([li[n]*abs(dot(bi,vi[n])) for n in [0,1,2]])    
    fk = f/(b[k][j]**2)/sum([lk[n]*abs(dot(bk,vk[n])) for n in [0,1,2]])    
    f2 = 1/(fi+fk)
    
    f = (f1+f2)/4
    #print ( f1  # print (ing f1 is verified by hess2ff.pl
    return f

def factor(i,j,k) : 
    # j is th central atom
    # i is the moving atom
    # k will be replaced with the terminal atoms in the sharing angles 
    U = uv(cross(B[j][k],B[j][i]))    
    b1 = cross(B[i][j],U) 
    directions = []    
    for angle in angles : 
        [a,c,b] = angles[angle]['num']
        if j == c and intercept([i,k],[a,b]) == [k]  : 
            k = difference([a,b,c],[i,j,k])[0]
            U = uv(cross(B[j][k],B[j][i]))    
            b2 = cross(B[i][j],U)
            directions += [b2]
    n = len(directions)
    if n == 0 : 
        return 1
    fact = 0
    for d in directions : 
        fact += abs(dot(b1,d))
    fact /= n
    return 1+fact


anglepars = {}
print ()
print ( '\n\nAngles (atom1   atom2   atom3   force constant)' )
for a in angles :
    [i,j,k] = angles[a]['num']
    fp = scale * projection(i,j,k)
    print ( i , '\t' , j , '\t' , k , '\t' , round(fp,2) )
    anglepars[a] = [i , j , k , fp , angles[a]['v']]


####################################################################################################
## calculating the force constants of the impropers by both Wilson and projection methods
####################################################################################################
def domain(t)   :  return max(-1,min(1,t))
improper_list = [impropers[i]['num'] for i in impropers] 
def projection(i,j,k,l):  # calculates the force constant of a bond by projection method
    kl = uv(disp(cc[k],cc[l]))
    ji = uv(disp(cc[j],cc[i]))
    jk = uv(disp(cc[j],cc[k]))
    jk = va(jk,dot(ji,jk))
    alt = disp(jk,ji)
    u = uv(alt)
    n = uv(cross(jk,kl))
    cost = domain(dot(u,n))
    sint = sqrt(1-cost**2)
    proj = va(alt,sint)

    n = [-a for a in n]
    Lj , Lk , Ll = L[j][i] , L[k][i] , L[l][i]
    vj , vk , vl = V[j][i] , V[k][i] , V[l][i]
    fj = sum([Lj[a]*abs(dot(n,vj[a])) for a in [0,1,2]])
    fk = sum([Lk[a]*abs(dot(n,vk[a])) for a in [0,1,2]])
    fl = sum([Ll[a]*abs(dot(n,vl[a])) for a in [0,1,2]])

    alt = norm(alt)
    f = alt**2*(fj + fk + fl)
    return f

#@@@@@@@@@@@@@@@@@@@@
improperpars = {}
print ( '\n\nImpropers (atom1   atom2   atom3   atom4   force constant)' )
for m in impropers :
    [i,j,k,l] = impropers[m]['num']
    f = scale * projection(i,j,k,l)
    print ( i , '\t' , j , '\t' , k , '\t' , l , '\t', round(f,2) )
    improperpars[m] = [i,j,k,l,f,0]
#===================================================================================================

### equivalencing the parameters
def equivalents(ics,icpars) :
    equics = equicoord (ics,symmetries)
    pool = []
    for b in icpars :
        if not b in pool :
            if len(equics[b]) > 1 : 
                print ( [ics[q]['num'] for q in equics[b]] )
            n = len(equics[b])
            avf = sum([icpars[q][-2] for q in equics[b]])/n
            avd = sum([icpars[q][-1] for q in equics[b]])/n
            pool += equics[b]
            for q in equics[b] :
                icpars[q][-2] = avf
                icpars[q][-1] = avd
    return icpars
print ( '\nThe program equivalencing scheme has detected the following equivalent bonds' )
bondpars = equivalents(bonds,bondpars)
print ( '\nThe program equivalencing scheme has detected the following equivalent Urey-Bradley interactions' )
ureypars = equivalents(ureys,ureypars)
print ( '\nThe program equivalencing scheme has detected the following equivalent angles' )
anglepars = equivalents(angles,anglepars)
print ( '\nThe program equivalencing scheme has detected the following equivalent improper angles' )
improperpars = equivalents(impropers,improperpars)
print ( '\nThe equivalences will be considered in the output parameter file' )

### writing the parameters in a temporary output file 
tf = open('tempars','a')

for b in bondpars :
    [i,j,f,v] = bondpars[b]
    tf.write('bond%3s%3s%10s%8s\n'%(str(i),str(j),str(round(f,3)),str(round(v,3))))
    
for u in ureypars :
    [i,j,f,v] = ureypars[u]
    tf.write('urey%3s%3s%10s%8s\n'%(str(i),str(j),str(round(f,3)),str(round(v,3))))

for a in anglepars :
    [i,j,k,f,v] = anglepars[a]
    tf.write('angle%3s%3s%3s%10s%8s\n'%(str(i),str(j),str(k),str(round(f,3)),str(round(v,2))))
    
for m in improperpars :
    [i,j,k,l,f,d] = improperpars[m]
    tf.write('improper%3s%3s%3s%3s%10s   0   0.00\n'%(str(i),str(j),str(k),str(l),str(round(f,3))))
    
tf.close()


print () 
print ( 'Harmonic parameterization is done  _______________________________________________' )
print ( '__________________________________________________________________________________' )


tf = open('hessdone','w')
tf.close()
