
# linear algebra module
# in this script all the functions work on lists not dictionaries


#       ________               ____________
#      /\__   __\             /\   ______  \
#      \/_/\  \_/             \ \  \____/\  \
#         \ \  \               \ \  \   \ \  \  
#          \ \  \               \ \  \___\_\  \
#           \ \  \               \ \   _____\  \
#            \ \  \        __     \ \  \____/\  \
#             \ \  \      /\ \     \ \  \   \ \  \
#              \_\  \_____\_\ \     \ \  \   \ \  \
#             /\_______________\     \ \__\   \ \__\
#             \/_______________/      \/__/    \/__/
#

from math import sqrt , acos , pi , cos , atan2 , copysign
from copy import deepcopy

####################################################################################################
## The required Linear Algebra methods
####################################################################################################

#___________________________________________
### Vector algebra
def disp(v1,v2):  # the displacement vector from atom1 to atom2 (v2-v1) 
    return [v2[i]-v1[i] for i in range(len(v1))]

def dot(v1,v2) :   # the dot product of two vectors
    return sum([v1[i]*v2[i] for i in range(len(v1))])


def norm(v) :  # = the magnitude of a vector   
    return sqrt(dot(v,v))

def dist(v1,v2):  # the distance between two atoms = the norm of the corresponding displacement vector, this gives the bond length
    # v1 and v2 are the cartesian coordinates of atom1 and atom2 
    v = disp(v1,v2)
    return norm(v)

def lstd(d) :  # convert a dictionary to a list
    return [d[i] for i in d]

def dicl(l) :  # converts a list to a dictionary
    return {i+1:l[i] for i in range(len(l))}

def va(v,a) : # multiplies a vector by a scaler
    return [i*a for i in v]


def uv(v) : # calculates the unit vector of a vector
    n = norm(v)
    return [i/n for i in v]

def cross (v1,v2) :  # v1 x v2
    return [v1[1]*v2[2]-v1[2]*v2[1] , v1[2]*v2[0]-v1[0]*v2[2] , v1[0]*v2[1]-v1[1]*v2[0]]


def outer (v,w) :  # returns the outer products of two vectors
    n = len(v)
    return [[v[i]*w[j] for j in range(n)] for i in range(n)]

#___________________________________________
### Matrix algebra
def zeros(m,n) :  # generates a matrix of zeros with size m*n
    return [[0 for j in range(n)] for i in range(m)]

def eye(n) : # return an identity matrix of size n
    A = zeros(n,n)
    for i in range(n) :
        A[i][i] = 1
    return A


def T(M):  # the transpose of a matrix = switching between the rows and columns 
    return [[M[i][j] for i in range(len(M))] for j in range(len(M[0]))]

def MM(m1,m2) :  # multiplying two matrices m1*m2
    m2 = T(m2)
    m , n = len(m1) , len(m2)
    M = [[0 for j in range(n)] for i in range(m)]    
    for i in range(m): 
        for j in range(n): 
            M[i][j] = dot(m1[i],m2[j])
    return M

def M_M(m1,m2) :  # subtracts m1 - m2 
    return [[m1[i][j]-m2[i][j] for j in range(len(m1[0]))]for i in range(len(m1))]

def Ma(m,a):   # multiplies a matrix with a constant
    return [[m[i][j]*a for j in range(len(m[0]))]for i in range(len(m))]

def Mv(m,v) :  # multiplies the matrix m by the vector v (m*v)
    return [dot(m[i],v) for i in range(len(m))]

def block(M,i1,i2, j1,j2):  # returns a submatrix from a matrix
    return [ [M[i][j] for j in range(j1,j2+1)] for i in range(i1,i2+1) ]    

def column (M,j):  # returns a column from a matrix
    return [M[i][j] for i in range(len(M))]

def pinv(v) :  # teh pseudo inverse of an (n*1) matrix (a special function)
    n = norm(v) 
    return [i/n**2 for i in v]

def trace(M) :  # returns the sum of the diagonal elements of a matrix
    return M[0][0] + M[1][1] + M[2][2]

def det(M) :  # the determinant of a matrix
    dt =   M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) \
         - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) \
         + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])
    return dt


def minors(M) :   # the minors of a matrix is the sum of the determinants of its submatrices 
    mnr =   M[1][1] * M[2][2] - M[1][2] * M[2][1] \
          + M[0][0] * M[2][2] - M[0][2] * M[2][0] \
          + M[0][0] * M[1][1] - M[0][1] * M[1][0]
    return mnr



#___________________________________________
### Numerical methods

def eigenval(M):
    a = - trace(M)
    b =   minors(M)
    c = - det(M)
    
    Q = (a**2-3*b) / 9
    R = (2*a**3-9*a*b+27*c) / 54
    f = (R**2-Q**3)
    if Q < 0.1 : 
        Q = 0 

    eps = 1e-8
    if Q >= 0 and f < eps :                
        r = R/sqrt(Q**3)
        r = max(-1,min(r,1))
        t = acos(r)/3
        pi2 = 2*pi/3
        com = - 2*sqrt(Q)
        a3 = a/3
        r1 = (com * cos(t)) - a3   
        r2 = (com * cos(t+pi2)) - a3  
        r3 = (com * cos(t-pi2)) - a3
        r = [r1 , r2 , r3]
        for i in range(3) : 
            if r[i] < 0 and r[i] > -0.1 : 
                r[i] = 0
        if not negative(r) : 
            return r

def negative (L) :
    c = 0
    for i in L :
        if i < -1e-3 :
            return 1 
        elif type(i) is complex : 
            c += 1
            if c > 1 : 
                return 1

def eigenvec (M , L) :
    v = []
    for i in L :
        A = deepcopy(M)   
        for j in range(3) :
            A[j][j] = A[j][j] - i 

        a  = (A[0][0]*A[1][1]) - (A[0][1]*A[1][0]) 
        b  = (A[0][1]*A[1][2]) - (A[0][2]*A[1][1])         
        c  = (A[0][2]*A[1][0]) - (A[0][0]*A[1][2])     
               
        vi = [b/a , c/a , 1.0]
        n = norm(vi)
        vi = [vi[j]/n for j in [0,1,2]]
        v += [vi]
        
    return v

'''
def QR(A) :
    m = len(A) 
    n = len(A[0])
    # initialize Q
    Q = eye(n)
    # loopin over the columns to be zeroed
    for j in range(n-1) :
        V = [A[i][j] for i in range(j,m)]
        normv = norm(V)
        if V[0] < 0 :
            V[0] -= normv
        elif V[0] > 0 :
            V[0] += normv   
    
        H  = outer(V,V)
        v  = dot(V,V)
        z = len(V)
        I = eye(z)


        H = [[I[i][j] - 2*H[i][j]/v for j in range(z)] for i in range(z)]
        
        if z < m :
            P = eye(m)
            k = m-z
            for i in range(z) :
                for j in range(z) :
                    P[i+k][j+k] = H[i][j]
        elif z == 3 :
            P = H
      
        A = MM(P,A)
        Q = MM(Q,P)

    R = A + []     
    return Q , R 


def eigenval(M) :
    A = deepcopy(M)
    for i in range(100) :
        Q , R = QR(A)
        A = MM(R,Q)
        L = [A[0][0] , A[1][1] , A[2][2]]
    return L
'''
#### least squres by QR giving the same results of pinv in matlab
## given A(m*n) and b(m)
## m >= n
## returns the solution x
def LSQR(A,b): 
    # A is not neccessarily to be symmetric or even square 
    m = len(A)
    n = len(A[0])
    r = m-n
    # the reduced QR factorization
    U = eye(m)
    S = deepcopy(A)
    for k in range(n):
        w = column(S,k)[k:m+1]
        w[0] += copysign(norm(w),w[0])    
        w = va(w,1/norm(w))
        o = M_M(eye(m-k),Ma(outer(w,w),2))
        s = MM(o,block(S,k,m-1,k,n-1))
        u = MM(block(U,0,m-1,k,m-1),o)
        for i in range(k,m):
            h = i-k
            for j in range(k,n):
                S[i][j] = s[h][j-k]
        for i in range(m) :
            for j in range(k,m):
                U[i][j] = u[i][j-k]                
    # Gaussian Elimination 
    x = deepcopy(b)
    x = Mv(T(U),x)
    for i in range(m):
        S[i] += [x[i]]    
    for i in range(n) :
        # make the diagonals = 1
        S[i]  = [S[i][j]/S[i][i] for j in range(n+1)]        
        for j in range(i):
            # zero out all elments above the diagonal element i
            S[j]  = [S[j][k]-S[j][i]*S[i][k] for k in range(n+1)]
    x = column(S,n)[0:n]
    return x


#___________________________________________
### geometry methods
def anglevalue(i,j,k):  # returns the value of the angle mad by three atoms given the cartesian coordinates of the three atoms
    ij , kj= uv(disp(i,j)) , uv(disp(k,j))
    cosa = dot(ij,kj)
    cosa = max(min(cosa,1),-1)   
    return acos(cosa)*180/pi


def dihedralvalue(i,j,k,l):  # returns the value of the dihedral mad by four atoms given the cartesian coordinates of the four atoms
        vij , vjk , vkl = uv(disp(j,i)) , uv(disp(k,j)) , uv(disp(l,k))
        ni = cross(vij,vjk)
        nk = cross(vjk,vkl)
        m  = cross(ni,vjk)
        x , y = dot(ni,nk) , dot(m ,nk)
        return atan2(y,x)*180/pi














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

def intercept(L1,L2) :
    return [i for i in L1 if i in L2]

def difference(L1,L2) :
    return [i for i in L1 if not i in L2]

def env1(i,atoms) :
    env  = [atoms[i]['a']] + sorted([atoms[j]['a'] for j in atoms[i]['cona']])
    return ''.join(env)


def atoms_bonds() :
    sf = open('tmp')
    line = sf.readline().split()
    atoms = {i+1:{'a':line[i] , 'ord':{} , 'cona':[] , 'ncon':0} for i in range(len(line))}
    bonds = {}
    k = 1
    for line in sf :
        line = line.split()
        if len(line) == 3 :
            bonds[k] = {}
            bonds[k]['num'] = sorted([int(line[0]) , int(line[1])])
            bo = line[2]
            if bo in ['ar','Ar','AR','co2','CO2','Co2','1.5'] :
                bo = 1.1
            elif bo == '2.5' :
                bo = 2.1
            else :
                bo = int(float(bo))
            bonds[k]['ord'] = bo
            [i,j] = bonds[k]['num']
            atoms[i]['ord'][j] = bo
            atoms[j]['ord'][i] = bo
            atoms[i]['cona'] += [j]
            atoms[j]['cona'] += [i]
            k += 1
    for i in atoms :
        atoms[i]['cona'] = sorted(atoms[i]['cona'])
        atoms[i]['ncon'] = len(atoms[i]['cona'])
        atoms[i]['con'] = env1(i,atoms)
    sf.close()
    return atoms , bonds

def get_angles(atoms):
    angles = {}
    a = 1
    for i in atoms :
        for j in atoms[i]['cona'] :
            for k in atoms[i]['cona'] :
                if k > j :
                    angles[a] = {}
                    angles[a]['num'] = [j , i , k]
                    a += 1
    return angles

def get_ureys(angles) :
    return {u:{'num':[angles[u]['num'][0],angles[u]['num'][2]]} for u in angles}
    

def get_dihedrals(atoms,bonds):
    dihedrals = {}
    d = 1
    for b in bonds :
        [i,j] = bonds[b]['num']
        ks , ls  = atoms[i]['cona'] , atoms[j]['cona']
        for k in ks :
            if k != j : 
                for l in ls :
                    if not l in [i,k] : 
                        dihedrals[d] = {}
                        dihedrals[d]['num'] = [k , i , j , l]
                        d += 1
    return dihedrals

def get_impropers(atoms) :
    impropers = {}
    m = 1
    for i in atoms :
        k = atoms[i]['cona']
        if atoms[i]['a'] == 'C' and atoms[i]['ncon'] == 3 :            
            for j in range (3) :
                if atoms[k[j]]['a'] in ['O','S','N','N+'] and atoms[i]['ord'][k[j]] == 2 :
                    impropers[m] = {}
                    impropers[m]['num'] = [i] + k
                elif atoms[k[j]]['a'] == 'C' and ('O' in atoms[k[j]]['con']  or 'N' in atoms[k[j]]['con']  ) : 
                    impropers[m] = {}
                    impropers[m]['num'] = [i] + k
        elif atoms[i]['a'] == 'N+' and atoms[i]['ncon'] == 3 and atoms[i]['con'].count('O') == 2: 
            impropers[m] = {}
            impropers[m]['num'] = [i] + k
    return impropers

def get_cycles(atoms,bonds,angles,dihedrals) :

    cycles = {}
    c = 1
    def addcycle(c,cys) :
        for s in cycles :
            cyc = cycles[s]['num']
            if intercept(cyc,cys) == cyc :
                return c
        cycles[c] = {}
        cycles[c]['num'] = cys
        cycles[c]['n'] = len(cys)
        
        n = cycles[c]['n']
        cyc = cycles[c]['num']
        cycles[c]['bon'] = []
        cycles[c]['ord'] = []
        for j in range(n-1) :
            bon = sorted([cyc[j],cyc[j+1]])
            for k in bonds :
                if bon == bonds[k]['num'] :
                    cycles[c]['bon'] += [k]
                    cycles[c]['ord'] += [bonds[k]['ord']] 
        return c + 1        


    # 3 mem cycles
    for a in angles :
        [i,j,k] = angles[a]['num']
        if i == k :
            cycles[c] = [i,j,k]
            c += 1

    # 4 mem cycles
    for d in dihedrals :
        [i,j,k,l] = dihedrals[d]['num']
        if i in atoms[l]['cona']:
            c = addcycle(c,[i,j,k,l])

    # 5-mem cycles :
    for d in dihedrals :
        da = dihedrals[d]['num']
        i , j = da[0] , da[3]
        ic1 = atoms[i]['cona']
        ic2 = atoms[j]['cona'] 
        fifth = intercept(ic1,ic2)
        fifth = difference(fifth , d)
        for k in fifth :
            cyc = da + [k]
            c = addcycle(cyc)


    # 6-mem and 7-mem cycles
    def cyc67 (ics1 , ics2):
        cycs = []
        for d in ics1:
            da = dihedrals[d]['num']  
            d1 = atoms[da[0]]['cona']
            d2 = atoms[da[3]]['cona']
            for i in ics2 :
                ia = ics2[i]['num']
                if not intercept(da,ia):
                    b = [ia[0] , ia[-1]]
                    j = intercept(d1,b)
                    k =  intercept(d2,b)
                    if j and k and j[0] != k[0] :
                        cyc = da + [k[0],j[0]]
                        if len(ia) == 3 :
                            cyc = da + [k[0],ia[1],j[0]]
                        cycs += [cyc]
        return cycs 
                   
    cycs67 = cyc67(dihedrals,bonds) + cyc67(dihedrals,angles)
    for cyc in cycs67 :
        c = addcycle(c,cyc)

    return cycles

def get_symmetries(atoms) :
    vascular = {}
    for a in atoms :
        pool = [a]+[i for i in atoms[a]['cona']]
        offs = [[[a]]]+[[[i for i in atoms[a]['cona']]]]    
        children = [i for i in atoms[a]['cona']]    
        while children :
            newoff = [[j for j in atoms[i]['cona'] if not j in pool] for i in children]
            if newoff :            
                offs += [newoff]       
            children = [j for i in newoff for j in i]                
            pool += children
        vascular[a] = [sorted([sorted([atoms[k]['con'] for k in j]) for j in i]) for i in offs]

    uniqs = {}   # unique atoms as keys and the symmetric atoms as definitions 
    pool = []
    for i in atoms:
        if not i in pool :
            temp = [j for j in atoms if j > i and vascular[i] == vascular[j]]
            uniqs[i] = temp
            pool += temp
    return uniqs   


def equicoord(coords,uniqs) :
    def equi(ic1,ic2) :
        n = len(ic1)
        k = 0
        uniqslist = [[i]+ uniqs[i] for i in uniqs]
        for i in range(n) : 
            for u in uniqslist : 
                if ic1[i] in u and ic2[i] in u : 
                    k += 1 
                    break
                elif ic1[i] in u and ic2[i-n] in u :    # to account for the reverse order
                    k += 1 
                    break
        if k == n : 
            return 1

    equics = {c:[] for c in coords}
    tmp = [c for c in coords]
    pool = []
    for i in tmp : 
        if i in coords : 
            equics[i] += [i]
            if not i in pool :
                for j in tmp : 
                    if j > i and not j in pool: 
                        if equi (coords[i]['num'],coords[j]['num']) :
                            equics[i] += [j]
                            pool += [j]
    return equics



def conic (i,atom) :
    # turns a list of atoms in numbers into a string of their labels
    # sometimes i need to coine a string to an array of atoms so tha
    # input: 2,[1,4,5,3]
    # output: 'CCNCO'
    ic = atom[i]['cona']
    conas = [atom[i]['a']] + sorted([atom[j]['a'] for j in ic ])
    return ''.join(conas)


def strcona (L,atom): 
    # turns a list of atoms in numbers into a string of their labels
    return ''.join(sorted([atom[j]['a'] for j in L]))


def shape(i,atom) : 
    return atom[i]['a']+str(atom[i]['ncon'])+str(int(atom[i]['nbs']))

def attr(i,atom) :
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


def got (i,L,x,atom) :
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

   
def hetero(L,atom) : 
    hets = []
    for i in L : 
        if not atom[i]['a'] in ['C','H'] :
            hets += [i]
    return hets
        



def onecyc(L,N,cycle): 
    # L is the list of atoms tested wether they are in one cycle with size n in the list n
    Ls = sorted(L)
    for c in cycle: 
        if cycle[c]['n'] in N : 
            if intercept(Ls,cycle[c]['num']) == Ls :
                return 1 

