

from os import walk , getcwd , system , remove
import sys
from functions import *

tf = open('tempars','w')
tf.close()

#_____________________________________________________
### exploring the input files 
hessf , espf , mol2 , scnfs , cfg = 0 , 0 , 0 , [] , 0
path = getcwd()
for root , dirs , files in walk(path) :
    for f in files :
        if '.fch' in f or '.fchk' in f :
            hessf = f
        elif '.esp' in f :
            espf = f
        elif '.scn' in f :
            scnfs += [f]
        elif '.mol2' in f :
            mol2 = f
        elif 'cfg' in f :
            cfg = f

alternative = 0
if cfg :
    sf = open(cfg)
    for line in sf :
        if line.startswith('alternative') :
            print (line)
            line = line.split()
            if len(line) > 1  :
                alternative = 1
    sf.close()


if not (hessf or espf or scnfs or mol2 or alternative) : 
    print ('\n\n')
    print ('  _____________________________________________________________________________')
    print ('      No input files (.fchk , .esp , .scn , .mol2 , .cfg) have been detected   ')
    print ('                           Program termination                                 '           )
    print ('  _____________________________________________________________________________')
    print ('\n\n')
    sys.exit()
	

#_____________________________________________________
### extracting the atom types and connectivities 
if hessf :
    tf = open('tmp','w')
    #_______________
    # atom symbols
    d = {1:'H' , 2:'M' ,  6:'C' ,  7:'N' ,  8:'O' ,  9:'F' , 10:'E' ,
                13:'A' , 15:'P' , 16:'a' , 17:'L' , 35:'B' , 53:'I' }
    sf = open(hessf)
    for line in sf :
        if 'Atomic numbers' in line :
            break
    nat = 0
    for line in sf :
        if 'Nuclear charges' in  line :
            break
        else :
            line = line.split()
            nat += len(line) 
            for s in line :            
                tf.write(d[int(s)] + ' ')
    tf.write('\n')

    #____________________________________________________________________________
    # bonds = {} : the bound atoms to each atom
    # maxb : the maximum number of bonds found for an atom
    # the atoms connected to each atom are extracted as atom[i]['con'] = [...]
    # bond orderlist is extracted alongside the bondlist
    atoms = [i for i in range(1,nat+1)]
    for line in sf : 
        if 'MxBond' in line : 
            line = line.split()
            maxb = int(line[2])
            break

    for line in sf : 
        if 'IBond' in line : 
            break 
    # bonds (connectivities)
    L = []
    for line in sf : 
        if 'RBond' in line : 
            break 
        else : 
            line = line.split()
            L += [int(i) for i in line]

    # bond orders  (may be needed for fixing the aromatic bonds)
    R = []
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
            if s[k] > i:
                tf.write(str(i)+ ' ' + str(s[k]) +  ' ' + str(t[k]) + '\n')
    sf.close()
    tf.close()    
else :
    print ('\n')
    print ('  _______________________________________________________________________________')
    print ('      WARNING::: The hessian has not been detected (.fchk or .fch is required)   ')
    print ('                 No harmonic force constants will be provided                    '           )
    print ('  _______________________________________________________________________________')
    print ('\n')


if not hessf and mol2 :
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


if scnfs and not (hessf or mol2): 
    print ('\n')
    print ('  _______________________________________________________________________________')
    print ('      WARNING::: Bond orders have not been detected (.fchk or .mol2 is required)  ')
    print ('                 No Lenard-Jones parameters will be provided unless              '   )
    print ('                 a MOL2 file is detected                                         ')
    print ('  _______________________________________________________________________________')
    print ('\n')
    scnf = scnfs[0]
    tf = open('tmp','w')
    sf = open(scnf)
    for line in sf : 
            if 'Symbolic Z-Matrix:' in line : 
                    break 
    for line in sf : 
            if 'Variables:' in line : 
                    break 
            else :
                    line = line.split()
                    tf.write(line[0] + ' ') 
    tf.write('\n')	
    
    for line in sf : 
            if '! A' in line : 
                    break 
            elif '! R' in line : 
                    b = line.split()[2].strip('R(').strip(')').split(',')
                    [i,j] = [b[0] , b[1]]
                    tf.write(i + ' ' + j + ' 0 \n')
    sf.close()
    tf.close()
elif not (scnfs or alternative) : 
    print ('\n')
    print ('  ____________________________________________________________________________________')
    print ('      WARNING::: PES scans have not been detected                                     ')
    print ('                 (.scn file is required or the alternative way has to be claimed)     ')
    print ('                 No torsional parameters will be provided                             '           )
    print ('  ____________________________________________________________________________________')
    print ('\n')
	

    
if espf and not (hessf or mol2 or scnfs):
    print ('\n')
    print ('  __________________________________________________________________')
    print ('      WARNING::: Connectivities are not introduced,                 ' )
    print ('                 (.fchk , .mol2 , or .scn is required)              ')
    print ('                 the equivalencing scheme is disabled               ')
    print ('                 proper charge fitting is not guranteed             ')
    print ('                 unless the charge options in the configuration     ')
    print ('                 file (.cfg) are carfully introduced                ')
    print ('  __________________________________________________________________')
    print ('\n')
    tf = open('tmp','w')
    sf = open(espf)
    for line in sf :
        if 'Charge =' in line :
            break
    for line in sf :
        line = line.split()
        if len(line) == 4 :
            tf.write(line[0] + ' ' )
        else :
            break
    sf.close()
    tf.close()
elif not espf : 
    print ('\n')
    print ('  ____________________________________________________')
    print ('      WARNING::: ESP file has not been detected       ')
    print ('                 No charges will be provided          '           )
    print ('  ____________________________________________________')
    print ('\n')

	


#__________________________________________________________
### Running the codes in succession

if hessf or espf or mol2 or scnfs : 
    print ('\n\n')
    print ('  _______________________________________________________________________________')
    print ('                        Onclickff program has started                            ')
    print ('      For detailed demonstration of the valid input files and entries            ')
    print ('      please check the {Readme} file                                             ')
    print ('                                                                                 ')
    print ('      For detailed information about the parameterization steps and the          ')
    print ('      decisions made by the program check the output log file                    ')
if espf or scnfs : 
    print ('\n               The entries in the configuration file are considered            ')
    print ('      Careful coordination between the entries in the configuration file         ')
    print ('      is strongly recommended,  otherwise the program will resolve               ')
    print ('      any contradictions his way and provide information and warning about that  ')
    print ('      whenever possible                                                          ')
if hessf or espf or mol2 or scnfs : 
    print ('  _______________________________________________________________________________')
    print ('\n\n')



def copyin(txt) :
    tf.write('\n\n')
    sf = open(txt)
    for line in sf :
        tf.write(line)
    sf.close()
    tf.write('\n\n')

tf = open('log.txt','w')    

#_____________
if hessf : 
    system('python hess.py ' + hessf + ' > harmonics')
    copyin('harmonics')

'''   
#_____________
if hessf or mol2 : 
    system('python typer.py > LJ')
    copyin('LJ')

'''
lj = '''
lj  -0.078 	2.05 	-0.01 	1.9   
lj  -0.024 	1.34 	-0.024 	1.34  
lj  -0.024 	1.34 	-0.024 	1.34  
lj  -0.024 	1.34 	-0.024 	1.34  
lj  -0.11 	2.0  	-0.11 	2.0   
lj  -0.12 	1.7  	-0.12 	1.4   
lj  -0.2 	1.85 	-0.2 	1.55  
lj  -0.046 	0.2245 	-0.046 	0.2245
lj  -0.078 	2.05 	-0.01 	1.9 
lj  -0.024 	1.34 	-0.024 	1.34
lj  -0.024 	1.34 	-0.024 	1.34
lj  -0.024 	1.34 	-0.024 	1.34
'''
tf = open('tempars','a')
tf.write(lj)
tf.close()
tf.close()

if espf or scnfs : 
    repeat = 1
    while repeat :
        tf = open('log.txt','a')
        # taking a copy of the configuration file
        tf.write ('\n\n_____________________________________________________________________\n')
        tf.write('A copy of the configuration file for the run number  ({})\n'.format(repeat))
        sf = open('cfg')
        marks = ['fc' , 'gc' , 'ec' , 'kc' , 'fp' , 'ft' , 'repeat' , 'phases' , 'kmax' , 'maxenergy' , 'nsteps' , 'weight']
        for line in sf :
            if not line.startswith('#') :                
                check = line.split()
                if intercept(check,marks) : 
                    tf.write(line)           
        sf.close()
        tf.write('_____________________________________________________________________\n\n')

        #_____________
        if espf :
            system('python resp.py ' + espf + ' > charges')
            copyin('charges')

        #_____________
        if scnfs or alternative :
            scans = ''.join([f+' ' for f in scnfs])
            system('python mcsa.py ' + scans + '   ' + str(alternative) + ' > torsionals')
            copyin('torsionals')
        tf.close()

        ask = input('Run charge and/or torsional parameterzation again? (y or n)   ')
        if not ask.replace(' ','').replace('\t','') == 'y' :
            repeat = 0
        else :
            repeat += 1
            
                
tf.close()            


    

remove('tmp')
remove('tempars')
remove('harmonics')
remove('charges')
remove('torsionals')

