import sys
import pyscsi


pyscsi.gv.globval.H_exact    = False; pyscsi.gv.globval.quad_fringe = False;
pyscsi.gv.globval.Cavity_on  = False; pyscsi.gv.globval.radiation   = False;
pyscsi.gv.globval.emittance  = False; pyscsi.gv.globval.IBS         = False;
pyscsi.gv.globval.pathlength = False; pyscsi.gv.globval.bpm         = 0;

pyscsi.Read_Lattice('/home/bengtsson/projects/in/lattice/commissioning')

pyscsi.Ring_GetTwiss(True, 0.0); pyscsi.printglob()

if False:
    import ctypes
    TotalTune = 2*ctypes.c_double
    TotalTune = TotalTune.from_address(int(pyscsi.gv.globval.TotalTune))
    print
    print TotalTune[0], TotalTune[1]

# Print list of global variables.
print
print pyscsi.gv

# Print list of attributs for globval.
print dir(pyscsi.gv.globval)

print
print pyscsi.gv.globval.H_exact
print pyscsi.gv.globval.TotalTune
print pyscsi.gv.globval.TotalTune

sys.stdout.write('\n')
sys.stdout.write('%6.3f %6.3f\n' %
                 (pyscsi.gv.globval.gvec('TotalTune')[0],
                  pyscsi.gv.globval.gvec('TotalTune')[1]))
sys.stdout.write('%6.3f %6.3f\n' %
                 (pyscsi.gv.globval.gvec('Chrom')[0],
                  pyscsi.gv.globval.gvec('Chrom')[1],
                  ))

sys.stdout.write('\n')
for k in range(0, 6):
    sys.stdout.write('%6.3f' % pyscsi.gv.globval.gvec('CODvect')[k])
sys.stdout.write('\n')

sys.stdout.write('\n')
for j in range(0, 6):
    for k in range(0, 6):
        sys.stdout.write('%10.3e' % 
                         pyscsi.gv.globval.gvec('OneTurnMat')[6*j+k])
    sys.stdout.write('\n')

#sys.stdout.write('\n')
#sys.stdout.write('%10.3e' % pyscsi.gv.globval.gvec('OneTurnMat')[0][0])

print
print pyscsi.gv.Cell
print pyscsi.gv.Cell[0]
print pyscsi.gv.Cell[1]
#print pyscsi.gv.Cell[0].A[0]
