import pyscsi


pyscsi.gv.globval.H_exact    = False; pyscsi.gv.globval.quad_fringe = False;
pyscsi.gv.globval.Cavity_on  = False; pyscsi.gv.globval.radiation   = False;
pyscsi.gv.globval.emittance  = False; pyscsi.gv.globval.IBS         = False;
pyscsi.gv.globval.pathlength = False; pyscsi.gv.globval.bpm         = 0;

pyscsi.Read_Lattice('/home/bengtsson/projects/in/lattice/commissioning')

pyscsi.Ring_GetTwiss(True, 0.0); pyscsi.printglob()

import ctypes
TotalTune = 2*ctypes.c_double
TotalTune = TotalTune.from_address(int(pyscsi.gv.globval.TotalTune))

# Print list of global variables.
print pyscsi.gv
print pyscsi.gv.globval
print pyscsi.gv.globval.H_exact
print TotalTune[0], TotalTune[1]
print type(pyscsi.gv.globval.TotalTune)
#print pyscsi.gv.globval.TotalTune[0], pyscsi.gv.globval.TotalTune[1]
#print pyscsi.gv.globval.CODvect[0]

print pyscsi.gv.Cell
print pyscsi.gv.Cell[0]
print pyscsi.gv.Cell[1]
#print pyscsi.gv.Cell[0].A[0]
