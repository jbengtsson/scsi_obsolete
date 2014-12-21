import pyscsi


pyscsi.gv.globval.H_exact    = False; pyscsi.gv.globval.quad_fringe = False;
pyscsi.gv.globval.Cavity_on  = False; pyscsi.gv.globval.radiation   = False;
pyscsi.gv.globval.emittance  = False; pyscsi.gv.globval.IBS         = False;
pyscsi.gv.globval.pathlength = False; pyscsi.gv.globval.bpm         = 0;

# Print list of global variables.
print pyscsi.gv

pyscsi.Read_Lattice('/home/bengtsson/projects/in/lattice/commissioning')

pyscsi.Ring_GetTwiss(True, 0.0); pyscsi.printglob()
