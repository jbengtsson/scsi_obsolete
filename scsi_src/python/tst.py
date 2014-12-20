import pyscsi


pyscsi.var.globval.H_exact    = False; pyscsi.var.globval.quad_fringe = False;
pyscsi.var.globval.Cavity_on  = False; pyscsi.var.globval.radiation   = False;
pyscsi.var.globval.emittance  = False; pyscsi.var.globval.IBS         = False;
pyscsi.var.globval.pathlength = False; pyscsi.var.globval.bpm         = 0;

print pyscsi.var.globval.H_exact

pyscsi.Read_Lattice('/home/bengtsson/projects/in/lattice/commissioning')

pyscsi.Ring_GetTwiss(True, 0.0); pyscsi.printglob()
