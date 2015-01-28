from scsi import *

pyscsi.Read_Lattice(home_dir+'/git_repos/scsi/scsi_src/src/glps/tracy_1')

pyscsi.Ring_GetTwiss(True, 0.0); pyscsi.printglob()

pyscsi.prt_lat('linlat.out', globval.bpm, True);

sys.stdout.write('\n')
for i in range(0, 6):
    for j in range(0, 6):
        sys.stdout.write('%14.6e' % globval.OneTurnMat[i][j])
    sys.stdout.write('\n')

print
print globval.radiation

sys.stdout.write('\n')
for k in range(0, 5):
    sys.stdout.write('%10s %14.6e %14.6e\n' %
                     (Cell[k].Elem.contents.name,
                      Cell[k].Beta[X_], Cell[k].Beta[Y_]))

Fnum = pyscsi.ElemIndex('SL1G2C01A'); loc = pyscsi.Elem_GetPos(Fnum, 1)
print
print Fnum, loc

print
print Cell[loc].Elem.contents.name, Cell[loc].Elem.contents.kind
print Cell[loc].Elem.deref('M').method, Cell[loc].Elem.deref('M').n
print Cell[loc].Elem.deref('M').order, Cell[loc].Elem.deref('M').n_design, \
      Cell[loc].Elem.deref('M').thick, \
      Cell[loc].Elem.deref('M').bnpar[HOMmax+Sext]

Fnum = pyscsi.ElemIndex('CAV'); loc = pyscsi.Elem_GetPos(Fnum, 1)

print
print Cell[loc].Elem.name, Cell[loc].Elem.kind
print Cell[loc].Elem.deref('C').volt, Cell[loc].Elem.deref('C').freq
