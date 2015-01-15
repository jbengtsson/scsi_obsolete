from scsi import *

pyscsi.Read_Lattice(home_dir+'/git_repos/scsi/scsi_src/glps/tracy_1')

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
    sys.stdout.write('%14.6e %14.6e\n' % (Cell[k].Beta[0], Cell[k].Beta[1]))


Fnum = pyscsi.ElemIndex('SL1G2C01A'); loc = pyscsi.Elem_GetPos(Fnum, 1)

print
print Cell[loc].Elem.PName, Cell[loc].Elem.Pkind
print Cell[loc].Elem.deref('M').Pmethod, Cell[loc].Elem.deref('M').PN

Fnum = pyscsi.ElemIndex('CAV'); loc = pyscsi.Elem_GetPos(Fnum, 1)

print
print Cell[loc].Elem.PName, Cell[loc].Elem.Pkind
print Cell[loc].Elem.deref('C').Pvolt, Cell[loc].Elem.deref('C').Pfreq
