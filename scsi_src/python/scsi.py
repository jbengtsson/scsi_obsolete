from ctypes import *

import sys
import os

#current_dir = os.getcwd()
home_dir = os.path.expanduser('~')
print home_dir+'/git_repos/scsi/scsi_src/lib'

sys.path.append(home_dir+'/git_repos/scsi/scsi_src/lib')
import pyscsi

#libc = cdll.LoadLibrary('libc.so.6')
#libc++ =  cdll.LoadLibrary('libstdc++.so.6')
#gslcblas = CDLL('libgslcblas.so', mode=RTLD_GLOBAL)
#gsl = CDLL('libgsl.so', mode=RTLD_GLOBAL)

#scsi = CDLL(home_dir+'/git_repos/scsi/scsi_src/lib/libscsi.so')
scsi = cdll.LoadLibrary(home_dir+'/git_repos/scsi/scsi_src/lib/libscsi.so')


PLANES = 2

ss_dim = 6
DOF = ss_dim/2

Vector2 = c_double*2
Vector  = c_double*6
Matrix  = Vector*6

partsName = c_char*150
PartsKind = c_long

class elemtype(Structure):
    _fields_ = [('PName', partsName),
                ('PL',    c_double),
                ('Pkind', PartsKind),
                ('UU',    c_void_p)]

class CellType(Structure):
    _fields_ = [('Fnum',     c_int),
                ('Knum',     c_int),
                ('S',        c_double),
                ('next_ptr', c_void_p),
                ('dS',       Vector2),
                ('dT',       Vector2),
                ('Elem',     elemtype),
                ('Nu',       Vector2),
                ('Alpha',    Vector2),
                ('Beta',     Vector2),
                ('Eta',      Vector2),
                ('Etap',     Vector2),
                ('BeamPos',  Vector),
                ('A',        Matrix),
                ('sigma',    Matrix),
                ('maxampl',  Vector2*PLANES)]

class globvalrec(Structure):
    _fields_ = [('dPcommon',    c_double),
                ('dPparticle',  c_double),
                ('delta_RF',    c_double),
                ('TotalTune',   Vector2),
                ('Omega',       c_double),
                ('U0',          c_double),
                ('Alphac',      c_double),
                ('Chrom',       Vector2),
                ('Energy',      c_double),
                ('Cell_nLoc',   c_long),
                ('Elem_nFam',   c_long),
                ('CODimax',     c_long),

                ('CODeps',      c_double),
                ('CODvect',     Vector),
                ('bpm',         c_int),
                ('hcorr',       c_int),
                ('vcorr',       c_int),
                ('qt',          c_int),
                ('gs',          c_int),
                ('ge',          c_int),
                ('OneTurnMat',  Matrix),
                ('Ascr',        Matrix),
                ('Ascrinv',     Matrix),
                ('Vr',          Matrix),
                ('Vi',          Matrix),

                ('MatMeth',     c_bool),
                ('Cavity_on',   c_bool),
                ('radiation',   c_bool),
                ('emittance',   c_bool),
                ('quad_fringe', c_bool),
                ('H_exact',     c_bool),
                ('pathlength',  c_bool),
                ('stable',      c_bool),
                ('Aperture_on', c_bool),
                ('EPU',         c_bool),
                ('wake_on',     c_bool),

                ('dE',          c_double),
                ('alpha_rad',   c_double*DOF),
                ('D_rad',       c_double*DOF),
                ('J',           c_double*DOF),
                ('tau',         c_double*DOF),
                ('IBS',         c_double),
                ('Qb',          c_bool),
                ('D_IBS',       c_double*DOF),
                ('wr',          Vector),
                ('wi',          Vector),
                ('eps',         c_double*DOF),
                ('epsp',        c_double*DOF),
                ('alpha_z',     c_double),
                ('beta_z',      c_double),
                ('RingType',    c_int)
                ]


#globval = cast(scsi.globval, POINTER(globvalrec))[0]
globval = globvalrec.in_dll(scsi, 'globval')

#Cell = POINTER(CellType).in_dll(scsi, 'Cell')
Cell = cast(scsi.Cell, POINTER(CellType))

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

print
print Cell[7].Elem.PName, Cell[7].Elem.Pkind

print
print pyscsi.ElemIndex('SL1G2C01A')
