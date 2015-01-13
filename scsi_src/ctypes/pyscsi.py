from ctypes import *

import sys
import os

#current_dir = os.getcwd()
home_dir = os.path.expanduser('~')

#libc = cdll.LoadLibrary('libc.so.6')
#libc++ =  cdll.LoadLibrary('libstdc++.so.6')
#gslcblas = CDLL('libgslcblas.so', mode=RTLD_GLOBAL)
#gsl = CDLL('libgsl.so', mode=RTLD_GLOBAL)

#scsi = CDLL(home_dir+'/git_repos/scsi/scsi_src/lib/libscsi.so')
scsi = cdll.LoadLibrary(home_dir+'/git_repos/scsi/scsi_src/lib/libscsi.so')

#sys.path.append(home_dir+'/git_repos/scsi/scsi_src/lib/_pyscsi.so')
#import pyscsi

Vector2 = c_double*2
Vector  = c_double*6
Matrix  = Vector*6

partsName = c_char*150
PartsKind = c_long

class elemtype(Structure):
    _fields_ = (('PName', partsName),
                ('PL',    c_double),
                ('PKind', PartsKind),
                ('UU',    c_void_p))

class CellType(Structure):
    _fields_ = (('Fnum',    c_long),
                ('Knum',    c_long),
                ('S',       c_double),
                ('dS',      Vector2),
                ('dT',      Vector2),
                ('Elem',    elemtype),
                ('Nu',      Vector2),
                ('Alpha',   Vector2),
                ('Beta',    Vector2),
                ('Eta',     Vector2),
                ('Etap',    Vector2),
                ('BeamPos', Vector))

class globvalrec(Structure):
    _fields_ = [('dPcommon',   c_double),
                ('dPparticle', c_double),
                ('delta_RF',   c_double),
                ('TotalTune',  Vector2),
                ('Omega',      c_double),
                ('U0',         c_double),
                ('Alphac',     c_double),
                ('Chrom',      Vector2),
                ('Energy',     c_double),
                ('Cell_nLoc',  c_long),
                ('Elem_nFam',  c_long),
                ('CODimax',    c_long),

                ('CODeps',     c_double),
                ('CODvect',    Vector),
                ('bpm',        c_int),
                ('hcorr',      c_int),
                ('vcorr',      c_int),
                ('qt',         c_int),
                ('gs',         c_int),
                ('ge',         c_int),
                ('OneTurnMat', Matrix),
                ('Ascr',       Matrix),
                ('Ascrinv',    Matrix),
                ('Vr',         Matrix),
                ('Vi',         Matrix)
                ]


#globval = cast(scsi.globval, POINTER(globvalrec))[0]
globval = globvalrec.in_dll(scsi, 'globval')

scsi._Z12Read_LatticePc(
    c_char_p(home_dir+'/git_repos/scsi/scsi_src/glps/tracy_1'))

scsi._Z13Ring_GetTwissbd(c_bool(True), c_double(0.0))
   
scsi._Z9printglobv()

sys.stdout.write('\n')
for i in range(0, 6):
    for j in range(0, 6):
        sys.stdout.write('%14.6e' % globval.OneTurnMat[i][j])
    sys.stdout.write('\n')
