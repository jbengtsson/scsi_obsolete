from ctypes import *

import sys
import os

#current_dir = os.getcwd()
home_dir = os.path.expanduser('~')

#libc = cdll.LoadLibrary('libc.so.6')
#libc++ =  cdll.LoadLibrary('libstdc++.so.6')

gslcblas = CDLL('libgslcblas.so', mode=RTLD_GLOBAL)
gsl = CDLL('libgsl.so', mode=RTLD_GLOBAL)

scsi = cdll.LoadLibrary(home_dir+'/git_repos/scsi/scsi_src/lib/libscsi.so')

#pyscsi = cdll.LoadLibrary(home_dir+'/git_repos/scsi/scsi_src/lib/_pyscsi.so')

#sys.path.append(home_dir+'/git_repos/scsi/scsi_src/lib/_pyscsi.so')
#import pyscsi

Vector2 = c_double*2
vector  = c_double*6

partsName = c_char*150
PartsKind = c_long

class elemtype(Structure):
    _fields_ = (('PName', partsName),
                ('PL', c_double),
                ('PKind', PartsKind),
                ('UU', c_void_p))

class CellType(Structure):
    _fields_ = (('Fnum', c_long),
                ('Knum', c_long),
                ('S', c_double),
                ('dS', Vector2),
                ('dT', Vector2),
                ('Elem', elemtype),
                ('Nu', Vector2),
                ('Alpha', Vector2),
                ('Beta', Vector2),
                ('Eta', Vector2),
                ('Etap', Vector2),
                ('BeamPos', vector))

class globvalrec(Structure):
    nLocMax1 = 10001
    _fields_ = (('dPcommon', c_double),
                ('dPparticle', c_double),
                ('maxamplH', Vector2 * nLocMax1),
                ('maxamplV', Vector2 * nLocMax1),
                ('RFAcceptance', c_double),
                ('TotalTune', Vector2),
                ('Omega', c_double),
                ('U0', c_double),
                ('Alphac', c_double),
                ('Chrom', Vector2),
                ('Energy', c_double),
                ('Cell_nLoc', c_long))

print scsi
globval = cast(scsi.globval, POINTER(globvalrec))[0]
print globval.dPcommon, globval.TotalTune[0], globval.TotalTune[1]

#scsi.Read_Lattice(home_dir+'/git_repos/scsi/scsi_src/glps/tracy_1')
   
