from ctypes import *

import sys
import os

#current_dir = os.getcwd()
home_dir = os.path.expanduser('~')

sys.path.append(home_dir+'/git_repos/scsi/scsi_src/lib')
import pyscsi

#libc = cdll.LoadLibrary('libc.so.6')
#libc++ =  cdll.LoadLibrary('libstdc++.so.6')
#gslcblas = CDLL('libgslcblas.so', mode=RTLD_GLOBAL)
#gsl = CDLL('libgsl.so', mode=RTLD_GLOBAL)

#scsi = CDLL(home_dir+'/git_repos/scsi/scsi_src/lib/libscsi.so')
scsi = cdll.LoadLibrary(home_dir+'/git_repos/scsi/scsi_src/lib/libscsi.so')


drift = 0; Mpole = 2; Quad = 2; Sext = 3
X_ = 0; Y_ = 1
x_ = 0; px_ = 1; y_ = 2; py_ = 3; delta_ = 4; ct_ = 5

c1 = 1.0/(2.0*(2.0-2.0**(1.0/3.0))); c2 = 0.5 - c1
d1 = 2.0*c1; d2 = 1.0 - 2.0*d1


PLANES = 2

ss_dim = 6; DOF = ss_dim/2

Vector2    = c_double*2
Vector     = c_double*6
Matrix     = Vector*6

partsName  = c_char*150
PartsKind  = c_long

HOMmax     = 21
mpolArray  = c_double*(2*HOMmax+1)
pthicktype = c_long

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
                ('RingType',    c_int)]

class DriftType(Structure):
    _fields_ = []

class MpoleType(Structure):
    _fields_ = [('Pmethod', c_int),
                ('PN',      c_int),

                ('PdSsys',  Vector2),
                ('PdSrms',  Vector2),
                ('PdSrnd',  Vector2),

                ('PdTpar',  c_double),
                ('PdTsys',  c_double),
                ('PdTrms',  c_double),
                ('PdTrnd',  c_double),

                ('PBpar',    mpolArray),
                ('PBsys',    mpolArray),
                ('PBrms',    mpolArray),
                ('PBrnd',    mpolArray),
                ('PB',       mpolArray),
                ('Porder',   c_int),
                ('n_design', c_int),
                ('Pthick',   pthicktype),
                ('PTx1',     c_double),
                ('PTx2',     c_double),
                ('Pgap',     c_double),
                ('Pirho',    c_double),
                ('Pc0',      c_double),
                ('Pc1',      c_double),
                ('Ps1',      c_double)]

class CavityType (Structure):
    _fields_ = [('Pvolt', c_double),
                ('Pfreq', c_double),
                ('phi',   c_double),
                ('Ph',    c_int)]

class elemtype(Structure):
    def deref(self, type):
        if type == 'D':
            return cast(self.U, POINTER(DriftType))[0]
        elif type == 'M':
            return cast(self.U, POINTER(MpoleType))[0]
        elif type == 'C':
            return cast(self.U, POINTER(CavityType))[0]
        else:
            print "deref: undef. type:", type
            exit(1)

    _fields_ = [('PName', partsName),
                ('PL',    c_double),
                ('Pkind', PartsKind),
                ('U',     POINTER(c_void_p))]

class CellType(Structure):
    _fields_ = [('Fnum',     c_int),
                ('Knum',     c_int),
                ('S',        c_double),
                ('next_ptr', POINTER(c_void_p)),
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

#globval = cast(scsi.globval, POINTER(globvalrec))[0]
globval = globvalrec.in_dll(scsi, 'globval')

Cell = cast(scsi.Cell, POINTER(CellType))
# Returns the same type and address but, somehow, does not work.
#Cell = POINTER(CellType).in_dll(scsi, 'Cell')

print pyscsi.get_code(pyscsi.gv.Cell[5])
