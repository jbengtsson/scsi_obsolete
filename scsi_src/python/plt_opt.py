from ctypes import *
from pylab  import *

import sys
import os

from scsi import *


twopi = 2.0*math.pi

drift = 0; Mpole = 2; Quad = 2; Sext = 3
X_ = 0; Y_ = 1
x_ = 0; px_ = 1; y_ = 2; py_ = 3; delta_ = 4; ct_ = 5

c1 = 1.0/(2.0*(2.0-2.0**(1.0/3.0))); c2 = 0.5 - c1
d1 = 2.0*c1; d2 = 1.0 - 2.0*d1


def getS():
    s = zeros(globval.Cell_nLoc+1)
    for k in range(0, globval.Cell_nLoc+1):
        s[k] = Cell[k].S
    return s


def getBetaX():
    betax = zeros(globval.Cell_nLoc+1)
    for k in range(0, globval.Cell_nLoc+1):
        betax[k] = Cell[k].Beta[X_];
    return betax


def getBetaY():
    betay = zeros(globval.Cell_nLoc+1)
    for k in range(0, globval.Cell_nLoc+1):
        betay[k] = Cell[k].Beta[Y_];
    return betay


def getPhiX():
    phix = zeros(globval.Cell_nLoc+1)
    for k in range(0, globval.Cell_nLoc+1):
        phix[k] = Cell[k].Nu[X_];
    return phix


def getPhiY():
    phiy = zeros(globval.Cell_nLoc+1)
    for k in range(0, globval.Cell_nLoc+1):
        phiy[k] = Cell[k].Nu[Y_];
    return phiy


def getEtaX():
    etax = zeros(globval.Cell_nLoc+1)
    for k in range(0, globval.Cell_nLoc+1):
        etax[k] = Cell[k].Eta[X_];
    return etax


def getEtaY():
    etay = zeros(globval.Cell_nLoc+1)
    for k in range(0, globval.Cell_nLoc+1):
        etay[k] = Cell[k].Eta[Y_];
    return etay


def get_code(k):
    if Cell[k].Elem.Pkind == drift:
        code = 0.0
    elif Cell[k].Elem.Pkind == Mpole:
        if Cell[k].Elem.M[0].Pirho != 0.0:
            code = 0.5
        elif Cell[k].Elem.M[0].n_design == Quad:
            (b2, a2) = pyscsi.get_bn_design_elem(
                Cell[k].Fnum, Cell[k].Knum, Quad)
            code = math.copysign(1, b2)
        elif Cell[k].Elem.M[0].n_design == Sext:
            (b3, a3) = pyscsi.get_bn_design_elem(
                Cell[k].Fnum, Cell[k].Knum, Sext)
            code = 1.5*math.copysign(1, b3)
        elif Cell[k].Fnum == globval.bpm:
            code = 2.0
        else:
            code = 0.0
    else:
        code = 0.0

    return code


def propagate_optics(i, n):
    ss = zeros(n+1); codes = zeros(n+1)
    alphas = zeros((n+1, 2)); betas = zeros((n+1, 2))
    nus = zeros((n+1, 2)); etas = zeros((n+1, 4))

    n1 = 1; j = 0 if i == 0 else i-1

    ss[0] = Cell[j].S; codes[0] = get_code(j)
    alphas[0] = array(Cell[j].Alpha)
    betas[0] = array(Cell[j].Beta)
    nus[0] = array(Cell[j].Nu)
    etas[0] = array([
        Cell[j].Eta[X_], Cell[j].Etap[X_],
        Cell[j].Eta[Y_], Cell[j].Etap[Y_]])

    L = Cell[i].Elem.PL
    if (i != 0) and \
           ((Cell[i].Elem.Pkind == drift) or
            ((Cell[i].Elem.Pkind == Mpole) and (L != 0.0))):
        n1 = n; s = Cell[i].S - L;

        alpha = Cell[i-1].Alpha; beta = Cell[i-1].Beta
        nu = Cell[i-1].Nu
        eta = Cell[i-1].Eta; etap = Cell[i-1].Etap

        h = L/n; A = get_A(alpha, beta, eta, etap)
        for j in range(0, n):
            s += h

            if Cell[i].Elem.Pkind == drift:
                A = Drift(h, A)
            elif Cell[i].Elem.Pkind == Mpole:
                Mp = Cell[i].Elem.M[0]

                if j == 1 and Mp.Pirho != 0.0:
                    A = EdgeFocus(Mp.Pirho, Mp.PTx1, Mp.Pgap, A)

                A = Drift(c1*h, A)
                A = thin_kick(
                    Mp.Porder, Mp.PB, d1*h, Mp.Pirho, Mp.Pirho, A)
                A = Drift(c2*h, A)
                A = thin_kick(
                    Mp.Porder, Mp.PB, d2*h, Mp.Pirho, Mp.Pirho, A)
                A = Drift(c2*h, A)
                A = thin_kick(
                    Mp.Porder, Mp.PB, d1*h, Mp.Pirho, Mp.Pirho, A)
                A = Drift(c1*h, A)

                if j == n and Mp.Pirho != 0.0:
                    A = EdgeFocus(Mp.Pirho, Mp.PTx2, Mp.Pgap, A)

            [alpha, beta, dnu, eta, etap] = get_ab(A)
            dnu = array(dnu)

            if L < 0.0:
                dnu = dnu - (1.0, 1.0)

            ss[j] = s; codes[j] = get_code(i)
            alphas[j] = alpha; betas[j] = beta; nus[j] = nu+dnu
            etas[j] = array([eta[X_], etap[X_], eta[Y_], etap[Y_]])

    return (n1, ss, codes, alphas, betas, nus, etas)


def prt_lat(fname, n):
    outf = open(fname, 'w')
    outf.write('#        name           s   code'
               '  alphax  betax   nux   etax   etapx')
    outf.write('  alphay  betay   nuy   etay   etapy\n')
    outf.write('#                      [m]'
               '                 [m]           [m]')
    outf.write('                   [m]           [m]\n')
    outf.write('#\n')

    for i in range(0, globval.Cell_nLoc+1):
        [n1, ss, codes, alphas, betas, nus, etas] = propagate_optics(i, n)

        for k in range(0, n1):
            outf.write('%4ld %15s %6.2f %4.1f'
                       ' %7.3f %6.3f %6.3f %6.3f %6.3f'
                       ' %7.3f %6.3f %6.3f %6.3f %6.3f\n' %
                       (i, Cell[i].Elem.PName, ss[k], codes[k],
                        alphas[k][X_], betas[k][X_], nus[k][X_],
                        etas[k][x_], etas[k][px_],
                        alphas[k][Y_], betas[k][Y_], nus[k][Y_],
                        etas[k][y_], etas[k][py_]))

    outf.close()


def get_codes():
    codes = []
    for k in range(0, globval.Cell_nLoc+1):
        codes.append(get_code(k))
    return codes


def plt_opt(displ):
    s = getS(); code = get_codes()
    beta = [getBetaX(), getBetaY()]
    nu = [getPhiX(), getPhiY()]
    eta = [getEtaX(), getEtaY()]

    plt.rcParams['savefig.dpi'] = 600 # For png.

    fig = plt.figure()
    fig.subplots_adjust(hspace=0.4) # default is 0.2.

    ax = fig.add_subplot(211)
    ax.grid(True)
    ax.set_title(r'$\beta$'); ax.set_xlabel('s [m]'); ax.set_ylabel('[m]')
    ax2 = ax.twinx(); ax2.set_ylim(-1.5, 20); ax2.set_autoscaley_on(False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax.plot(s, beta[0], label=r'$\beta_x$', color='blue')
    ax.plot(s, beta[1], label=r'$\beta_y$', color='red')
    ax2 = plt.step(s, code, color='black')
    # Legend must be created first.
    ax.legend(loc='upper right')

    ax = fig.add_subplot(212)
    ax.grid(True)
    ax.set_title(r'$\eta$'); ax.set_xlabel('s [m]'); ax.set_ylabel('[m]')
    ax2 = ax.twinx(); ax2.set_ylim(-1.5, 20); ax2.set_autoscaley_on(False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax.plot(s, eta[0], label=r'$\eta_x$', color='blue')
    ax.plot(s, eta[1], label=r'$\eta_y$', color='red')
    ax2 = plt.step(s, code, color='black')
    # Legend must be created first.
    ax.legend(loc='upper right')

##     fig.savefig('optics.png')
    fig.savefig('optics.ps', orientation='landscape')

    if displ:
        plt.ion(); plt.show(); plt.ioff()
        raw_input('<ret> to continue>')


sys.stdout.write('\n')
pyscsi.Read_Lattice('/home/bengtsson/git_repos/scsi/scsi_src/glps/tracy_1')

pyscsi.Ring_GetTwiss(True, 0)

sys.stdout.write('\n')
sys.stdout.write('C   = %7.5f [m]\n'
                 % (Cell[globval.Cell_nLoc].S))
sys.stdout.write('nu  = [%7.5f, %7.5f]\n' % \
                     (globval.TotalTune[0], globval.TotalTune[1]))
sys.stdout.write('ksi = [%5.3f, %5.3f]\n' % \
                     (globval.Chrom[0], globval.Chrom[1]))

#pyscsi.prt_lat('linlat.out', 10)

plt_opt(True)
