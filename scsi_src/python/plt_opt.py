from scsi  import *
from pylab import *


def get_code(k):
    if Cell[k].Elem.Pkind == drift:
        code = 0.0
    elif Cell[k].Elem.Pkind == Mpole:
        if Cell[k].Elem.deref('M').Pirho != 0.0:
            code = 0.5
        elif Cell[k].Elem.deref('M').n_design == Quad:
            (b2, a2) = pyscsi.get_bn_design_elem(
                Cell[k].Fnum, Cell[k].Knum, Quad)
            code = math.copysign(1, b2)
        elif Cell[k].Elem.deref('M').n_design == Sext:
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


def get_opt():
    s = zeros(globval.Cell_nLoc+1)
    code = zeros(globval.Cell_nLoc+1)
    beta = zeros((2, globval.Cell_nLoc+1))
    nu = zeros((2, globval.Cell_nLoc+1))
    eta = zeros((2, globval.Cell_nLoc+1))
    for k in range(0, globval.Cell_nLoc+1):
        s[k] = Cell[k].S;
        code[k] = get_code(k)
#        code[k] = pyscsi.get_code(byref(Cell[k]))
        beta[X_, k] = Cell[k].Beta[X_]; beta[Y_, k] = Cell[k].Beta[Y_]
        nu[X_, k] = Cell[k].Nu[X_]; nu[Y_, k] = Cell[k].Nu[Y_]
        eta[X_, k] = Cell[k].Eta[X_]; eta[Y_, k] = Cell[k].Eta[Y_]

    return(s, code, beta, nu, eta)


def plt_opt(displ):
    (s, code, beta, nu, eta) = get_opt()

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

pyscsi.prt_lat('linlat1.out', globval.bpm, True)
pyscsi.prt_lat('linlat.out', globval.bpm, True, 10)

plt_opt(True)
