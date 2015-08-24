import math
import numpy
import re
import StringIO

# Module to translate from ELEGANT to Tracy-2,3 lattice.

# Elegant -> Tracy-2,3 dictionary.
ele2tracy = {
    'charge'    : 'Marker',
    'mark'      : 'Marker',
    'watch'     : 'Marker',
    'ematrix'   : 'Marker',
    'rfca'      : 'Cavity',
    'drif'      : 'Drift',
    'drift'     : 'Drift',
    'csrdrif'   : 'Drift',
    'rcol'      : 'Drift',
    'csrcsbend' : 'Bending',
    'quad'      : 'Quadrupole',
    'sext'      : 'Sextupole',
    'freq'      : 'Frequency'
    }


def parse_rpnc(stack):
    # Reverse Polish Notation calculator.
    arg = stack.pop()
    if arg in '+-*/':
        op = arg
        arg = parse_rpnc(stack)
        return parse_rpnc(stack) + (' %s %s' % (op, arg))
    else:
        return ('%s' % (arg))


def parse_decl(decl):
    stack = (decl.strip('%')).split()
    lhs = stack.pop()
    op = stack.pop()
    return ('%s = ' % (lhs)) + parse_rpnc(stack)


def get_index(tokens, token):
    return tokens.index(token) if token in tokens else None


def get_arg(str):
    if str.find('"') == -1:
        arg = str
    else:
        arg = parse_rpnc((str.strip('"')).split())
    return arg


def parse_definition(line, tokens):
    for k in range(len(tokens)):
        # Remove white space; unless a string.
        if not tokens[k].startswith('"'):
            tokens[k] = re.sub('[\s*]', '', tokens[k])
    str = ''
    if tokens[1] == 'twiss':
        pass
    elif tokens[1] == 'charge':
        str = '%s: %s; { %s }' % (tokens[0], ele2tracy[tokens[1]], line)
    elif tokens[1] == 'mark':
        str = '%s: %s;' % (tokens[0], ele2tracy[tokens[1]])
    elif tokens[1] == 'watch':
        str = '%s: %s;' % (tokens[0], ele2tracy[tokens[1]])
    elif tokens[1] == 'drift' or  tokens[1] == 'drif':
        loc_l = tokens.index('l')
        str = '%s: %s, L = %s;' % \
            (tokens[0], ele2tracy[tokens[1]], get_arg(tokens[loc_l+1]))
    elif tokens[1] == 'csrdrif' or tokens[1] == 'csrdrift':
        loc_l = tokens.index('l')
        str = '%s: %s, L = %s;' % \
            (tokens[0], ele2tracy[tokens[1]], get_arg(tokens[loc_l+1]))
    elif tokens[1] == 'rcol':
        loc_l = tokens.index('l')
        str = '%s: %s, L = %s;' % \
            (tokens[0], ele2tracy[tokens[1]], get_arg(tokens[loc_l+1]))
    elif tokens[1] == 'csrcsbend':
        loc_l = tokens.index('l')
        loc_phi = tokens.index('angle')
        loc_e1 = get_index(tokens, 'e1')
        loc_e2 = get_index(tokens, 'e2')
        loc_k = get_index(tokens, 'k')
        str = '%s: %s, L = %s, T = %s' % \
            (tokens[0], ele2tracy[tokens[1]], get_arg(tokens[loc_l+1]),
             get_arg(tokens[loc_phi+1]))
        if loc_e1:
            str += ', T1 = %s' % (get_arg(tokens[loc_e1+1]))
        if loc_e2:
            str += ', T2 = %s' % (get_arg(tokens[loc_e2+1]))
        if loc_k:
            str += ', K = %s' % (get_arg(tokens[loc_k+1]))
        str += ', N = Nbend, Method = 4;'
    elif tokens[1] == 'quad':
        loc_l = tokens.index('l')
        loc_k = tokens.index('k1')
        str = '%s: %s, L = %s, K = %s, N = Nquad, Method = 4;' % \
            (tokens[0], ele2tracy[tokens[1]], get_arg(tokens[loc_l+1]),
             get_arg(tokens[loc_k+1]))
    elif tokens[1] == 'sext':
        loc_l = tokens.index('l')
        loc_k = tokens.index('k2')
        str = '%s: %s, L = %s, K = %s, N = Nsext, Method = 4;' % \
            (tokens[0], ele2tracy[tokens[1]], get_arg(tokens[loc_l+1]),
             get_arg(tokens[loc_k+1]))
    elif tokens[1] == 'rfca':
        loc_l = tokens.index('l')
        loc_f = tokens.index('freq')
        loc_v = tokens.index('volt')
        loc_phi = tokens.index('phase')
        cav_name = tokens[0]+'_c'
        str = '%s: %s, Frequency = %s, Voltage = %s;\n' % \
            (cav_name, ele2tracy[tokens[1]], get_arg(tokens[loc_f+1]),
        get_arg(tokens[loc_v+1]))
        drift_name = tokens[0]+'_d'
        str += '%s: Drift, L = %s;\n' % \
            (drift_name, get_arg(tokens[loc_l+1])+'/2')
        str += '%s: %s, %s, %s;' % \
            (tokens[0], drift_name, cav_name, drift_name)
    elif tokens[1] == 'line':
        str = '%s: %s' % (tokens[0], tokens[2].strip('('))
        n = len(tokens)
        for k in range(3, n-1):
            str += ', %s' % (tokens[k])
        str += ', %s;' % (get_arg(tokens[n-1].strip(')')))
    elif tokens[1] == 'ematrix':
        pass
    else:
        print '*** undefined token'
        print line
        print tokens
        exit(1)
    return str


def parse_line(line, outf):
    line_lc = line.lower()
    if not line_lc.rstrip():
        # Blank line.
        outf.write('\n')
    elif line_lc.startswith('!'):
        # Comment.
        outf.write('{ %s }\n' % (line.strip('!')))
    elif line_lc.startswith('%'):
        # Declaration.
        outf.write('%s;\n' % (parse_decl(line_lc.strip('%'))))
    else:
        tokens = re.split(r'[,:=]', line_lc)
        if line_lc.find(':') != -1:
            # Definition.
            outf.write('%s\n' % (parse_definition(line_lc, tokens)))


def prt_decl(outf):
    outf.write('define lattice; ringtype = 1;\n')
    outf.write('\n')
    outf.write('Energy = 3.0;\n')
    outf.write('\n')
    outf.write('dP = 1e-8; CODeps = 1e-14;\n')
    outf.write('\n')
    outf.write('Meth = 4; Nbend = 4; Nquad = 4; Nsext = 2;\n')
    outf.write('\n')
    outf.write('pi = 4.0*arctan(1.0); c0 = 2.99792458e8;\n')
    outf.write('\n')


def transl_file(file_name):
    str = file_name.split('.')[0]+'.lat'
    inf = open(home_dir+file_name, 'r')
    outf = open(str, 'w')
    prt_decl(outf)
    line = inf.readline()
    while line:
        line = line.strip('\r\n')
        while line.endswith('&'):
            # Line
            line = line.strip('&')
            line += (inf.readline()).strip('\r\n')
        parse_line(line, outf)
        line = inf.readline()
    outf.write('\n')
    outf.write('cell: ring, symmetry = 1;\n')
    outf.write('\n')
    outf.write('end;\n')


home_dir = '/home/bengtsson/vladimir/'

transl_file('lattice2p0_v1_20150522.lte')
