import math
import numpy
import re
import StringIO

# Module to translate from ELEGANT to Tracy-2,3 lattice.

home_dir = '/home/bengtsson/vladimir/'

# Elegant -> Tracy-2,3 dictionary.
ele2tracy = {
    'charge'    : 'N/A',
    'mark'      : 'marker',
    'watch'     : 'marker',
    'ematrix'   : 'marker',
    'rfca'      : 'cavity',
    'drif'      : 'drift',
    'drift'     : 'drift',
    'csrdrif'   : 'drift',
    'csrcsbend' : 'bending',
    'quad'      : 'quadrupole',
    'sext'      : 'sextupole',
    'line'      : ''
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
        outf.write('%s;\n' % (parse_decl(line.strip('%'))))
    else:
#        tokens = re.split(r'[,\s]\s*', line_lc)
        tokens = re.split(r'[,:\s]\s*', line_lc)
#        print tokens
        if line_lc <= ':':
            # Definition.
            if tokens[1] == 'twiss':
                # Ignore.
                a = 1
        elif tokens[1] == 'charge':
            # Charge definition.
            outf.write('{ %s }\n' % (line_lc))
        elif tokens[1] == 'mark':
            outf.write('%s %s;\n' % (tokens[0], ele2tracy[tokens[1]]))
        elif tokens[1] == 'watch':
            outf.write('%s %s;\n' % (tokens[0], ele2tracy[tokens[1]]))
        elif tokens[1] == 'rfca':
            outf.write('%s %s;\n' % (tokens[0], ele2tracy[tokens[1]]))
            print line_lc
            print tokens
            print tokens.index('rfca')
            exit()
        elif tokens[1] == 'line':
            str = line_lc.split('=')[1]
            outf.write('%s %s;\n' % (tokens[0], str.strip('()')))
        else:
            print line_lc
            print tokens
            exit()
            outf.write('not defined\n')
            outf.write('%s;\n' % (line_lc))


def rd_lines(file_name):
    str = file_name.split('.')[0]+'.lat'
    inf = open(home_dir+file_name, 'r')
    outf = open(str, 'w')
    line = inf.readline()
    while line:
        line = line.strip('\r\n')
        while line.endswith('&'):
            # Line
            line = line.strip('&')
            line += (inf.readline()).strip('\r\n')
            print line
#        print '%s' % (line)

        parse_line(line, outf)

        line = inf.readline()


rd_lines('lattice2p0_v1_20150522.lte')
