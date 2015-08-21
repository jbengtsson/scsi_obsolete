import math
import numpy
import re

# Module to translate from ELEGANT to Tracy-2,3 lattice.

home_dir = '/home/bengtsson/vladimir/'

def rd_lines(file_name):
    for line in open(home_dir+file_name):
        line = line.rstrip('\n')
#        print '%s' % (line)
        if not line.rstrip():
            # Blank line.
            print '# blank line'
            continue
        elif line.startswith('!'):
            # Comment.
            print '# comment'
            continue
        elif line.startswith('%'):
            # Declaration.
            token = re.split(r'[;,=\s]\s*', line)
            print '# definition: %-10s %-10s' % (token[1], token[3])
            continue
        else:
            token = re.split(r'[;,=\s]\s*', line)
            if token[0].endswith(':'):
                print 'assignment', token[1]


rd_lines('lattice2p0_v1_20150522.lte')
