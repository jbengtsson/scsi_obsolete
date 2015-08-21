import math
import numpy
import re
import StringIO

# Module to translate from ELEGANT to Tracy-2,3 lattice.

home_dir = '/home/bengtsson/vladimir/'

# Elegant -> Tracy-2,3 dictionary.
ele2tracy = { 'watch' : 'marker', 'charge' : 'N/A'}


def rd_lines(file_name):
    str = file_name.split('.')[0]+'.lat'
    f = open(str, 'w')
    for line in open(home_dir+file_name):
        line = line.strip('\n').lower()
#        print '%s' % (line)
        if not line.rstrip():
            # Blank line.
            f.write('\n')
            continue
        elif line.startswith('!'):
            # Comment.
            f.write('{ %s }\n' % (line.strip('!')))
            continue
        elif line.startswith('%'):
            # Declaration.
            token = re.split(r'[;,=\s]\s*', line)
            f.write('%s = %s;\n' % (token[3], token[1]))
        else:
            token = re.split(r'[;,=\s]\s*', line)
            if token[0].endswith(':'):
                if token[0] == 'c:':
                    f.write('{ %s }' % (line))
                else:
                    f.write('%s %s;\n' % (token[0], ele2tracy[token[1]]))
                    continue


rd_lines('lattice2p0_v1_20150522.lte')
