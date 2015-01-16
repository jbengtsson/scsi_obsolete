/*
 * Copyright (C) 2008 Lingyun Yang.
 *
 * For more information please contact lingyun.yang@gmail.com
 *
 */

/**
 * This is an example of how to use the general lattice parser, i.e. glps
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <cstring>

#include "glps.h"
#include "glpsutil.h"
extern int yy_flex_debug;
extern int yydebug;

using namespace std;

int main(int argc, char* argv[])
{
    if (argc == 1) {
        cout << "teslaparse [-p|-s] lattice.lat" << endl;
        return 0;
    }

    // get the version
    int vmaj, vmin, vpatch, vrevision;
    glps_version(&vmaj, &vmin, &vpatch, &vrevision);
    
    for (int i = 1; i < argc - 1; ++i) {
        if (argv[i] == string("-s")) yy_flex_debug = 1;
        if (argv[i] == string("-p")) yydebug = 1;
    }

    // parse the input lattice from command line.
    int ret = parse_lattice(argv[argc - 1]);
    if ( ret ) {
        cerr << "ERROR: " << ret << endl;
        free_lattice();
        return ret;
    }

    cout << "# glps version " << vmaj << '.' << vmin << '.' << vpatch << endl;

    list<string> elems = glps_elements();
    string str;
    for (list<string>::iterator i = elems.begin(); i != elems.end(); ++i) {
        glps_read(*i, str);
        cout << *i << ": " << str << "  ";
        list<string> prpt = glps_properties(*i);
        
        for (list<string>::iterator k=prpt.begin(); k != prpt.end(); ++k) {
            cout << *k << ", ";
        }
        cout << ";" << endl;
    }
    
    
    list<string> lines = glps_lines();
    for (list<string>::iterator i = lines.begin(); i != lines.end(); ++i) {
        cout << *i << ": " << endl;
        list<string> elems = glps_elements(*i, true);
        for (list<string>::iterator j = elems.begin(); j != elems.end(); ++j) {
            size_t pt = (*j).rfind('.');
            string tmp = (*j).substr(pt+1);
            list<string> prpt = glps_properties(tmp);
            cout << "    " << *j << " " << tmp;            
            for (list<string>::iterator k=prpt.begin(); k != prpt.end(); ++k) {
                cout << " " << *k;
            }
            cout << endl;
        }
        cout << ";" << endl;
    }

    free_lattice();
}
