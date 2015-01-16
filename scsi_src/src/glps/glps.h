/*
 * Copyright (C) 2008 Lingyun Yang.
 *
 * This software may be used and distributed according to the terms
 * of the GNU General Public License, incorporated herein by reference.
 *
 * For more information please contact lyyang@bnl.gov
 *
 */


/* Author: Lingyun Yang, lyyang@bnl.gov */

#include <stdlib.h>
#include <stdio.h>
/* #include "elemtree.h" */


#ifndef GLPS_H_
#define GLPS_H_

#ifdef __cplusplus
extern "C" {
 #endif 

 #define GLPS_VERSION "2.0.0"
 #define glps_VERSION_MAJOR 2   
 #define glps_VERSION_MINOR 0
 #define glps_VERSION_PATCH 0

    struct elemtree_el;
    /* typedef struct elemtree_el elemtree; */

    enum GLPS_DEF_TYPE {
        GLPS_NONE = -1,         /* */
        GLPS_NODE_UNDEF = -2,   /* not defined */
        GLPS_ASSIGN = 128,      /* assignment */
        GLPS_LINE,              /* beam line */
        GLPS_ACTION,            /* action */
        GLPS_ELEMENT,           /* element */
        GLPS_NODE_TAG,          /* a tag, no value */
        GLPS_NODE_REAL,         /* float/double value */
        GLPS_NODE_FUNC1,        /* function, f(x) */
        GLPS_NODE_FUNC2,        /* function, f(x,y) */
        GLPS_NODE_FUNC3,        /* function, f(x,y,z) */
        GLPS_NODE_VEC,          /* vector */
        GLPS_NODE_STR,          /* string */
        GLPS_ERROR              /* */
    } ;
    

    /*! Global symbol table, store every variable, function defined. It also
     *  stores the property of magnets and actions. But since the magnet
     *  properties are copied into a new structure after matched in Yacc
     *  grammar, and freed after that, they are seperated from the global
     *  linked list.
     */
    struct symbol_table_rec {
        int            nodetype;
        double     (*funcptr)();    /* cos, sin, ... */
        char              *name;    /* PI, L, ... */
        size_t             size;    /* */
        union {
            double            value;    /* */
            double             *vec;    /* */
            char               *str;    /* */
            struct elemtree_el *line;
        };

        struct symbol_table_rec *property;
        struct symbol_table_rec *next;
    };
    typedef struct symbol_table_rec GLPS_SYMB;


    /*!\brief parse the lattice file into internal data structure 
     */
    int parse_lattice(const char* f);

    /*! \brief release the internal structure storing parsed lattice 
     */
    int free_lattice();

    char *get_path_root(char *f);
    int parse_lattice_flat(char* f, char* flat);
    int parse_lattice_xml(char* lat, char* xml);
    int parse_lattice_ini(char* f);
    int print_flat_lattice(FILE* pf, char* bl);


    /* ----------------------------------------- */

    /*! \brief get the definition of a statement, including action, element
     *   and beamline 
     */
    /*! \brief compare two string, case-incensitive */
    int str_case_cmp(const char* ps1, const char* ps2);
    void glps_version(int* vmaj, int* vmin, int* vpat, int* vrev);
    void glps_error (int glps_errno, const char * reason);

    /* symb linked list operation */
    GLPS_SYMB* symb_lookup(GLPS_SYMB* head, const char* name);

    /* read data in symb linked list */
    int symb_get_double(char* name, double* val);
    int symb_get_vector(char* name, double** val, size_t *n);
    int symb_get_string(char* name, char** val);

#ifdef __cplusplus
}
#endif

#endif
