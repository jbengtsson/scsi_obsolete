/*
 * Copyright (C) 2008 Lingyun Yang.
 *
 * Author: Lingyun Yang, lyyang@bnl.gov
 */

#include "glps.h"
#include "parser.h"
#include "glpserror.h"
#include "elemtree.h"
#include <string.h>

extern GLPS_SYMB* symbtab;
extern char* lat_title;
extern char glps_err_msg[];
extern int glps_err_code;
extern YYLTYPE yylloc;

#define DBG_MSG(msg, args ...) \
 fprintf(stderr, "[%s:%d:%s]: " msg, __FILE__, __LINE__, __func__, ##args)


/** \brief Get version of parser library
 */
void glps_version(int* vmaj, int* vmin, int* vpat, int* vrev)
{
    if (vmaj != NULL) *vmaj = glps_VERSION_MAJOR;
    if (vmin != NULL) *vmin = glps_VERSION_MINOR;
    if (vpat != NULL) *vpat = glps_VERSION_PATCH;
    if (vrev != NULL) *vrev = 0;
}

/*!
 * \brief case-insensitive compare, similar return convention as strcmp.
 */
int str_case_cmp(const char* ps1, const char* ps2)
{
    char c1, c2;
    int v;
    if (ps1 == NULL && ps2 == NULL) return 0;
    else if (ps1 == NULL && ps2) return 1;
    else if (ps1 && ps2 == NULL) return 1;

    do {
        c1 = *ps1++;
        c2 = *ps2++;
        if (c1 >='A' && c1 <= 'Z') c1 += 'a' - 'A';
        if (c2 >='A' && c2 <= 'Z') c2 += 'a' - 'A';        
        v = (int) c1 - (int) c2;
    }while ((v==0) && (c1 != '\0'));

    return v;
}

GLPS_SYMB* symb_lookup_stmt(GLPS_SYMB* head, const char* name)
{
    /* search for statement: variable/action/element */
    while (head) {
        if (head->name && str_case_cmp(name, head->name) == 0) break;
        head = head->next;
    }

    return head;
}

GLPS_SYMB* symb_lookup_field(
    GLPS_SYMB* head,
    const char* name, const char* field)
{
    GLPS_SYMB* stmt = symb_lookup_stmt(head, name);

    /* */        /* did not found */
    if (!stmt) return NULL;
 
    GLPS_SYMB* p = NULL;
    /* found and need subfield */
    if (stmt->nodetype == GLPS_ACTION || stmt->nodetype == GLPS_ELEMENT) {
        p = stmt->property;
        while (p) {
            if (p->name && str_case_cmp(field, p->name) == 0) break;
            p = p->next;
        }

        if (!p) {
            /* GLPS_ERROR_MSG(GLPS_EGRAMMAR, "no field '%s' for '%s'\n",
                           field, name);
            */
            sprintf(glps_err_msg, "no field '%s' for '%s'\n", field, name);
            glps_error(GLPS_EGRAMMAR, glps_err_msg);
        }
    } 

    return p;
}

GLPS_SYMB* symb_lookup(GLPS_SYMB* head, const char *name)
{
    int ndot = 0;
    char *stmt, *field, *idx, *err;
    char *str = strdup(name);
    stmt  = strtok(str, ".");
    field = strtok(NULL, ".");
    idx   = strtok(NULL, ".");
    err   = strtok(NULL, ".");

    if ( idx || err ) {
        GLPS_WARN("parser doesnot support the form of 'elem.vecfield.index' "
                  "as %s. It will be reduced to 'elem.vecfield'",
            name);
    } 

    GLPS_SYMB *p = NULL;
    if (field) {
        p = symb_lookup_field(head, stmt, field);
    } else {
        p = symb_lookup_stmt(head, stmt);
    }

    free(str);
    return p;
}

void symb_print(GLPS_SYMB* p)
{
    int i = 0;
    static int depth = 0;
    GLPS_SYMB *r;
    while (p) {
        switch(p->nodetype) {
        case GLPS_NODE_REAL:
            fprintf(stdout, "%s: %f ", p->name, p->value);
            break;
        case GLPS_NODE_VEC:
            fprintf(stdout, "%s: ", p->name);
            for (i = 0; i < p->size; ++i)
                fprintf(stdout, " %f", p->vec[i]);
            break;
        case GLPS_NODE_STR:
            fprintf(stdout, "%s: \"%s\" ", p->name, p->str);
            break;
        case GLPS_NODE_FUNC1:
            fprintf(stdout, "F: %s(x) ", p->name);
            break;
        case GLPS_NODE_FUNC2:
            fprintf(stdout, "F: %s(x,y)", p->name);
            break;
        case GLPS_NODE_FUNC3:
            fprintf(stdout, "F: %s(x,y,z)", p->name);
            break;
        case GLPS_NODE_TAG:
            fprintf(stdout, "");
            break;
        case GLPS_ACTION:
            fprintf(stdout, "%s ", p->name);
            ++depth;
            symb_print(p->property);
            --depth;
            break;
        case GLPS_ELEMENT:
            fprintf(stdout, "%s: %s ", p->name, p->str);
            ++depth;
            symb_print(p->property);
            --depth;
            break;
        case GLPS_LINE:
            et_print(p->line);
            fprintf(stdout, ":%s", p->name);
            break;
        default:
            fprintf(stdout, "UNKNOW type: %d %d", p->nodetype, GLPS_ELEMENT);
        }
        
        if (depth == 0) fprintf(stdout, "\n");

        p = p->next;
    }
}

GLPS_SYMB* symb_new(char* name)
{
    /* DBG_MSG("new %s\n", name); */
    GLPS_SYMB* p = (GLPS_SYMB*) malloc(sizeof(GLPS_SYMB));
    p->name = strdup(name);
    p->nodetype = GLPS_NODE_UNDEF;
    p->size     = 0;
    p->value    = 0;
    p->funcptr  = NULL;
    p->next     = NULL;
    p->property = NULL;
}

/**
 * assign the node contents, not the position in a linked list
 */
void symb_assign_node(GLPS_SYMB* dst, const GLPS_SYMB* src)
{
    /* do not touch dst->next */
    int i = 0;
    dst->nodetype = src->nodetype;
    dst->size     = src->size;
    switch(src->nodetype) {
    case GLPS_NODE_REAL:
        dst->value    = src->value;
        break;
    case GLPS_NODE_STR:
        dst->str = strdup(src->str);
        break;
    case GLPS_NODE_VEC:
        dst->vec = (double*) malloc(src->size * sizeof(double));
        for (i = 0; i < src->size; ++i) dst->vec[i] = src->vec[i];
        break;
    case GLPS_NODE_FUNC1:
    case GLPS_NODE_FUNC2:
    case GLPS_NODE_FUNC3:
        dst->funcptr = src->funcptr;
        break;
    default:
        fprintf(stderr, "=glps=: Unknown data type %d when copying %s\n",
                src->nodetype, src->name);
    }

    if (src->property) {
        fprintf(stderr, "=glps=: Did not copy properties\n");
        dst->property = NULL;
    }
}

GLPS_SYMB* symb_copy_node(const GLPS_SYMB* src)
{
    GLPS_SYMB* p = symb_new(src->name);
    symb_assign_node(p, src);
    return p;
}

int symb_free_list(GLPS_SYMB* p)
{
    int n = 0;
    GLPS_SYMB* head = p;
    while (p) {
        /* if (p->name) DBG_MSG("freeing: %s\n", p->name); */
        head = p->next;
        switch(p->nodetype) {
        case GLPS_ELEMENT:
            free(p->str); break;
        case GLPS_NODE_STR:
            free(p->str); break;
        case GLPS_NODE_VEC:
            if (p->size > 0) free(p->vec);
            break;
        case GLPS_LINE:
            et_delete(p->line);
            break;
        case GLPS_ACTION:
        case GLPS_NODE_TAG:
        case GLPS_NODE_REAL:
        case GLPS_NODE_FUNC1:
        case GLPS_NODE_FUNC2:
        case GLPS_NODE_FUNC3:
            break;
        case GLPS_NODE_UNDEF:
            break;
        default:
            fprintf(stderr, "free an unknown type in node: %s code=%d\n",
                    p->name, p->nodetype);
        }

        if (p->property) symb_free_list(p->property);

        if (p->name) free(p->name);
        free(p);
        ++n;
        p = head;
    }
    return n;
}


GLPS_SYMB* symb_find_append(char* name)
{
    if (symbtab == NULL) {
        symbtab = symb_new(name);
    }

    GLPS_SYMB* p = symbtab;
    while (p) {
        if (str_case_cmp(p->name, name) == 0) {
            return p;
        } else if (p->next == NULL) {
            p->next = symb_new(name);
            p = p->next;
            break;
        } else {
            p = p->next;
        }
    }

    return p;
}


/*! \brief Add function of single variable into expressions.
 *
 * \param Function name, such as "sin", "cos", "tan".
 * \param pointer to the function, can be sin, cos, tan, and also a user
 * defined function.
 *
 * \brief add function support, sin, exp, cos, .....
 */
void symb_addfunc(char* name, double (*func) (double))
{
    GLPS_SYMB *sp = symb_find_append(name);
    typedef double (*func1)();

    sp->nodetype = GLPS_NODE_FUNC1;
    sp->funcptr = (func1)func;
}


/*! \brief Add predefined constant
 *
 * \param Constant variable name, such as "PI".
 * \param The definition
 *
 */
void symb_addval(char* name, double val)
{
    GLPS_SYMB *sp = symb_find_append(name);
    sp->nodetype = GLPS_NODE_REAL;
    sp->value = val;
    /* symb_print(sp); */
}


GLPS_SYMB* symb_append_node(GLPS_SYMB* head, GLPS_SYMB* p)
{
    GLPS_SYMB* tail = head;
    while(tail->next) {
        tail = tail->next;
    }
    
    tail->next = p;
    return tail;
}

/*! \brief Get the property table (linked list) of an element. If it does not
 * exist create one.
 *
 * \param elem pointer to an element.
 * \param propty property name.
 * \brief Get the property of an element.
 */
GLPS_SYMB *get_statement_property(
    GLPS_SYMB *elem,
    const char* propty, int append)
{
}



/*! \brief get the path and root file name 
 *
 * allocated new mem, need free by caller.
 */
char *get_path_root(char *f)
{
    size_t n = 0; //
    char *r;
    char* s = f + strlen(f);
    while(s > f && *s != '.') {
        /* fprintf(stdout, "%c", *s); */
        --s;
    }
    /* fprintf(stdout, "\n"); */

    if (s == f) return NULL;
    n = s - f;

    r = (char*) malloc(sizeof(char)*(n+1));
    strncpy(r, f, n);
    r[n] = '\0';

    return r;
}



int symb_get_double(char* name, double* val)
{
    GLPS_SYMB* p = symb_lookup(symbtab, name);
    
    if (!p) {
        return GLPS_ENOVAR;
    } else if (p->nodetype != GLPS_NODE_REAL) {
        return GLPS_ENOTREAL;
    } else {
        *val = p->value;
        return GLPS_SUCCESS;
    }
}

int symb_get_vector(char* name, double** val, size_t *n)
{
    GLPS_SYMB* p = symb_lookup(symbtab, name);
    
    if (!p) {
        return GLPS_ENOVAR;
    } else if (p->nodetype != GLPS_NODE_VEC) {
        return GLPS_ENOTVEC;
    } else {
        *val = p->vec;
        *n = p->size;
        return GLPS_SUCCESS;
    }
}

int symb_get_string(char* name, char** val)
{
    GLPS_SYMB* p = symb_lookup(symbtab, name);
    
    if (!p) {
        return GLPS_ENOVAR;
    } else if (p->nodetype != GLPS_NODE_STR) {
        return GLPS_ENOTSTR;
    } else {
        *val = p->str;
        return GLPS_SUCCESS;
    }
}


/** use symb_lookup for element/beamline */
elemtree *beamline_defined(const char* name)
{
    GLPS_SYMB* p = symb_lookup(symbtab, name);
    if (!p || p->nodetype != GLPS_LINE ) return NULL;
    /* DBG_MSG("found beam line: %p %s\n", p, p->name); */
    return p->line;
}
