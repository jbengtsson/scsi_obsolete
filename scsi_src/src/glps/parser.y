/* -*- Mode: c; tab-width:4; indent-tabs-mode: nil; c-basic-offset:4 -*-  */
%{
 /*
  * Copyright (C) 2008-2011 Lingyun Yang.
  *
  * This software may be used and distributed according to the terms
  * of the GNU General Public License, incorporated herein by reference.
  *
  * For more information please contact lyyang@lbl.gov
  *
  */

 /* Author: Lingyun Yang, lyyang@lbl.gov */
 /* date: 2011-05-20 09:40 */

 #include <stdio.h>
 #include <string.h>
 #include <math.h>
 #include <ctype.h>

 #include "glps.h"
 #include "latio.h"
 #include "elemtree.h"
 #include "glpserror.h"

 #ifndef M_PI
 #define M_PI 3.14159265358979323846
 #endif

 #define XML 1

    /* int yydebug=1; */

 /*! global symbol table (linked list) */
 /*  static SP_SYMB_LST* symbtab = NULL; */
 GLPS_SYMB* symbtab = NULL;

 extern unsigned int lineno;
 extern char linebuf[];
 extern char strbuf[];
 extern char* glps_current_latfile;

 extern int lat_file_stack_ptr;

 /*! lattice title */
 static char* lat_title = NULL;

 static int echo = 1;
 /*! Flat file */
 /* static FILE* pflat; */

 /* ["key_0", "elements_0", "key_1", "elements_1", ... ] */
 char** beam_lines; 
 /* a linked list of beamlines(elemtree type) */
 bltree *bl_list, *bl_list_tail;
 bltree* beamline_defined(char* beamline);

 /*! size of beam_lines array, twice of the number of beam lines */
 int n_beam_lines = 0;

 /* save the error message */
 char glps_err_msg[20480];
 int glps_err_code = GLPS_SUCCESS;
 FILE * glps_err_stream = NULL ;

 /* if yywrap returns 0, the scanner continues scanning, while if it returns 1,
   the scanner returns a zero token to report the end of file. If yywrap()
   returns 0 to indicate that there is more input, it needs first to adjust
   yyin to point to a new file, probably using fopen()
 */
 extern FILE* yyin;
 int yywrap(void) { 
     return 1; 
 };
 int yylex();
 int yyerror(char* s);


 bltree *beamline_defined(char* name);

%}

%locations

%union {
    double  dval;   /* plain number */
    char*    str;   /* string */
    GLPS_SYMB *symp;     /* variable */
    struct elemtree_el *bl;
}

%token <str> STRING ELEMENT ACTION
%token <dval> REAL
%token <dval> NUMBER
%token <symp> ELEM_PROP ACT_PROP
%token <symp> ID
%token <dval> FUNC

%token <str> BL

%token SET TITLE SHOW INCLUDE
%token BEND BPM CAVITY CORR DRIFT MARKER MULTIPOLE INV LINE QUAD SEXT WIGGLER 

%left ','
%left '-' '+'
%left '*' '/'
%nonassoc UMINUS

%type <symp> expression
%type <symp> expression_list
%type <symp> stmt_property
%type <symp> stmt_property_list
/* %type <str> beamline */
%type <bl> beamline

%%

statement_list: statement 
        |       statement_list statement
        ;

statement: ID '=' expression ';' {
    /* DBG_MSG("line: %d,  assign %s %s %f\n", yyloc.first_line,
       $1->name, $3->name, $3->value); */
    /* symb_print(symbtab); */
    GLPS_SYMB* p = symb_lookup(symbtab, $1->name);
    if (p == NULL) {
        if (glps_err_code) {
            GLPS_ERROR_MSG(glps_err_code, "=glps=: %s\n", glps_err_msg);
            YYABORT;
        }
        /* DBG_MSG("append var %s\n", $1->name); */
        symb_assign_node($1, $3);
        symb_free_list($3);
        symb_append_node(symbtab, $1);
    } else if (p->nodetype == $3->nodetype) {
        symb_assign_node(p, $3);
        /* $1 has UNDEF type, if in 'elem.field = 0.5' */
        $1->nodetype = p->nodetype;
        symb_free_list($1);
        symb_free_list($3);
    } else {
        symb_free_list($3);
        GLPS_ERROR_MSG(GLPS_EASSIGN, "incompatible type\n");
        YYABORT;
    }
    /* symb_print(symbtab); */
 }
| ID '=' STRING ';' {
    GLPS_SYMB* p = symb_lookup(symbtab, $1->name);
    if (p == NULL) {
        /* DBG_MSG("append string %s\n", $1->name); */
        $1->nodetype = GLPS_NODE_STR;
        $1->str = $3;
        symb_append_node(symbtab, $1);
    } else if (p->nodetype == GLPS_NODE_STR) {
        if (p->str) free(p->str);
        p->str = $3;
        p->size = strlen($3);
        p->nodetype = GLPS_NODE_STR;
    } else {
        free($3);
        GLPS_ERROR_MSG(GLPS_EGRAMMAR, 
                       "can not assign string to non-string '%s'\n", $1->name);
        symb_free_list($1);
        YYABORT;
    }
   }
| ID '=' '(' expression_list ')' ';' { 
    GLPS_SYMB* p = symb_lookup(symbtab, $1->name);

    if (p == NULL) {
        if (glps_err_code) {
            GLPS_ERROR_MSG(glps_err_code, "=glps=: %s\n\n", glps_err_msg);
            YYABORT;
        }
        symb_assign_node($1, $4);
        symb_free_list($4);

        symb_append_node(symbtab, $1);
    } else {
        /* symb_print(symbtab); */
        /* re-assign to new value */
        if (p->nodetype != GLPS_NODE_VEC) {
            GLPS_ERROR_MSG(GLPS_ENOTVEC, "'%' is not a vector\n", p->name);
            YYABORT;
        }
        free(p->vec);
        p->size = $4->size;
        p->vec = $4->vec;
        $4->size = 0;
        symb_free_list($4);
    }
    /* symb_print(symbtab); */
   }
| ACTION ',' ID ';' {
    /* GLPS_ERROR_MSG(GLPS_EGRAMMAR, "'act, id;' not supported yet\n"); */
    GLPS_SYMB* p = symb_lookup(symbtab, $3->name);
    if (!p) {
        GLPS_ERROR_MSG(GLPS_ENOVAR, "undefined variable '%s'\n", $3->name);
        symb_free_list($3);
        YYABORT;
    } else {
        GLPS_SYMB* act = symb_new($1);
        symb_assign_node($3, p);
        act->nodetype = GLPS_ACTION;
        act->property = $3;
        symb_append_node(symbtab, act);
        free($1);
    }
 }
| ACTION ',' STRING ';' {
    /* GLPS_ERROR_MSG(GLPS_EGRAMMAR, "'act, str;' not supported yet\n"); */
    GLPS_SYMB* act = symb_new($1);
    GLPS_SYMB* prpt = symb_new("string");
    act->property = prpt;
    act->nodetype = GLPS_ACTION;
    prpt->nodetype = GLPS_NODE_STR;
    prpt->str = $3;
    symb_append_node(symbtab, act);
    free($1);
 }
| ACTION ',' stmt_property_list ';' {
    /* GLPS_ERROR_MSG(GLPS_EGRAMMAR, "'act, str;' not supported yet\n"); */
    GLPS_SYMB* act = symb_new($1);
    act->property = $3;
    act->nodetype = GLPS_ACTION;
    symb_append_node(symbtab, act);
    free($1);
  }
|  ELEMENT ':' ID ';' {
    /* DBG_MSG("A simple element\n"); */
    $3->nodetype = GLPS_ELEMENT;
    $3->str = $3->name;
    $3->name = $1;
    symb_append_node(symbtab, $3);
   }
| ELEMENT ':' ID ',' stmt_property_list ';' {
    /* DBG_MSG("Complicated element\n"); */
    /* symb_print($3); */
    /* symb_print($5); */
    $3->nodetype = GLPS_ELEMENT;
    $3->str = $3->name;
    $3->name = $1;
    $3->property = $5;
    symb_append_node(symbtab, $3);
  }
|  ID ':' LINE '=' '(' beamline ')' ';'  {
    /* DBG_MSG("Now the final beam line part !!\n"); */

    int nleaf = 0, nnode = 0;
    et_count($6, &nleaf, &nnode);
    $1->nodetype = GLPS_LINE;
    $1->line = et_new_node($1->name);
    $1->line->left = $6;
    $1->line->right = NULL;
    symb_append_node(symbtab, $1);
   }
;

expression_list: expression ',' expression {
    /* new a symp */
    GLPS_SYMB *p = symb_new("expression_list");
    p->nodetype = GLPS_NODE_VEC;
    p->size = 2;
    p->vec = (double*) malloc(2*sizeof(double));
    p->vec[0] = $1->value;
    p->vec[1] = $3->value;

    symb_free_list($1);
    symb_free_list($3);
    $$ = p;
 }
| expression_list ',' expression {
    int i = 0;
    double *v = $1->vec;
    size_t n = $1->size;
    $1->vec = (double*) malloc((n + 1) * sizeof(double));
    for ( i = 0; i < n; ++i) $1->vec[i] = v[i];
    $1->vec[n] = $3->value;
    free(v);
    symb_free_list($3);
    $1->size = n + 1;
    $$ = $1;
}
;

expression: expression '+' expression {
    if ($1->nodetype != $3->nodetype) {
        GLPS_ERROR_MSG(GLPS_EGRAMMAR, "can not do '+' for different data type");
    } else {
        GLPS_SYMB *p = symb_new("expression");
        p->value = $1->value + $3->value;
        p->nodetype = GLPS_NODE_REAL;
        $$ = p;
    }
    symb_free_list($1);
    symb_free_list($3);
 }
|  expression '-' expression { 
    if ($1->nodetype != $3->nodetype) {
        GLPS_ERROR_MSG(GLPS_EGRAMMAR, "can not do '-' for different data type");
    } else {
        GLPS_SYMB *p = symb_new("expression");
        p->value = $1->value - $3->value;
        p->nodetype = GLPS_NODE_REAL;
        $$ = p;
    }
    symb_free_list($1);
    symb_free_list($3);
}
|  expression '*' expression {
    if ($1->nodetype != $3->nodetype) {
        GLPS_ERROR_MSG(GLPS_EGRAMMAR, "can not do + for different data type");
    } else {
        GLPS_SYMB *p = symb_new("expression");
        p->value = $1->value * $3->value;
        p->nodetype = GLPS_NODE_REAL;
        $$ = p;
    }
    symb_free_list($1);
    symb_free_list($3);
   }
|  expression '/' expression {
    if ($1->nodetype != $3->nodetype) {
        GLPS_ERROR_MSG(GLPS_EGRAMMAR, "can not do + for different data type");
    } else {
        GLPS_SYMB *p = symb_new("expression");
        p->value = $1->value / $3->value;
        p->nodetype = GLPS_NODE_REAL;
        $$ = p;
    }
    symb_free_list($1);
    symb_free_list($3);
    }
|  '-' expression %prec UMINUS { 
    if ($2->nodetype = GLPS_NODE_REAL) $2->value = -($2->value);
    $$ = $2; 
   }
|  '(' expression ')' { $$ = $2; }
|  NUMBER { 
    char s[128];
    sprintf(s, "number %f", $1);
    GLPS_SYMB *p = symb_new(s);
    p->nodetype = GLPS_NODE_REAL;
    p->value = $1;
    $$ = p;
   }
|  REAL { 
    char s[128];
    sprintf(s, "number %f", $1);
    GLPS_SYMB *p = symb_new(s);
    p->nodetype = GLPS_NODE_REAL;
    p->value = $1;
    $$ = p;
   }
|  ID { 
    /* COPY instead of lookup and return. This ID will appears at the RHS
       and will be freed */
    GLPS_SYMB* p = symb_lookup(symbtab, $1->name);
    if (!p) {
        GLPS_ERROR_MSG(GLPS_ENOVAR, "undefined variable %s\n", $1->name);
        /* could have memleak when error */
        /* symb_print(symbtab); */
        YYABORT;
    } else {
        symb_assign_node($1, p);
        $$ = $1;
    }
   }
|  ID '(' expression ')' { /* for function like sin(pi) */
    GLPS_SYMB* p = symb_lookup(symbtab, $1->name);
    if (!p) {
        GLPS_ERROR_MSG(GLPS_ENOFUNC, "undefined function '%s'\n", $1->name);
        YYABORT;
    }
    /* DBG_MSG("new function %s\n", $1->name); */
    if ( p->funcptr ) {
        $1->nodetype = GLPS_NODE_REAL;
        $1->value = (p->funcptr)( $3->value );
        symb_free_list($3);
        $$ = $1;
    } else {
        GLPS_ERROR_MSG(GLPS_ECALLABLE, "'%s' is not callable(not a function)\n",
                   $1->name);
        $$ = NULL;
    }
   }
;

stmt_property_list: stmt_property {
    $$ = $1;
 }
| stmt_property_list ',' stmt_property {
    symb_append_node($1, $3);
    $$ = $1;
  }
;

stmt_property: ELEM_PROP '=' expression {
    symb_assign_node($1, $3);
    symb_free_list($3);
    $$ = $1;
 }
| ELEM_PROP '=' STRING {
    $1->nodetype = GLPS_NODE_STR;
    $1->str = $3;
    $$ = $1;
  }
| ELEM_PROP '=' '(' expression_list ')' {
    symb_assign_node($1, $4);
    symb_free_list($4);
    $$ = $1;
  }
| ACT_PROP '=' expression {
    symb_assign_node($1, $3);
    symb_free_list($3);
    $$ = $1;
  }
| ACT_PROP '=' STRING {
    /* assume it is initial */
    $1->str = $3;
    $1->nodetype = GLPS_NODE_STR;
    $$ = $1;
  }
| ACT_PROP '=' '(' expression_list ')' {
    symb_assign_node($1, $4);
    symb_free_list($4);
    $$ = $1;
  }
;

beamline: BL { 
    elemtree *bl = beamline_defined($1);
    if ( bl == NULL) {
        /* a new element */
        $$ = et_new_node($1);
    } else {
        /* copy tree bl->bl, use NULL as its parent, no reverse order */
        $$ = et_copy_tree(bl, NULL, 1);
    }

    /* DBG_MSG("beamline: ( %s )\n", $$->val); */
    
    free($1);
 }
|  beamline ',' beamline {
    $$ = et_combine($1, $3);
   }
|  NUMBER '*' '(' beamline ')' {
    int i = 0, inv = 1;
    int ndup = (int)($1);
    elemtree *et = NULL, *tmp;
    if (ndup == 0) {
        et = NULL;
    } else {
        inv = ndup < 0 ? -1 : 1;
        for (i = 0; i < abs(ndup); ++i) {
            tmp = et_copy_tree($4, NULL, inv);
            et = et_combine(et, tmp);
        }
    }
    et_delete($4);
    /* fprintf(stderr, "---+-> %d*%s\n", ndup, $4->val); */
    $$ = et;
   }
|  NUMBER '*' BL {
    /* could be a beamline name, could be an element */
    /* DBG_MSG("expand beamline: %f %s\n", $1, $3); */
    int ndup = abs($1), i = 0;
    int inv = $1 < 0 ? -1 : 1;

    elemtree *ret = NULL, *tmp;
    elemtree *bl = beamline_defined($3);
    if (bl == NULL) {
        for (i = 0; i < ndup; ++i) {
            tmp = et_new_node($3);
            et_push_back(&ret, tmp);
        }
    } else {
        /* DBG_MSG("expand beamline: %f %s inv=%d\n", $1, $3, inv); */
        for (i = 0; i < ndup; ++i) {
            ret = et_combine(ret, et_copy_tree(bl, NULL, inv));
        }
    }

    free($3);
    ret->factor = 1;
    /* fprintf(stderr, "N*BL: %s\n", $3); */
    $$ = ret;
   }
|  INV BL { 
    elemtree *ret = NULL;
    elemtree* bl = beamline_defined($2);
    if (bl == NULL) {
        ret = et_new_node($2);
        ret->factor = -1;
    } else {
        ret = et_copy_tree(bl, NULL, -1);
    }
    /* fprintf(stderr, "Reverse Beam Line: %s\n", $2); */
    free($2);
    $$ = ret;
 }
|  INV '(' beamline ')' {
    elemtree* ret = et_copy_tree($3, NULL, -1);
    /* ret->factor = -1; */
    et_delete($3);
    $$ = ret;
   }
| '(' beamline ')' {
    $$ = $2;
  }
;
%%





/*! \brief Error report with line number.
 *
 * If the error (syntax error) happens at the first line, the position (line
 * number) is not given correctly.
 */
int yyerror(char* s)
{
    fprintf(stderr, "Error: %s\nLine %u, %s\n  %s\n",
            glps_current_latfile, lineno, s, linebuf);
    fprintf(stderr, "ERROR: line: %d - %d, column: %d - %d\n",
            yylloc.first_line, yylloc.last_line, yylloc.first_column,
            yylloc.last_column);
    sprintf(glps_err_msg, "ERROR: %s\n  at line: %d - %d, column: %d - %d\n",
            s, yylloc.first_line, yylloc.last_line, yylloc.first_column,
            yylloc.last_column);
    /* GLPS_ERROR_MSG("unknown error", GLPS_EERR) ; */
}

int lyyerror(YYLTYPE t, char* s)
{
    fprintf(stderr, "Error: %s\nLine %u, %s\n  %s\n",
            glps_current_latfile, lineno, s, linebuf);
    fprintf(stderr, "ERROR: line: %d - %d, column: %d - %d\n",
            t.first_line, t.last_line, t.first_column, t.last_column);
    /* GLPS_ERROR_MSG("unknown error", GLPS_EERR) ; */
}


/*! \brief Add a property list to a element.
 *
 * append to the end.
 */

int parse_lattice(const char* f)
{
    /* add mathematical function support */
    extern double sqrt(), exp(), log(), atan(), sin(), cos(), fabs();

    /* initialize the random generator */
    srand(0);

    /* initialize the stack of "include" */
    lat_file_stack_ptr = 0;
    glps_err_code = 0;

    if (!(yyin = fopen(f, "r"))) {
        // fprintf(stderr, "ERROR: can  not open file %s\n", f);
        // return -1;
        GLPS_ERROR_MSG(GLPS_ENOLAT, "can not open lattice file '%s'", f) ;
    }
    glps_current_latfile = strdup(f);

    symb_addval("PI", M_PI);
    symb_addval("true", 1.0);
    symb_addval("false", 0.0);

    symb_addfunc("sqrt", sqrt);
    symb_addfunc("exp", exp);
    symb_addfunc("log", log);
    symb_addfunc("log10", log10);

    symb_addfunc("sin", sin);
    symb_addfunc("cos", cos);
    symb_addfunc("tan", tan);
    symb_addfunc("asin", asin);
    symb_addfunc("acos", acos);
    symb_addfunc("atan", atan);
    symb_addfunc("cosh", cosh);
    symb_addfunc("sinh", sinh);
    symb_addfunc("tanh", tanh);

    symb_addfunc("abs", fabs);

    /* yypush_buffer_state(yycreate_buffer(yyin, 10240)); */
    int errcode = 0;
    while(!feof(yyin)) {
        /* the detailed code is in glps_err_code */
        errcode = yyparse();
        /* fprintf(stdout, "Returned from parser!\n"); */
    }
    /* DBG_MSG("closing file %s code:%d\n", f, glps_err_code); */
    fclose(yyin);
    free(glps_current_latfile);
    glps_current_latfile = NULL;

 #ifdef DEBUG
    fprintf(stdout, "\n#------------------------------------------------#\n");
    fprintf(stdout, "# DEBUG: printing all statements\n");
    fprintf(stdout, "#------------------------------------------------#\n\n");
    if (!errcode) symb_print(symbtab);
 #endif

    return glps_err_code;
}


int free_lattice()
{
    /* Free all the memory for lattice parsing */
    symb_free_list(symbtab);
    /* symb_print(symbtab); */
    symbtab = NULL;
    lat_file_stack_ptr = 0;
    lineno = 0;

    /* clear the error messages */
    glps_err_code = GLPS_SUCCESS;
    glps_err_msg[0] = '\0';
}


void
glps_error (int glps_errno, const char * reason)
{
    if (glps_err_stream == NULL) {
        fflush(stdout);
        fprintf(stderr, reason);
    } else {
        fflush(glps_err_stream);
        fprintf(glps_err_stream, reason);
        fflush(glps_err_stream);
    }
}
