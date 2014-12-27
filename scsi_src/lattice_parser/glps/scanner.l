/* -*- Mode: c; tab-width:4; indent-tabs-mode: nil; c-basic-offset:4 -*-  */
%{
/* 
 *
 * Copyright (C) 2008 Lingyun Yang.
 *
 * This software may be used and distributed according to the terms
 * of the GNU General Public License, incorporated herein by reference.
 *
 * For more information please contact lingyun.yang@gmail.com
 *
 */


/* Author: L. Yang, lingyun.yang@gmail.com */
/* date: 2010-09-01 14:37 */

#include "glps.h"
#include "parser.h"
#include "glpserror.h"
#include <math.h>
#include <string.h>
#include <ctype.h>

#define MAX_LINE_BUF_SIZE 1024

/* #define DEBUG 1 */

char* element_name;
unsigned int lineno = 1;

char linebuf[MAX_LINE_BUF_SIZE];
char strbuf[MAX_LINE_BUF_SIZE];

extern char** beam_lines;
extern int yydebug;
extern GLPS_SYMB* symbtab;

/* open multiple lattice file, push/pop the stack. */
#define MAX_INCLUDE_DEPTH 10
int lat_file_stack_ptr = 0;
static YY_BUFFER_STATE lat_file_stack[MAX_INCLUDE_DEPTH]; 
static int lat_file_lineno[MAX_INCLUDE_DEPTH];
static char* glps_lat_file[MAX_INCLUDE_DEPTH];

char* glps_current_latfile;

#define MAX_ERRMSG 1024
static char errmsg[MAX_ERRMSG];

 /* processing property statement in an definition 
  * e.g. L = 0.1, property=(1,2,3,4,5); 
  * every time, when meet an identifier, we record its name,
  * in case we need it on the right hand side. we have to keep its value later.
  * i.e. vec/vec_size/value 
  */

int yycolumn = 1;
#define YY_USER_ACTION yylloc.first_line = yylloc.last_line = yylineno; \
    yylloc.first_column = yycolumn; yylloc.last_column = yycolumn+yyleng-1; \
    yycolumn += yyleng;
%}

%option yylineno
/* %option debug */
/* %option reentrant */

  /* definitions */

  /* NUM   -?[0-9]+ */
  /* REAL  -?(([0-9]+)|(([0-9]*\.)?[0-9]+)([eE][-+]?[0-9]+)?)  */
NUM   [0-9]+
REAL  (([0-9]+)|(([0-9]*\.)?[0-9]+)([eE][-+]?[0-9]+)?) 
WS    [ \t\n]*
IDENTIFIER  [A-Za-z][A-Za-z0-9_\.]*

  /* x exclusive start-condition: */
%x sc_comment sc_str sc_beamline sc_element sc_action sc_assign

%%

  /* --------------------------------------- */

  /* rules */
<INITIAL><<EOF>>  { 
 #ifdef DEBUG
    fprintf(stderr, "--- Reach the End of Lattice without an \"END\" statement. --\n"); 
 #endif

    if ( --lat_file_stack_ptr < 0) {
        /* DBG_MSG("Thats the end %d\n", lat_file_stack_ptr); */
        /* yy_delete_buffer(YY_CURRENT_BUFFER); */
        yy_flush_buffer(YY_CURRENT_BUFFER);
        yyterminate();
    }else {
     #ifdef DEBUG
        DBG_MSG("==glps== [INFO] Switching to buffer ...\n");
     #endif
        yy_delete_buffer(YY_CURRENT_BUFFER);
        yy_switch_to_buffer(lat_file_stack[lat_file_stack_ptr]);
        lineno = lat_file_lineno[lat_file_stack_ptr];
    }
}

<sc_element,sc_action,sc_beamline><<EOF>>  {
    sprintf(errmsg, "statement is incomplete at the end of file");
    DBG_MSG("statement is incomplete at the end of file\n");
    GLPS_LAT_ERROR(errmsg, glps_current_latfile,
                   lineno, linebuf, GLPS_EEOF);
    /*
    fprintf(stderr, "Reach the end of lattice file,"
            "but the statement is not finished with ;\n");
    return 0;
    */
     }

  /* Syntax error at first line will also be copied into linebuf */
<*>\n.*  { 
    ++lineno;
    yycolumn = 1;
    strncpy(linebuf, yytext+1, MAX_LINE_BUF_SIZE-1);
    linebuf[MAX_LINE_BUF_SIZE-1] = '\0';
    yyless(1); 
 #ifdef DEBUG
    fprintf(stdout, "%d: %s\n", lineno, linebuf); 
 #endif
   }

  /* Comment */ 
  /*
  "{"                    { BEGIN(sc_comment); }
  <sc_comment>[^"}"\n]+  ;
  <sc_comment>\n         {++lineno;}
  <sc_comment>"}"+       { 
      BEGIN(INITIAL); 
      fprintf(stderr, "--- COMMENT ---\n");
   }
   */

#.*\n  { ++lineno; yycolumn = 1; /* fprintf(stderr, "\n---- Comment ----\n"); */} 

<*>\"[^\"\n]*\"  { 
    /* String, doesn't span > 1 line */
    /* TODO: remove this constraint */
    /* assuming "string", here remove the " */
    yylval.str = (char* )malloc(strlen(yytext));
    strncpy(yylval.str, yytext+1, strlen(yytext)-2);
    yylval.str[strlen(yytext)-2] = '\0';
    return STRING; 
   }

<*>\"[^\"\n]*$ {
    fprintf(stderr, "Unterminated string\n");
    yylval.str = strdup(yytext);
    return STRING;
   }

<*>{REAL} { yylval.dval = atof(yytext); return NUMBER; }
<*>{NUM}+ { yylval.dval = atof(yytext); return NUMBER; }

include[ \t\n]*","[ \t\n]*[^;\t\n]+";" {
    /* include the lattice input */
    char *head, *tail, *inc, *incfile;
    size_t sz;
    if (lat_file_stack_ptr >= MAX_INCLUDE_DEPTH) {
        sprintf(errmsg, "lattice file nested too deep (>%d)", MAX_INCLUDE_DEPTH);
        GLPS_LAT_ERROR(errmsg, glps_current_latfile,
                           lineno, linebuf, GLPS_ENOELEMT);
        fprintf(stderr, "Includes nested too deeply\n");
        exit(-1);
    }
 #ifdef DEBUG
    fprintf(stderr, "#[INFO] Include file: %s\n", yytext);
 #endif

    head = strchr(yytext, '"');
    tail = strrchr(yytext, '"');
    if (head && tail) sz = tail - head; /* strspn(head, tail); */
    else {
        sprintf(errmsg, "file name should be quoted, e.g. include, \"lattice.lat\"");
        GLPS_LAT_ERROR(errmsg, glps_current_latfile,
                   lineno, linebuf, GLPS_ENOELEMT);
        fprintf(stderr, "[ERROR]: %s is not a quoted file name\n", yytext);
    }

    inc = (char*) malloc((sz+2)*sizeof(char));
    incfile = inc;
    ++head;
    while(head != tail){
        *inc = *head;
        ++inc;
        ++head;
    }
    *inc = '\0';

 #ifdef DEBUG
    fprintf(stderr, "File1: %s\n", incfile);
 #endif

    glps_lat_file[lat_file_stack_ptr] = glps_current_latfile;
    lat_file_lineno[lat_file_stack_ptr] = lineno;
    lat_file_stack[lat_file_stack_ptr++] = YY_CURRENT_BUFFER;

    yyin = fopen(incfile, "r");

    if (!yyin) {
        sprintf(errmsg, "can not open lattice file: \"%s\"", incfile);
        GLPS_LAT_ERROR(errmsg, glps_current_latfile,
                       lineno, linebuf, GLPS_EFOPEN);
    }
    
    glps_current_latfile = strdup(incfile);
    yy_switch_to_buffer(yy_create_buffer(yyin, YY_BUF_SIZE));
    lineno = 0;

    free(incfile);
   }

<INITIAL>end{WS}";" { /* logical end of file */
    #ifdef DEBUG
        fprintf(stderr, "-- End of Lattice --\n"); 
    #endif

    if ( lat_file_stack_ptr == 0) {
        /* DBG_MSG("LOGICAL END of lattice file\n"); */
        /* DBG_MSG("T END of lattice file\n"); */
        yy_flush_buffer(YY_CURRENT_BUFFER);
        /* yy_delete_buffer(YY_CURRENT_BUFFER); */
        yyterminate();
    }else {
        --lat_file_stack_ptr;
     #ifdef DEBUG
        fprintf(stderr, "#[INFO] Switching... to buffer %d\n", lat_file_stack_ptr);
     #endif

        yy_delete_buffer(YY_CURRENT_BUFFER);
        yy_switch_to_buffer(lat_file_stack[lat_file_stack_ptr]);
        lineno = lat_file_lineno[lat_file_stack_ptr];
        free(glps_current_latfile);
        glps_current_latfile = glps_lat_file[lat_file_stack_ptr];
    }
    } 

<INITIAL>{IDENTIFIER}/{WS}":"{WS}"line"{WS}"=" {
    GLPS_SYMB* sp;
    BEGIN(sc_beamline); 
    /* fprintf(stderr, "\n----> Line: %s\n", yytext); */
    /* NEW_PROPERTY_SYMBOL(yytext); */
    sp = symb_new(yytext);
    yylval.symp = sp;
    
    return ID;}


<INITIAL>{IDENTIFIER}/{WS}: {
    BEGIN(sc_element);
    /* fprintf(stderr, "\n---> element: %s\n", yytext); */

    yylval.str = strdup(yytext);
    /* yylval.str = strdup("hello"); */

    return ELEMENT;
 }

<INITIAL>{IDENTIFIER}/{WS}=  {
    BEGIN(sc_assign);
    yyless(0);
}
<INITIAL>{IDENTIFIER}/{WS},  { 
    BEGIN(sc_action); 
    /* fprintf(stderr, "Action: %s\n", yytext); */
    yylval.str = strdup(yytext);
    return ACTION; 
 }

<INITIAL>{IDENTIFIER}[ \t\n]+{IDENTIFIER}/{WS},  { 
    BEGIN(sc_action); 
    /* fprintf(stderr, "Action: %s\n", yytext); */
    yylval.str = strdup(yytext);
    return ACTION; 
 }

<sc_beamline>LINE/{WS}=   { 
    /* fprintf(stderr, "\n---> beamline: \n", yytext); */
    return LINE; }

<sc_action>{IDENTIFIER}/{WS}= {
    /* convert_key_case(yytext); */
 #ifdef DEBUG
    fprintf(stderr, "Action: %s\n", yytext);
 #endif  

    yylval.symp = symb_new(yytext);
    return ACT_PROP; }

<sc_element>{IDENTIFIER}/{WS}= { 
    /* general element */ 
 #ifdef DEBUG
    fprintf(stderr, "Magnet Property: %s\n", yytext);
 #endif

    yylval.symp = symb_new(yytext);
    return ELEM_PROP; 
}

<sc_beamline>{IDENTIFIER} { 
    /* BL: Line=(ID1, ID2, ...) */
    yylval.str = strdup(yytext);
    
    return BL; 
 }

<*>{IDENTIFIER} {
    /* always a new copy */
    yylval.symp = symb_new(yytext);
    return ID; 
}

<sc_beamline>-?[0-9]+ { yylval.dval = atoi(yytext); return NUMBER;}
<sc_beamline>"-" { return INV; }

<*>"+"|"-"|"*"|"/" { return yytext[0]; }

<*>":"|","|"="|"("|")" { return yytext[0]; }

<*>";" { BEGIN(INITIAL); return yytext[0]; }

<*>[ \t\r] ;

<*>"&" ;

<*>. {
    sprintf(errmsg, "unknown character: \"%s\", skipped", yytext);
    fprintf(stderr, "file \"%s\", line %d:\n%s\n%s\n", glps_current_latfile,
            lineno, linebuf, errmsg);
    /*
      GLPS_LAT_ERROR(errmsg, glps_current_latfile,
      lineno, linebuf, GLPS_ECHAR);
    */
    /* fprintf(stderr, "\nUnknown:|%c|, Skip\n", yytext[0]);  */
    return yytext[0];
   }

%%

