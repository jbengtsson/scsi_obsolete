// glpsutil.cc
//
// Copyright (C) 2011 Lingyun Yang (lyyang at bnl dot gov)
//

#include <cstring>
#include <iostream>
#include "glpsutil.h"
#include "elemtree.h"
#include "glpserror.h"
#include "glps.h"

using namespace std;

extern GLPS_SYMB* symbtab;


extern "C" GLPS_SYMB* symb_lookup_stmt(GLPS_SYMB* head, const char* name);

std::list<GlpsStatement> glps_statements()
{
    std::list<GlpsStatement> stmtlst;
    GlpsStatement node;
    GLPS_SYMB* p = symbtab;
    while (p) {
        node.name = p->name;
        node.type = p->nodetype;
        stmtlst.push_back(node);
    }
    return stmtlst;
}

std::list<string> glps_lines()
{
    std::list<string> lines;
    GLPS_SYMB* p = symbtab;
    while (p) {
        if (p->nodetype == GLPS_LINE) {
            lines.push_back(p->name);
        }
        p = p->next;
    }
    return lines;
}

std::list<string> glps_elements()
{
    std::list<string> elem;
    GLPS_SYMB* p = symbtab;
    while (p) {
        if (p->nodetype == GLPS_ELEMENT) {
            elem.push_back(p->name);
        }
        p = p->next;
    }
    return elem;
}

std::list<string> glps_elements(const string& line, bool fullname, bool revtag)
{
    GLPS_SYMB* p = symb_lookup(symbtab, line.c_str());
    std::list<string> ret;
    if (p) {
        int n = 1024;
        char *str = (char*)calloc(n, sizeof(char));
        et_list_elements(&str, &n, p->line, NULL);

        char *tok = strtok(str, " ,");
        while (tok) {
            string s(tok);
            size_t pos = s.rfind(".");
            if (fullname) {
                ret.push_back(tok);
            } else if (pos == string::npos) {
                ret.push_back(tok);
            } else {
                size_t i = 0;
                for(i = 0; i < s.length(); ++i)
                    if (isalpha(s[i])) break;
                string tag = s.substr(0, i);
                if (revtag) ret.push_back(tag + s.substr(pos+1));
                else ret.push_back(s.substr(pos+1));
            }

            tok = strtok(NULL, " ,");
        }
        free(str);
    }
    return ret;
}

/** \brief return the core name for an element in beamline.
 *
 * skip the 'reverse' sign.
 */
string glps_corename(const string& fullname)
{
    string core = glps_trimindex(fullname);
    size_t pos = core.rfind(".");
    size_t i = 0;
    for (i = 0; i < core.length(); ++i) {
        if (isalpha(core[i])) break;
    }

    if (pos == string::npos) return core.substr(i);
    else {
        return core.substr(i, pos);
    }
}

string glps_trimindex(const string& name)
{
    size_t pos = name.rfind(".");
    if (pos == string::npos) return name;
    else {
        string tail = name.substr(pos+1);
        return name.substr(0, pos);
    }
}

std::list<string> glps_properties(const std::string& name)
{
    GLPS_SYMB* p = symb_lookup(symbtab, name.c_str());
    std::list<string> prpt;
    
    if (p && p->property) {
        p = p->property;
        if(p->name) prpt.push_back(p->name);
        while(p->next) {
            p = p->next;
            prpt.push_back(p->name);
        }
    }

    return prpt;
}

int glps_read(const string& param, double& val)
{
    string name = glps_trimindex(param);
    int idx = -1;
    if (name.length() < param.length()) {
        string tmp = param.substr(name.length() + 1);
        idx = atoi(tmp.c_str());
    }

    GLPS_SYMB* p = symb_lookup(symbtab, param.c_str());

    if (!p) return GLPS_ENOVAR;
    else if (p->nodetype == GLPS_NODE_REAL) {
        val = p->value;
        return GLPS_SUCCESS;
    } else if (p->nodetype == GLPS_NODE_VEC) {
        if (idx >= 0 && idx < p->size) {
            val = p->vec[idx];
            return GLPS_SUCCESS;
        } else if (idx >= p->size) {
            return GLPS_EBOUNDARY;
        } else {
            return GLPS_ENOTREAL;
        }
        return GLPS_ENOVAR;
   }

    return GLPS_ENOVAR;
}

int glps_read(const string& param, string& val)
{
    GLPS_SYMB* p = symb_lookup(symbtab, param.c_str());

    if (!p) {
        return GLPS_ENOVAR;
    } else {
        switch( p->nodetype ) {
        case GLPS_NODE_REAL:
        case GLPS_NODE_VEC:
            return GLPS_ENOTSTR;
        case GLPS_NODE_STR:
        case GLPS_ELEMENT:
        case GLPS_ACTION:
            val = p->str;
            return GLPS_SUCCESS;
        default:
            return GLPS_ENOVAR;
        }
    }
}

int glps_read(const string& param, vector<double>& val)
{
    GLPS_SYMB* p = symb_lookup(symbtab, param.c_str());

    if (!p) {
        return GLPS_ENOVAR;
    } else if (p->nodetype == GLPS_NODE_VEC) {
        for (size_t i = 0; i < p->size; ++i)
            val.push_back(p->vec[i]);
    
        return GLPS_SUCCESS;
    } else if (p->nodetype == GLPS_NODE_STR) {
        return GLPS_ENOTVEC;
    } else if (p->nodetype == GLPS_NODE_REAL) {
        return GLPS_ENOTVEC;
    }
}


int glps_read(const string& stmt, const string& prop, double &val)
{
    string name = glps_corename(stmt);
    GLPS_SYMB* p = symb_lookup_stmt(symbtab, name.c_str());
    if (!p) return GLPS_ENOVAR;

    GLPS_SYMB* q = symb_lookup_stmt(p->property, prop.c_str());
    if (!q) return GLPS_ENOVAR;

    if (q->nodetype == GLPS_NODE_REAL) {
        val = q->value;
        return GLPS_SUCCESS;
    } else {
        return GLPS_ENOTREAL;
    }
}

int glps_read(const string& stmt, const string& prop, string& val)
{
    string name = glps_corename(stmt);
    GLPS_SYMB* p = symb_lookup_stmt(symbtab, name.c_str());
    if (!p) return GLPS_ENOVAR;

    GLPS_SYMB* q = symb_lookup_stmt(p->property, prop.c_str());
    if (!q) return GLPS_ENOVAR;

    if (q->nodetype == GLPS_NODE_STR) {
        val = q->str;
        return GLPS_SUCCESS;
    } else {
        return GLPS_ENOTREAL;
    }

    return GLPS_ENOVAR;
}

int glps_read(const string& stmt, const string& prop, vector<double>& val)
{
    string name = glps_corename(stmt);
    GLPS_SYMB* p = symb_lookup_stmt(symbtab, name.c_str());
    if (!p) return GLPS_ENOVAR;

    GLPS_SYMB* q = symb_lookup_stmt(p->property, prop.c_str());
    if (!q) return GLPS_ENOVAR;

    if (q->nodetype == GLPS_NODE_VEC) {
        val.assign(q->size, 0.0);
        for (size_t i = 0; i < q->size; ++i) 
            val[i] = q->vec[i];
        return GLPS_SUCCESS;
    } else {
        return GLPS_ENOTREAL;
    }

    return GLPS_ENOVAR;
}
