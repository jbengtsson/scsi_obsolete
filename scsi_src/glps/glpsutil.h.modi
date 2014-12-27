// glpsutil.h
//
// Copyright (C) 2011 Lingyun Yang (lyyang at bnl dot gov)
// 

#ifndef GLPSUTIL_H_
#define GLPSUTIL_H_

#include <string>
#include <list>
// #include <vector>

//using std::list;
//using std::string;
//using std::vector;
using namespace std;

struct GlpsStatement {
    string name;
    int type;
};

list<GlpsStatement> glps_statements();

list<string> glps_lines();

list<string> glps_elements();
list<string> glps_elements(
    const string& line, bool fullname = false, bool revtag=true);

string glps_corename(const string& fullname);

list<string>
glps_properties(const string& name);

string glps_trimindex(const string& name);

int glps_read(const string& param, double& val);
int glps_read(const string& param, string& val);
int glps_read(const string& param, std::vector<double>& val);

int glps_read(const string& stmt, const string& prop, double &val);
int glps_read(const string& stmt, const string& prop, string& val);
int glps_read(const string& stmt, const string& prop, std::vector<double>& val);

#endif
