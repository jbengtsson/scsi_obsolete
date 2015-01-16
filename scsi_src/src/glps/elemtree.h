/**
 * Lingyun YANG (lingyun.yang@gmail.com)
 */

#include <stdlib.h>
#include <stdio.h>

#ifndef ELEMTREE_H_
#define ELEMTREE_H_

#ifdef __cplusplus
extern "C" {
 #endif 

struct elemtree_el {
    int factor;
    char* val;
    struct elemtree_el *left, *right;
    struct elemtree_el *parent;
};

typedef struct elemtree_el elemtree;

struct beamline_tree_el {
    elemtree* bl;
    struct beamline_tree_el *next;
};
typedef struct beamline_tree_el bltree;

struct elemstack_el {
    char* line;
    struct elemstack_el* next;
};
typedef struct elemstack_el elemstack;

void et_push_front(elemtree** root, elemtree *node);
void et_push_back(elemtree** root, elemtree *node);
elemtree* et_combine(elemtree* left, elemtree *right);
elemtree* et_copy_tree(elemtree* root, elemtree *parent, int rev);
elemtree *et_find(elemtree* root, char* element);
elemtree* et_new_node(const char* val);
void et_set_beamline_name(elemtree* root, char* name);
void et_delete(elemtree* root);
void et_print(const elemtree* root);
/* void et_delete_beamline_dependency(bltree* bl, elemtree *elem); */
void et_clean_tree(elemtree* root);
void et_reset_tree(elemtree* parent, elemtree *root, const char* elem);
void et_list_elements(char** str, int* maxlen, elemtree* root, elemstack* stk);

#ifdef __cplusplus
}
 #endif 

#endif /* ELEMTREE_H_ */
