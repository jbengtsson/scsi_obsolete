#include "elemtree.h"
#include <string.h>

void et_push_front(elemtree** root, elemtree *node)
{
    // if (!(*root)) return;
    if (*root == NULL) {
        *root = node;
        return;
    }

    elemtree* elem = (elemtree*) malloc(sizeof(elemtree));
    elem->right = (*root)->right;
    if(elem->right) elem->right->parent = elem;
    elem->left = (*root)->left;
    if(elem->left) elem->left->parent = elem;

    (*root)->left = node;
    node->parent = (*root);
    (*root)->right = elem;
    elem->parent = (*root);

    if (elem->left == NULL && elem->right == NULL) {
        // root pointed to a leaf, do not touch it's name
        elem->val = (*root)->val;
        (*root)->val = NULL;
    } else {
        // root was pointing non-leaf
    }
    /*
    elem = (*root)->left;
    printf("%p: %p %p %s\n", elem, elem->left, elem->right, elem->val);
    elem = (*root)->right;
    printf("%p: %p %p %s\n", elem, elem->left, elem->right, elem->val);
    */
}

void et_push_back(elemtree** root, elemtree *node)
{
    // if (!(*root)) return;
    if (*root == NULL) {
        *root = node;
        return;
    }

    elemtree* elem = et_new_node(NULL);
    elem->right = (*root)->right;
    if (elem->right) elem->right->parent = elem;
    elem->left = (*root)->left;
    if (elem->left) elem->left->parent = elem;

    (*root)->right = node;
    node->parent = (*root);
    (*root)->left = elem;  
    elem->parent = (*root);
    if (elem->left == NULL && elem->right == NULL) {
        elem->val = (*root)->val;
        (*root)->val = NULL;
    } else {
        elem->val = NULL;
    }
}

elemtree* et_combine(elemtree* left, elemtree *right)
{
    if (!left && !right) return NULL;
    else if (!left) return right;
    else if (!right) return left;

    elemtree *elem = (elemtree*) malloc(sizeof(elemtree));
    elem->left = left;
    elem->right = right;
    elem->val = NULL;
    left->parent = elem;
    right->parent = elem;

    return elem;
}

/* copy the whole binary tree, attach it to parent */
elemtree* et_copy_tree(elemtree* root, elemtree* parent, int rev)
{
    if (root == NULL) return NULL;

    elemtree * node = (elemtree*) malloc(sizeof(elemtree));

    node->factor = root->factor*rev;
    if (root->val) node->val = strdup(root->val);
    else node->val = NULL;
    node->parent = parent;

    if (rev == -1) {
        // reverse order
        node->left = et_copy_tree(root->right, node, rev);
        node->right = et_copy_tree(root->left, node, rev);
        return node;
    } else if (rev == 1) {
        node->left = et_copy_tree(root->left, node, rev);
        node->right = et_copy_tree(root->right, node, rev);
        return node;
    } else {
        /* ERROR: */
        return NULL;
    }
}

elemtree* et_find(elemtree* root, char* element)
{
    if (root == NULL) return NULL;
    else if (root->val && strcmp(root->val, element) == 0) return root;

    elemtree *elem = et_find(root->left, element);
    if (elem) return elem;

    return et_find(root->right, element);
}

elemtree* et_new_node(const char* val)
{
    elemtree *elem = (elemtree *)malloc(sizeof(elemtree));
    elem->left = elem->right = NULL;
    if (val) {
        elem->val = (char*)malloc(strlen(val)+1);
        strcpy(elem->val, val);
    } else {
        elem->val = NULL;
    }

    elem->factor = 1;
    elem->parent = NULL;

    return elem;
}

// can not set the name before having any node
void et_set_beamline_name(elemtree* root, char* name)
{
    if (root == NULL) {
        return;
    } 

    if (root->left == NULL && root->right == NULL) {
        // it is a leaf
        elemtree *elem = et_new_node(root->val);
        free(root->val);
        root->val = strdup(name);
        root->left = elem;
    } else {
        if (root->val) free(root->val);
        root->val = strdup(name);
    }
}

void et_delete(elemtree* root)
{
    
    if (root == NULL) return;
    if (root->left) et_delete(root->left);
    if (root->right) et_delete(root->right);

    if (root->val) free(root->val);
    free(root);
}


void et_count(elemtree *root, int *leaf, int *node)
{
    if (root == NULL) return;

    ++(*node);
    if (root->left == NULL && root->right == NULL) ++(*leaf);

    et_count(root->left, leaf, node);
    et_count(root->right, leaf, node);
}

void et_print(const elemtree* root)
{
    static int indent = 0;
    int i, isline = 0;
    if (!root) return;
    if (root->left == NULL && root->right == NULL) {
        if (root->val) printf("%s ", root->val);
    } else {
        // nonleaf
        /* fprintf(stderr, "val= %p\n", root->val); */
       if (root->val != NULL) {
            isline = 1;
            /* printf("\n"); */
            for (i = 0; i < indent; ++i) printf("  ");
            printf("%s: ", root->val);
            ++indent;
        }
    }

    if (root->left) et_print(root->left);
    if (root->right) et_print(root->right);
    
    if (isline) {
        /* for (i = 0; i < indent; ++i) printf("  "); */
        printf("\n");
        --indent;
    }
}

elemtree *glroot;

// reset pointers which point to other beamline
void et_clean_tree(elemtree* root)
{
    // delete this tree, if it has a non-null val
    if (root == NULL) return;
    else if (root->left == NULL && root->right == NULL) return;
    else 
    if (root->left != NULL && root->left->val != NULL &&
        (root->left->left != NULL || root->left->right != NULL)) {
        root->left = NULL;
    } 
    if (root->right != NULL && root->right->val != NULL &&
        (root->right->left != NULL || root->right->right !=NULL)) {
        root->right = NULL;
    }

    et_clean_tree(root->left);
    et_clean_tree(root->right);
}

// assumes the beamline are declared before using, therefore we can delete them in order.
void et_reset_tree(elemtree* parent, elemtree *root, const char* elem)
{
    if (root == NULL) return;
    else if (root->val && strcmp(root->val, elem) == 0) {
        /*fprintf(stderr, "FOUND: \n");
        et_print(parent);
        et_print(root);
        */
        /*
        fprintf(stderr, "%p %p %p\n", parent, parent->left, parent->right);
        fprintf(stderr, "%p %p %p\n", root, root->left, root->right);
        */

        if (parent->left == root) parent->left = NULL;
        else if (parent->right == root) parent->right = NULL;
        else {
            /* fprintf(stderr, "What is wrong?\n"); */
        }
        /*
        fprintf(stderr, "%p %p %p\n", parent, parent->left, parent->right);
        fprintf(stderr, "%p %p %p\n", root, root->left, root->right);
        et_print(parent);
           et_print(root); */
        return;
    } else {
        et_reset_tree(root, root->left, elem);
        et_reset_tree(root, root->right, elem);
    }
}


void et_list_elements(char** str, int* maxlen, elemtree* root, elemstack* stk)
{
    /* if it is a leaf */
    char * buf;
    const int NINC = 512;
    elemstack *p = stk;
    if (root->left == NULL && root->right == NULL) {
        /* fprintf(stderr, "Stack: "); */
        /* count length */
        p = stk;
        int n = 0, n2;
        while (p) {
            n += strlen(p->line) + 1;
            p = p->next;
        }
        n += strlen(root->val) + 5;
        while ( strlen(*str) + n >= *maxlen) {
            buf = (char*) malloc((*maxlen+n+NINC)*sizeof(char));
            strncpy(buf, *str, *maxlen);
            free(*str);
            *str = buf;
            *maxlen += n+NINC;
            /* fprintf(stderr, "inc: %d, ", n+NINC); */
        }
        /* fprintf(stderr, "\n"); */
        if (root->factor == -1) {
            strcat(*str, "-");
        } else if (root->factor != 1) {
            /* ERROR */
            strcat(*str, "?????");
        }
        p = stk;
        while(p) {
            strcat(*str, p->line);
            strcat(*str, ".");
            /* fprintf(stderr, "%s, ", p->line); */
            p = p->next;
        }
        strcat(*str, root->val);
        strcat(*str, ", ");
        /* fprintf(stderr, "\nLeaf: %s\n", root->val); */
        return;
    }

    /* it is a root */
    if (root->val) {
        if(p) {
            while(p->next) p = p->next;
            p->next = (elemstack*)malloc(sizeof(elemstack));
            p = p->next;
        } else {
            p = (elemstack*)malloc(sizeof(elemstack));
            stk = p;
        }
        p->line = strdup(root->val);
        p->next = NULL;
    }

    if (root->left) et_list_elements(str, maxlen, root->left, stk);
    if (root->right) et_list_elements(str, maxlen, root->right, stk);

    // finished both children
    // Printing
    /* fprintf(stderr, "Leaving: %s\n", root->val); */
    /* fprintf(stderr, "stack: "); */
    /* p = stk; */
    /* while(p) { */
    /*     fprintf(stderr, "%s, ", p->line); */
    /*     p = p->next; */
    /* } */
    /* fprintf(stderr, "\n"); */

    // if this is a line, remove this from stack, since we are leaving it.
    if (root->val) {
        // remove the tail element of stack
        p = stk;
        if (p->next) {
            // not the last "tag" in stack
            while (p->next->next != NULL) p = p->next;
            if (p->next->line) free(p->next->line);
            free(p->next);
            p->next = NULL;
        }else {
            // 
            /* fprintf(stderr, "Free %s from stack\n", p->line); */
            free(p->line);
            free(p);
            stk = p = NULL;
        }
    }

    // 
    /* fprintf(stderr, "clean stack: stk= %p p= %p\n    ", stk, p); */
    /* p = stk; */
    /* while(p) { */
    /*     fprintf(stderr, "%s, ", p->line); */
    /*     p = p->next; */
    /* } */
    /* fprintf(stderr, "Gone \n\n"); */
}

#ifdef DEBUG_ELEMTREE
int main()
{
    elemtree* root1 = NULL, *root2 = NULL;
    
    elemtree* elem, *pt;
    elem = et_new_node("a001");

    et_push_front(&root1, elem);
    et_set_beamline_name(root1, "ROOT1");

    et_print(root1);

    pt = et_new_node("a002");
    et_push_front(&root1, pt);

    et_print(root1);
    printf("%p\n", root1);
    pt = et_new_node("a004");
    et_push_back(&root1, pt);
    printf("%p\n", root1);
    
    et_print(root1);
    printf("\n");
    printf("%p: %p %p %s\n", root1, root1->left, root1->right, root1->val);
    elem = root1->left;
    printf("%p: %p %p %s\n", elem, elem->left, elem->right, elem->val);
    elem = root1->right;
    printf("%p: %p %p %s\n", elem, elem->left, elem->right, elem->val);

    
    elem = et_new_node("b001");
    et_push_front(&root2, elem);
    et_set_beamline_name(root2, "ROOT2");

    pt = et_new_node("b002");
    et_push_front(&root2, pt);

    printf("%p\n", root2);
    pt = et_new_node("b003");
    et_push_front(&root2, pt);

    printf("%p\n", root2);
    pt = et_new_node("b004");
    et_push_back(&root2, pt);
    printf("%p\n", root2);

    elemtree *ring = NULL;

    elem = et_new_node("mid0000");
    ring = et_combine(root1, elem);
    et_set_beamline_name(ring, "RING");
    printf("ring: %s %p %p\n", ring->val, ring->right, elem);
    et_print(ring);

    elem = et_combine(ring, root2);
    et_set_beamline_name(elem, "FULL");

    printf("full ring\n");
    et_print(elem);
    ring = elem;

    int n1 = 0, n2 = 0;
    et_count(ring, &n1, &n2);
    printf("Leaf: %d, node: %d\n", n1, n2);

    et_print(ring);
    /* printf("clean...\n");
       et_clean_tree(ring->left);
       et_clean_tree(ring->right);
    */

    int nlen = 4;
    char* str = (char*)malloc(nlen*sizeof(char));
    char* line = NULL;
    et_list_elements(&str, &nlen, ring, NULL);
    printf("LINE: %s\n", str);
    free(str);

    nlen = 4;
    str = (char*)malloc(nlen*sizeof(char));
    elemtree *rev = et_copy_tree(ring, NULL, -1);
    et_list_elements(&str, &nlen, rev, NULL);
    printf("LINE: %s\n", str);
    free(str);

    et_delete(ring);
    printf("size: %d\n", sizeof(elemtree));

    return 0;
}
#endif
