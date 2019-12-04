#ifndef _MYNODE_H
#define _MYNODE_H

#include <lemon/list_graph.h>

struct node {
    lemon::ListDigraph::Node data;
    node *next; // pointer called "next" points to node
};

#endif