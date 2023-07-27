// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file lltree.cpp
 * @brief Extraction of tree of level lines from an image
 * 
 * (C) 2011-2014, 2019, Pascal Monasse <pascal.monasse@enpc.fr>
 */

#include "lltree.h"
#include <algorithm>
#include <stack>
#include <cassert>

/// Constructor
LLTree::iterator::iterator(LLTree::Node* node, TreeTraversal o)
: n(node), order(o) {
    if(n && o==PostOrder)
        goBottom();
}

/// Go to left-most leaf of current node.
void LLTree::iterator::goBottom() {
    for(LLTree::Node* b=n->child; b; b=n->child)
        n=b;
}

LLTree::Node& LLTree::iterator::operator*() const {
    return *n;
}

LLTree::Node* LLTree::iterator::operator->() const {
    return n;
}

bool LLTree::iterator::operator==(const iterator& it) const {
    return (n==it.n);
}
bool LLTree::iterator::operator!=(const iterator& it) const {
    return (n!=it.n);
}

/// Increment iterator
LLTree::iterator& LLTree::iterator::operator++() {
    if(order==PreOrder) {
        LLTree::Node* next=n->child;
        if(!next)
            while((next=n->sibling) == 0)
                if((n=n->parent) == 0)
                    break;
        n=next;
    } else { // PostOrder
        LLTree::Node* next=n->sibling;
        if(next) {
            n = next;
            goBottom();
        } else
            n = n->parent;
    }
    return *this;
}

/// Build tree structure of level lines: [2]Algorithm 4.
LLTree::LLTree(const std::vector<LLTree::Node>& n)
: nodes_(n), root_(0) {
    root_ = &nodes_.front();
}

/// Destructor
LLTree::~LLTree() {
    return; // Deactivate deletion!
    for(std::vector<Node>::iterator it=nodes_.begin(); it!=nodes_.end(); ++it)
        delete it->ll;
}
