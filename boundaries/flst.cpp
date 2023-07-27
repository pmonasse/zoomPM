/**
 * SPDX-License-Identifier: MPL-2.0+
 * @file flst.cpp
 * @brief Top-down extration of the tree of shapes
 * @author Pascal Monasse <monasse@imagine.enpc.fr>
 *
 * Copyright (c) 2018 Pascal Monasse
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the Mozilla Public License as
 * published by the Mozilla Foundation, either version 2.0 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * Mozilla Public License for more details.
 *
 * You should have received a copy of the Mozilla Public License
 * along with this program. If not, see <https://www.mozilla.org/en-US/MPL/2.0/>
 */

#include "tree.h"
#include <cassert>

struct cimage {
    int nrow, ncol;
    const unsigned char* gray;
};
typedef cimage* Cimage;
inline unsigned char gray(Cimage im, LsPoint pt)
{ return im->gray[pt.y*im->ncol+pt.x]; }

/// Strict comparison between numbers
#define COMPARE(t,a,b) (t==LsShape::INF? (a<b): (a>b))

/// Return connectivity for level set of type \a t.
inline int connectivity(LsShape::Type t) {
    return ((t == LsShape::INF)? 4: 8);
}

/// Direction of an edgel.
typedef unsigned char DirEdgel;
static const DirEdgel EAST  = 0;
static const DirEdgel NORTH = 1;
static const DirEdgel WEST  = 2;
static const DirEdgel SOUTH = 3;
static const DirEdgel DIAGONAL = 4;
static const DirEdgel NE = 4;
static const DirEdgel NW = 5;
static const DirEdgel SW = 6;
static const DirEdgel SE = 7;

/// Make a 180 turn compared to direction \a dir.
static DirEdgel turn_180(DirEdgel dir) {
    if(dir >= DIAGONAL) {
        dir -= 2;
        if(dir < DIAGONAL)
            dir += DIAGONAL;
    } else if(dir < 2)
        dir += 2;
    else
        dir -= 2;
    return dir;
}

/// Edgel, vertical or horizontal boundary between adjacent pixels.
class Edgel {
public:
    Edgel(short int x, short int y, DirEdgel d);

    bool operator==(const Edgel& e) const
    { return (pt.x == e.pt.x && pt.y == e.pt.y && dir == e.dir); }
    bool operator!=(const Edgel& e) const
    { return ! (e == *this); }

    bool inverse(Cimage im);
    LsPoint origin() const;
    bool exterior(LsPoint& ext, Cimage im) const;
    bool go_straight(Cimage im);
    void next(Cimage im, LsShape::Type type, int level);

    LsPoint pt; ///< Interior pixel coordinates (left of edgel direction)
    DirEdgel dir; ///< Direction of edgel
private:
    void turn_left(int connect);
    void turn_right(int connect);
    void finish_turn(Cimage im, int connect);
};

/// Constructor.
Edgel::Edgel(short int x, short int y, DirEdgel d)
: pt(), dir(d) {
    pt.x = x;
    pt.y = y;
}

/// Change to inverse edgel
bool Edgel::inverse(Cimage im) {
    if(! exterior(pt, im))
        return false;
    dir = turn_180(dir);
    return true;
}

/// Exterior pixel of edgel.
/// Return \a false if we are an image boundary edgel.
bool Edgel::exterior(LsPoint& ext, Cimage im) const {
    ext = pt;
    switch(dir) {
    case EAST:  return (++ext.y < im->nrow);
    case NORTH: return (++ext.x < im->ncol);
    case WEST:  return (--ext.y >= 0);
    case SOUTH: return (--ext.x >= 0);
    case NE: return (++ext.y < im->nrow && ++ext.x < im->ncol);
    case NW: return (++ext.x < im->ncol && --ext.y >= 0);
    case SW: return (--ext.y >= 0 && --ext.x >= 0);
    case SE: return (--ext.x >= 0 && ++ext.y < im->nrow);
    default: assert(false);
    }
    return false;
}

/// Go straight along current direction.
/// Return \c false if we end up outside the image.
bool Edgel::go_straight(Cimage im) {
    switch(dir) {
    case EAST:  return (++pt.x < im->ncol);
    case NORTH: return (--pt.y >= 0);
    case WEST:  return (--pt.x >= 0);
    case SOUTH: return (++pt.y < im->nrow);
    default: assert(false);
    }
    return false;
}

/// Begin a left turn.
void Edgel::turn_left(int connect) {
    if(connect == 8)
        dir += DIAGONAL;
    else if(++dir == DIAGONAL)
        dir = 0;
}

/// Begin a right turn.
void Edgel::turn_right(int connect) {
    if(connect == 8) {
        if(dir == 0)
            dir = DIAGONAL;
        --dir;
    } else {
        dir += DIAGONAL-1;
        if(dir < DIAGONAL)
            dir += DIAGONAL;
    }
}

/// Finish a left or right turn.
inline void Edgel::finish_turn(Cimage im, int connect) {
    dir -= DIAGONAL;
    if(connect == 4)
        go_straight(im);
    else if(++dir == DIAGONAL)
        dir = 0;
}

/// Move to next edgel along the level line.
void Edgel::next(Cimage im, LsShape::Type type, int level) {
    int connect = connectivity(type);
    if(dir >= DIAGONAL) {
        finish_turn(im, connect);
        return;
    }
    Edgel left(*this), right(*this);
    bool bLeftIn = left.go_straight(im), bRightIn = false;
    if(bLeftIn) {
        unsigned char v = gray(im, left.pt);
        bLeftIn = COMPARE(type, v, level);
        bRightIn = left.exterior(right.pt, im);
        if(bRightIn) {
            v = gray(im,right.pt);
            bRightIn = COMPARE(type, v, level);
        }
    }
    if(bLeftIn && ! bRightIn) // Go straight
        *this = left;
    else if(! bLeftIn && (! bRightIn || connect == 4))
        turn_left(connect);
    else {
        *this = (connect==4)? left: right;
        turn_right(connect);
    }
}

/// Return coordinates of origin of edgel.
inline LsPoint Edgel::origin() const {
    assert(dir<DIAGONAL);
    LsPoint p = pt;
    if(dir == EAST || dir == NORTH)
        ++p.y;
    if(dir == NORTH || dir == WEST)
        ++p.x;
    return p;
}

/// Initialize shape \a s, whose edgel \a e is on the boundary. One pixel of
/// the private area is found. \a level is the gray level of the parent.
static void init_shape(Cimage im, LsTree& tree,
                       LsShape& s, const Edgel& e, int level) {
    s.type = (gray(im,e.pt) < level)? LsShape::INF: LsShape::SUP;
    s.gray = (s.type==LsShape::INF)? 0: 255;
    s.bIgnore = false;
    s.bBoundary = false;
    s.area = 1;

    Edgel cur = e;
    do {
#ifdef BOUNDARY
        if(cur.dir < DIAGONAL)
            s.contour.push_back(cur.origin());
#endif
        int j = cur.pt.y * im->ncol + cur.pt.x;
        unsigned char v = im->gray[j];
        if(! COMPARE(s.type, v, s.gray)) {
            s.gray = v;
            s.pixels[0] = cur.pt;
        }
        assert(!tree.smallestShape[j] || tree.smallestShape[j] == s.parent);
        tree.smallestShape[j] = 0;
        cur.next(im, s.type, level);
    } while(cur != e);

    int i = s.pixels[0].y*im->ncol+s.pixels[0].x;
    tree.smallestShape[i] = &s;
}

/// Follow boundary of a child of shape \a s, starting at edgel \a e. Pixels
/// on the immediate exterior at the gray level of \a s are added to the
/// private area. The pixels on the immediate interior are marked as if they
/// were in the private area of \a s, to avoid following again the boundary.
static void find_child(Cimage im, LsTree& tree, LsShape& s, const Edgel& e) {
    LsShape::Type type = (gray(im,e.pt) < s.gray)? LsShape::INF: LsShape::SUP;

    Edgel cur = e;
    do {
        int i = cur.pt.y * im->ncol + cur.pt.x;
        assert(COMPARE(type, im->gray[i], s.gray));
        assert(tree.smallestShape[i] == 0 || tree.smallestShape[i] == &s);
        tree.smallestShape[i] = &s;
        LsPoint pt;
        if(cur.exterior(pt, im)) {
            i = pt.y * im->ncol + pt.x;
            if(tree.smallestShape[i] == 0 && im->gray[i] == s.gray) {
                s.pixels[s.area++] = pt;
                tree.smallestShape[i] = &s;
            }
        }
        cur.next(im, type, s.gray);
    } while(cur != e);
}

inline bool edge8(unsigned char vi, unsigned char ve) {
    if(vi == ve)
        return false;
    return (connectivity((vi<ve)? LsShape::INF: LsShape::SUP)==8);
}

/// Consider the exterior pixel of edgel \a e. If it is at the level of
/// shape \a s, add it to the private area. Otherwise, follow the boundary of
/// the child shape, adding to the private area the pixels on its immediate
/// exterior at level of \a s.
/// Return whether the edge belongs to the shape and is on its boundary.
static bool add_neighbor(Cimage im, LsTree& tree, LsShape& s, Edgel e,
                         std::vector<Edgel>& children) {
    if(! e.inverse(im)) {
        s.bBoundary = true;
        return false;
    }
    int i = e.pt.y*im->ncol + e.pt.x;
    if(! tree.smallestShape[i]) {
        if(im->gray[i] == s.gray) {
            s.pixels[s.area++] = e.pt;
            tree.smallestShape[i] = &s;
        } else {
            children.push_back(e);
            find_child(im, tree, s, e);
        }
    }
    return edge8(s.gray, im->gray[i]);
}

/// Fill the private area of shape \a s and find its children.
/// Put in \a children one seed edgel per child.
static void find_pp_children(Cimage im, LsTree& tree, LsShape& s,
                             std::vector<Edgel>& children) {
    for(int i = 0; i < s.area; i++) {
        const LsPoint& pt = s.pixels[i];
        assert(tree.smallestShape[pt.y*im->ncol+pt.x] == &s);
        Edgel e(pt.x, pt.y, EAST);

        e.dir = EAST;   bool E = add_neighbor(im, tree, s, e, children);
        e.dir = NORTH;  bool N = add_neighbor(im, tree, s, e, children);
        e.dir = WEST;   bool W = add_neighbor(im, tree, s, e, children);
        e.dir = SOUTH;  bool S = add_neighbor(im, tree, s, e, children);

        e.dir = NE;   if(N && E) add_neighbor(im, tree, s, e, children);
        e.dir = NW;   if(N && W) add_neighbor(im, tree, s, e, children);
        e.dir = SW;   if(S && W) add_neighbor(im, tree, s, e, children);
        e.dir = SE;   if(S && E) add_neighbor(im, tree, s, e, children);
    }
}

/// Add a new child to shape \a parent.
/// Fields other than family pointers are initialized later in \c init_shape().
static LsShape* add_child(LsTree& tree, LsShape& parent) {
    LsShape* old = parent.child;
    parent.child = &tree.shapes[tree.iNbShapes++];
    parent.child->parent = &parent;
    parent.child->sibling = old;
    parent.child->child = 0;
    return parent.child;
}

/// Extract tree of shapes rooted at \a root.
/// \param im the input image.
/// \param tree the output tree, where newly extracted shapes are appended.
/// \param root the current root of the tree.
/// \param e an edgel at the boundary of \a root.
/// \param level gray level of parent.
static void create_tree(Cimage im, LsTree& tree, LsShape& root,
                        const Edgel& e, int level) {
    init_shape(im, tree, root, e, level);

    std::vector<Edgel> children;
    find_pp_children(im, tree, root, children);

    std::vector<Edgel>::const_iterator it = children.begin();
    for(; it != children.end(); ++it) {
        LsShape* child = add_child(tree, root);
        child->pixels = root.pixels + root.area;
        create_tree(im, tree, *child, *it, root.gray);
        root.area += child->area;
    }
}

/// Top-down FLST algorithm.
void LsTree::flst_td(const unsigned char* gray) {
    cimage image = {nrow, ncol, gray};
    int area = ncol * nrow;

    for(int i = area-1; i >= 0; i--)
        smallestShape[i] = 0;

    shapes[0].type = LsShape::SUP;
    shapes[0].pixels = new LsPoint[area];
    Edgel e(0, 0, SOUTH);
    create_tree(&image, *this, shapes[0], e, -1);
    assert(area == shapes[0].area);
}
